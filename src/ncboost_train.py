# To run this script, make sur that you've downloaded and setup NCBoost annotation file.
import numpy as np
import polars as pl
import sys, os
sys.path.append('.')


# ## 1 - Load and annotate variants
# This first section load the training positive & negative variant sets.
# It requires the user to provide a set of positive and negative variants, containing at least:
# 
# chr: str, chromosome bearing the variant, [1:22], X, Y
# 
# pos: int, hg38 genomic position of the variant
# 
# ref: str, reference allele at the position
# 
# alt: str, alternative allele at the position

from src.ncboost_functions import add_ncboost_features

positives =  pl.read_csv(source='data/training/positives_for_training.tsv',
                         separator = '\t',
                         schema_overrides = {'chr' : str},
                         null_values = 'NA'
                         )
positives = add_ncboost_features(positives, db_path='data/ncboost_v2_prescored/ncboost_v2_hg38_20260202_full')
positives = positives.drop(['NCBoost', 'NCBoost_chr_rank_perc'])
positives = positives.with_columns(label=1)
print(positives.shape)
positives.head()

negatives =  pl.read_csv(source=f'data/training/negatives_for_training.tsv', 
                         separator='\t',
                         schema_overrides={'chr' : str},
                         null_values='NA'
                         )
negatives = add_ncboost_features(negatives, db_path='data/ncboost_v2_prescored/ncboost_v2_hg38_20260202_full')
negatives = negatives.drop(['NCBoost', 'NCBoost_chr_rank_perc'])
negatives = negatives.with_columns(label=0)
print(negatives.shape)
negatives.head()

positives.write_csv(file='data/training/positives_annotated.tsv',
                    separator="\t", 
                    include_header=True,
                    null_value='NA'
                    )
negatives.write_csv(file='data/training/negatives_annotated.tsv', 
                    separator="\t", 
                    include_header=True,
                    null_value='NA'
                    )

# ## 2 - Training NCBoost framework


from src.ncboost_functions import get_feature_list

A, B, C, D = get_feature_list()
variables = ['chr', 'pos', 'ref', 'alt', 'region', 'closest_gene_name', 'closest_gene_ENSG', 
             'label', 'partition', 'matching_index', 'CADD', 'ReMM']
region_list = ['upstream', 'downstream', 'intronic', 'UTR5', 'UTR3', 'intergenic']
features = A + B + C + D + region_list



from src.ncboost_functions import ncboost_train
model_name = f'toy_models'
save_path = f'models/{model_name}'
if not os.path.exists(save_path):
    os.makedirs(save_path)


positives =  pl.read_csv(source = 'data/training/positives_annotated.tsv', 
                   separator = '\t',
                   schema_overrides = {'chr' : str},
                   null_values='NA',
                   )

negatives =  pl.read_csv(source = 'data/training/negatives_annotated.tsv', 
                   separator = '\t',
                   schema_overrides = {'chr' : str},
                   null_values='NA',
                   )


data = pl.concat([negatives.select(variables + features), positives.select(variables + features)], how = 'vertical')

model_dict, feature_importance, annotated_data = ncboost_train(data, features, save_path)


annotated_data.write_csv(file = f'{save_path}/training_data_scored.tsv', 
                        separator = "\t", 
                        include_header = True,
                        null_value = 'NA'
                        )


feature_importance =  pl.read_csv(source = f'{save_path}/feature_importance.tsv', 
                   separator = '\t',
                   null_values='NA',
                   )
from src.ncboost_functions import plot_feature_importance
plot_feature_importance(feature_importance, save_path)


annotated_data =  pl.read_csv(source = f'{save_path}/training_data_scored.tsv', 
                   separator = '\t',
                   null_values='NA',
                   schema_overrides={'chr':str}
                   )
from src.ncboost_functions import plot_roc_prc
plot_roc_prc(annotated_data, save_path, scores = ['NCBoost', 'CADD', 'ReMM'])





