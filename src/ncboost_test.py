import polars as pl
import os, sys

inF = sys.argv[1]

name, extension = os.path.splitext(inF)
annotated_path = f'{name}_annotated.tsv'
output_path = f'{name}_scored.tsv'

# Load and annotate variants

variants = pl.read_csv(source = f'{inF}', 
                       separator = '\t',
                       null_values='NA',
                       schema_overrides={'chr':str, 'pos':int}
                       )

variants.head()

from ncboost_functions import add_ncboost_features, ncboost_score

variants = add_ncboost_features(variants, db_path='data/ncboost_v2_prescored/ncboost_v2_hg38_20260202_full')
variants = variants.drop('NCBoost')

variants.write_csv(file=annotated_path, 
                   separator="\t", 
                   include_header=True
                   )


# Score Variants

variants = pl.read_csv(source=annotated_path, 
                       separator='\t',
                       null_values='NA',
                       schema_overrides={'chr':str}
                       )
variants.head()

model_folder = 'ncboost_models'
variants = ncboost_score(variants, model_name='ncboost_models')
variants.head()

variants.write_csv(file=output_path,
                   separator = "\t", 
                   include_header = True
                   )
