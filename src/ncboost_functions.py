import polars as pl
import xgboost as xgb
import numpy as np
import pandas as pd
import subprocess
import tabix
from functools import partial
from tqdm import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import precision_recall_curve, auc, roc_curve
from multiprocessing.pool import ThreadPool as Pool


def get_feature_list() -> tuple[list]:
    """
    Return a tuple of 4 lists, respectively containing the name of the 4 feature subsets used by NCBoost: A, B, C and D.

    Returns:
        tuple(list) : tuple of 4 lists containing respectively the names of A, B, C and D sets of features  
    """
    A = ['GerpN', 'priPhCons', 'mamPhCons', 'verPhCons', 'priPhyloP', 'mamPhyloP', 'verPhyloP', 'GerpRS', 'GerpRSpval',
      'GerpS', 'ZooPriPhyloP', 'ZooVerPhyloP', 'ZooRoCC', 'ZooUCE']
    B = ['bStatistic', 'Roulette-AR',  'CDTS', 'mean_MAF', 'mean_MAF_afr', 'mean_MAF_ami', 'mean_MAF_amr', 'mean_MAF_asj', 'mean_MAF_eas',
    'mean_MAF_fin', 'mean_MAF_mid', 'mean_MAF_nfe', 'mean_MAF_sas']
    C = ['slr_dnds', 'gene_age', 'pLI', 'zscore_mis', 'zscore_syn', 'loeuf', 'GDI', 'ncRVIS',
    'ncGERP', 'RVIS_percentile', 'pcGERP']
    D = ['GC', 'CpG', 'SpliceAI']
    return(A, B, C, D)


def get_chr_list() -> list:
    """
    Return a list containing the 24 chromosomes [1.-22] X Y.

    Returns:
        (list) : list containing the 24 chromosomes  
    """
    return([f'{x}' for x in range(1,23)] + ['X', 'Y'])



def get_ncboost_header(db_path: str) -> list[str]:
    """
    Return a list containing the column names of NCBoost prescored file.
    
    Parameters:
        db_path (str) : the path of the folder containing the prescored files. Path be absolute or relative.

    Returns:
        list(str) : list containing the column names of NCBoost prescored file
    """
    path_to_file = f"{db_path}/WG_chr1.tsv.gz"
    command = f"zgrep -m 1 chr {path_to_file}"
    proc = subprocess.Popen([command], stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell = True)
    out, err = proc.communicate()
    out = out.decode('UTF-8').strip().split('\t')
    return(out)


def ncboost_query(l_chr: str, l_pos: int, l_ref: str, l_alt: str, tb: tabix.open) -> pl.DataFrame:
    """
    Return NCBoost annotations and score contained in the NCBoost prescored file for a specific chr:pos-pos:ref>alt SNV.
    
    Parameters:
        l_chr (str) : chromosome index of the variant, 1-22, X or Y
        l_pos (int) : absolute genomic position of the variant in GRCh38 genome assembly coordinates.
        l_ref (str) : reference allele at the position of the variant
        l_alt (str) : alternative allele representing the variant
        tb (tabix.open) : tabix handle of the indexed file

    Returns:
        pl.DataFrame : polars df containing the annotation and scores present in the NCBoost prescored file for the queried variant.
    """
    query = f"{l_chr}:{l_pos}-{l_pos}"
    records = tb.querys(query)
    res = [x for x in records]
    try:
        out = pl.DataFrame(res, orient='row')
        out = out.filter((pl.col('column_2') == l_ref) & (pl.col('column_3') == l_alt))
        return(out)
    except:
        return(pl.DataFrame(None))


def ncboost_query_chr(data: pl.DataFrame, l_chr: str, db_path: str) -> pl.DataFrame:
    """
    Query NCBoost annotations and score for all variants mapped to as specific chromosome in GRCh38 genomic coordinates.
    
    Parameters:
        data (pl.DataFrame) : polars DataFrame of N rows, each row describing a variant with at least the following fields: chr, pos, ref & alt
        l_chr (str) : chromosome index of the corresponding batch of variant
        db_path (str) : the path of the folder containing the prescored files. Path be absolute or relative

    Returns:
        pl.DataFrame : polars df containing the annotation and scores present in the NCBoost prescored file for all variants on the specified chromosome.
    """
    data = data.filter(pl.col('chr') == l_chr)
    ncboost_header = get_ncboost_header(db_path)
    old_dtypes = data.dtypes
    old_header = set(data.columns)
    ncboost_path = f"{db_path}/WG_chr{l_chr}.tsv.gz"
    tb = tabix.open(ncboost_path)
    out = []
    for var in data.iter_rows(named=True):
        l_df = ncboost_query(l_chr=var['chr'], l_pos=var['pos'], l_ref=var['ref'], l_alt=var['alt'], tb=tb)
        if not l_df.is_empty():
            out.append(l_df)
    if len(out) != 0:
        df = pl.concat(out)
        df.columns  = ncboost_header
        df = df.with_columns(pl.col('pos').cast(int))
        extra_columns = [x for x in old_header if x not in df.columns]
        data = data.select(['chr', 'pos', 'ref', 'alt'] + extra_columns)
        data = data.join(df, how='left', on=['chr', 'pos', 'ref', 'alt'])
        data = data.with_columns([pl.col(l_col).cast(pl.String) for l_col in data.columns])
        data = data.with_columns(pl.col('pos').cast(int))
    else:
        empty_header = ['chr', 'pos', 'ref', 'alt'] + [x for x in old_header if x not in ncboost_header] + [x for x in ncboost_header if x not in ['chr', 'pos', 'ref', 'alt']]
        data = pl.DataFrame(None, schema=empty_header, schema_overrides={x : str for x in empty_header})
        data = data.with_columns(pl.col('pos').cast(int))
    return(data)


def add_ncboost_features(data: pl.DataFrame, db_path: str) -> pl.DataFrame:
    """
    Query NCBoost annotations and score for all SNVs described in a polars DataFrame.
    
    Parameters:
        data (pl.DataFrame) : polars DataFrame of N rows, each row describing a variant with at least the following fields: chr, pos, ref & alt
        db_path (str) : the path of the folder containing the prescored files. Path be absolute or relative

    Returns:
        pl.DataFrame : polars df containing the annotation and scores present in the NCBoost prescored file for all input variants.
    """
    results = []
    for l_chr in tqdm(get_chr_list(), total=24, desc='Chromosome'):
        results.append(ncboost_query_chr(data, l_chr, db_path))
    return(pl.concat(results))


def ncboost_query_score(l_chr: str, l_pos: int, l_ref: str, l_alt: str, tb: tabix.open) -> pl.DataFrame:
    """
    Return NCBoost chromosome rank percentile contained in the NCBoost prescored file for a specific chr:pos-pos:ref>alt SNV.
    
    Parameters:
        l_chr (str) : chromosome index of the variant, 1-22, X or Y
        l_pos (int) : absolute genomic position of the variant in GRCh38 genome assembly coordinates.
        l_ref (str) : reference allele at the position of the variant
        l_alt (str) : alternative allele representing the variant
        tb (tabix.open) : tabix handle of the indexed file

    Returns:
        pl.DataFrame : polars df containing the annotation and scores present in the NCBoost prescored file for the queried variant.
    """
    query = f"{l_chr}:{l_pos}-{l_pos}"
    records = tb.querys(query)
    res = [x for x in records]
    try:
        out = pl.DataFrame(res, orient='row')
        out = out.filter((pl.col('column_2') == l_ref) & (pl.col('column_3') == l_alt))
        return(out[0,-1])
    except:
        return(None)


def ncboost_split_train_test(data: pl.DataFrame, l_partition: int) -> tuple[pl.DataFrame]:
    """
    Split polars DataFrame containing variants into training and testing sets based on the specified partition (from 1 to 10)

    Parameters:
        data (pl.DataFrame) : polars DataFrame of N rows, each row describing a variant with at least the following fields: chr, pos, ref, alt & partition
        l_partition (int) : index of the partition to use as train/test split

    Returns:
        tuple(pl.DataFrame) : tuple returning the training and testing sets, respecetively, as defined by the input partition id.
    """
    l_training_set = data.filter(pl.col('partition') != l_partition)
    l_testing_set = data.filter(pl.col('partition') == l_partition)
    return(l_training_set, l_testing_set)


def get_model_feature_importance(model: xgb.Booster, l_partition: int) -> pl.DataFrame:
    """
    Return the feature importance of the xgboost model trained on the specified partition (from 1 to 10)

    Parameters:
        model (xbg.Booster) : trained xgboost model
        l_partition (int) : index of the partition

    Returns:
        pl.DataFrame : polars df containing the importance of each feature in the specified model.
    """
    feature_importance_dict = model.get_score(fmap = '', importance_type = 'total_gain')
    total = sum(feature_importance_dict.values(), 0.0)
    feature_importance_dict = {k: np.round(v / total, 5) for k, v in feature_importance_dict.items()}
    l_feature_importance_df = pl.DataFrame(feature_importance_dict)
    l_feature_importance_df = l_feature_importance_df / l_feature_importance_df[:,:-1].sum_horizontal()
    l_feature_importance_df = l_feature_importance_df.with_columns(pl.lit(l_partition).cast(pl.Int64).alias('partition'))
    return(l_feature_importance_df)
    

def ncboost_train(data: pl.DataFrame, features: list, save_path: str) -> tuple[dict[xgb.Booster], pl.DataFrame, pl.DataFrame]:
    """
    Train 10 xgboost models on the specified data as train/test sets and the specified features. Feature importance will be computed, saved and returned. 

    Parameters:
        data (pl.DataFrame) : polars DataFrame containing annotated variants
        features : list of features to be used for training/testing
        save_path : path to the folder where models, scored data and training figures/metrics will be saved

    Returns:
        tuple(dict(xgb.Booster), pl.DataFrame, pl.DataFrame) : Return a dict containing the 10 models, the feature importance across the 10 models and the scored data
    """
    df_schema = {l_feature : float for l_feature in features}
    df_schema['partition'] = int
    feature_importance_df = pl.DataFrame(schema = df_schema)
    model_dict = {}
    annotated_data = pl.DataFrame()
    xgb.config_context(verbosity = 1)
    for l_partition in range(1,11):
        print(f'Training model {l_partition}')
        l_training_set, l_testing_set = ncboost_split_train_test(data, l_partition)
        xgb_params = {
            'max_depth' : 25, 
            'learning_rate' : 0.01, 
            'objective' : 'binary:logistic',
            'scale_pos_weight' : 0.91, # (N - n_pos) / N
            'nthread' : 10,
            'colsample_bytree' : 0.5,
            'min_child_weight' : 1/np.sqrt(l_training_set.shape[0]/len(features)),
            'gamma' : 10,
            'subsample' : 0.75,
            'seed' : 1,
            'eval_metric' : 'auc', #auc, aucpr, logloss
            'validate_parameters' : True,
            'tree_method' : 'exact'
            }
        l_train_features = l_training_set.select(features)
        l_train_labels = l_training_set['label']
        l_test_features = l_testing_set.select(features)
        l_test_labels = l_testing_set['label']

        dtrain = xgb.DMatrix(l_train_features, label = l_train_labels, feature_names = features)
        dtest = xgb.DMatrix(l_test_features, label = l_test_labels, feature_names = features)
        watchlist = [(dtest, "eval"), (dtrain, "train")]
        evals_result = {}
        model = xgb.train(
            xgb_params, 
            dtrain, 
            num_boost_round = 3000, 
            verbose_eval = False, 
            evals_result=evals_result, 
            evals = watchlist
        )
        y_pred = model.predict(dtest)
        l_testing_set = l_testing_set.with_columns(pl.lit(y_pred).alias('NCBoost'))
        annotated_data = pl.concat([annotated_data, l_testing_set])

        l_feature_importance_df = get_model_feature_importance(model, l_partition)
        feature_importance_df = pl.concat([feature_importance_df, l_feature_importance_df], how = 'diagonal')

        # Saving models
        model_dict[f'model_{l_partition}'] = model
        model.save_model(f'{save_path}/model_{l_partition}.json')
    annotated_data.write_csv(file = f'{save_path}/scored_data.tsv', 
                            separator = "\t", 
                            include_header = True,
                            null_value = 'NA'
                            )
    t_feature_importance_df = feature_importance_df.drop('partition')
    t_feature_importance_df = t_feature_importance_df.transpose()
    t_feature_importance_df = t_feature_importance_df.with_columns(features = pl.Series(feature_importance_df.columns[:-1]))
    t_feature_importance_df = t_feature_importance_df.with_columns(pl.mean_horizontal(pl.all()).alias('feature_mean'))
    t_feature_importance_df = t_feature_importance_df.sort('feature_mean', descending = True, nulls_last = True)
    t_feature_importance_df = t_feature_importance_df.drop('feature_mean')
    t_feature_importance_df.write_csv(file = f'{save_path}/feature_importance.tsv', 
                            separator = "\t", 
                            include_header = True,
                            null_value = 'NA'
                            )
    return(model_dict, t_feature_importance_df, annotated_data)


def plot_feature_importance(feature_importance_df: pl.DataFrame, save_path: str) -> None:
    """
    Plot and save the feature importance (+/- std) of the feature set used to train the 10 models 

    Parameters:
        feature_importance_df : polars DataFrame containing the feature importance of the M features across the 10 models
        save_path : path to the folder where models, scored data and training figures/metrics were saved

    Returns:
        (None)
    """
    feature_names = feature_importance_df['features'].to_list()
    feature_importance_df = feature_importance_df.unpivot(index = 'features')
    fig, ax = plt.subplots(figsize=(5,8))
    sns.barplot(feature_importance_df.to_pandas(), y = 'features', x = 'value', order = feature_names, errorbar = 'sd', orient = 'h')
    plt.tight_layout()
    plt.savefig(f'{save_path}/Feature_importance.png', bbox_inches='tight', format='png', dpi=300)
    # plt.savefig(f'{save_path}/Feature_importance.pdf', bbox_inches='tight', format='pdf')
    plt.show()


def get_fpr_tpr_auc(df: pl.DataFrame, score: str) -> tuple[list, list, float]:
    """
    Compute the fpr, tpr and auroc of the dataset for the queried score

    Parameters:
        df (pl.DataFrame) : polars DataFrame containing at least the label and the queried score for each variant
        score (str) : name of score to be considered to compute the tpr, fpr and auroc

    Returns:
        tuple(list, list, float) : returns fpr and tpr at each split and the corresponding auroc
    """
    df = df.filter(pl.col(score).is_not_null())
    fpr, tpr, _ = roc_curve(df['label'], df[score])
    auroc = auc(fpr, tpr)
    return fpr, tpr, auroc


def get_pre_rec_auc(df: pl.DataFrame, score: str) -> tuple[list, list, float]:
    """
    Compute the precision, recall and aupr of the dataset for the queried score

    Parameters:
        df (pl.DataFrame) : polars DataFrame containing at least the label and the queried score for each variant
        score (str) : name of score to be considered to compute the precision, recall and aupr

    Returns:
        tuple(list, list, float) : returns precision and recall at each split and the corresponding aurpr
    """
    df = df.filter(pl.col(score).is_not_null())
    precision, recall, _ = precision_recall_curve(df['label'], df[score])
    auprc  = auc(recall, precision)
    return precision, recall, auprc


def plot_roc_prc(annotated_data: pl.DataFrame, save_path: str, scores: list[str]= ['NCBoost', 'CADD', 'ReMM'], figure_name: str= 'ROC_PRC') -> None:
    """
    Plot the ROC and PR curves of the input data for the specified scores and save the corresponding figure

    Parameters:
        annotated_data (pl.DataFrame) : polars DataFrame containing the variants' label and scores
        save_path (str) : path to the folder where the figure will be saved
        scores (list(str)) : list of names of scores to be added to the ROC and PR curves
        figure_name (str) : name of the file containing the plot - default is 'ROC_PRC'

    Returns:
        (None)
    """    
    fig, axs = plt.subplots(1, 2, figsize=(8,4))
    lw = 2
    color_dict = {'NCBoost' : '#b5179e', 
                  'CADD' : 'orange',
                  'ReMM' : 'blue'}
    for l_score in scores:
        fpr, tpr, auroc = get_fpr_tpr_auc(annotated_data, l_score)
        axs[0].plot(fpr, tpr, color=color_dict[l_score], lw=lw, label=f'{l_score} : {auroc:0.2f}')
        precision, recall, auprc = get_pre_rec_auc(annotated_data, l_score)
        axs[1].plot(recall, precision, color=color_dict[l_score], lw=lw, label=f'{l_score} : {auprc:0.2f}')

    axs[0].plot([0, 1], [0, 1], color='navy', lw=1, linestyle='--')
    axs[0].set_xlim([-0.02, 1.0])
    axs[0].set_ylim([0.0, 1.05])
    axs[0].set_xlabel('FPR')
    axs[0].set_ylabel('TPR')
    axs[0].spines['top'].set_visible(False)
    axs[0].spines['right'].set_visible(False)
    axs[0].title.set_text('ROC curve')
    axs[0].legend(loc="lower right")

    axs[1].set_xlim([0.0, 1.0])
    axs[1].set_ylim([0.0, 1.05])
    axs[1].set_xlabel('Recall')
    axs[1].set_ylabel('Precision')
    axs[1].spines['top'].set_visible(False)
    axs[1].spines['right'].set_visible(False)
    axs[1].title.set_text('Precision-Recall curve')
    axs[1].legend(loc="upper right")

    plt.tight_layout()
    plt.savefig(f'{save_path}/{figure_name}.png', format='png', bbox_inches='tight', dpi=300)
    # plt.savefig(f'{save_path}/{figure_name}.pdf', format='pdf', bbox_inches='tight')
    plt.show()


def ncboost_score(data: pl.DataFrame, model_name: str= 'ncboost_models') -> pl.DataFrame:
    """
    Score the input data using the specific model bundle

    Parameters:
        data (pl.DataFrame) : polars DataFrame containing the annotated variants to be scored
        model_name : name of the folder containing the bundle of 10 models to be used to infer NCBoost score

    Returns:
        pl.DataFrame : polars dataframe containing variants and their corresponding NCBoost score
    """  
    geneDB_path = 'data/geneDB_ncboost2.tsv'
    partition_data = pl.read_csv(geneDB_path, 
                                    separator='\t', 
                                    has_header=True, 
                                    null_values = 'NA',
                                    schema_overrides={'chr' : str},
                                    columns=['ENSG', 'chr', 'partition']
                                    )
    data = data.join(partition_data, left_on=['closest_gene_ENSG', 'chr'], right_on = ['ENSG', 'chr'], how='left')
    data = data.drop('partition_right', strict=False)
    models = load_models(f'models/{model_name}')
    data_model_tuples = [(data.filter(pl.col('partition') == l_id), models[l_id]) for l_id in range(1,11)]
    scored_list = []
    with Pool(processes = 10) as pool:
        for l_partition in pool.map(apply_model, data_model_tuples):
            scored_list.append(l_partition)
    scored_data = pl.concat(scored_list, how = 'vertical')
    scored_data = scored_data.sort(by=['chr', 'pos'])
    return(scored_data)


def load_models(path) -> dict[xgb.Booster]:
    """
    Load the bundle of trained 10 models.

    Parameters:
        path : path to the folder containing the bundle of 10 models to be loaded

    Returns:
        (dict(xgb.Booster)) : dictionnary containing the 10 trained ncboost models
    """ 
    model_dict = {}
    for l_id in range(1,11):
        model = xgb.XGBRegressor()
        model.load_model(f'{path}/model_{l_id}.json')
        model_dict[l_id] = model
    return(model_dict)


def apply_model(data_model_tuple: tuple[pl.DataFrame, xgb.Booster]) -> pl.DataFrame:
    """
    Use one specified NCBoost model to infer the score of the input data 

    Parameters:
        data_model_tuple (tuple(pl.DataFrame, xgb.Booster)) : tuple containing a polars DataFrame with the variant features and the ncboost model of the corresponding partition.

    Returns:
        (pl.DataFrame): polars DataFrame containing the input data scored by the provided NCBoost model.
    """ 
    data, model = data_model_tuple
    y_pred = model.predict(np.asarray(data[model.feature_names_in_]))
    data = data.with_columns(pl.lit(y_pred).cast(float).alias('NCBoost'))
    return(data)
