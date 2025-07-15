import polars as pl
import xgboost as xgb
import numpy as np
import pandas as pd
import subprocess
import tabix
from functools import partial
from tqdm.notebook import tqdm
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.metrics import precision_recall_curve, auc, roc_curve
from multiprocessing.pool import ThreadPool as Pool


def get_feature_list():
    A = ['GerpN', 'priPhCons', 'mamPhCons', 'verPhCons', 'priPhyloP', 'mamPhyloP', 'verPhyloP', 'GerpRS', 'GerpRSpval',
      'GerpS', 'ZooPriPhyloP', 'ZooVerPhyloP', 'ZooRoCC', 'ZooUCE', 'Roulette-AR']
    B = ['bStatistic',  'CDTS', 'mean_MAF', 'mean_MAF_afr', 'mean_MAF_ami', 'mean_MAF_amr', 'mean_MAF_asj', 'mean_MAF_eas',
    'mean_MAF_fin', 'mean_MAF_mid', 'mean_MAF_nfe', 'mean_MAF_sas']
    C = ['slr_dnds', 'gene_age', 'pLI', 'zscore_mis', 'zscore_syn', 'loeuf', 'GDI', 'ncRVIS',
    'ncGERP', 'RVIS_percentile', 'pcGERP']
    D = ['GC', 'CpG', 'SpliceAI']
    return(A, B, C, D)

def get_ncboost_header(db_path):
    path_to_file = f"{db_path}/WG_chr1.tsv.gz"
    command = f"zgrep -m 1 chr {path_to_file}"
    proc = subprocess.Popen([command], stdin=subprocess.PIPE, stdout=subprocess.PIPE, shell = True)
    out, err = proc.communicate()
    out = out.decode('UTF-8').strip().split('\t')
    return(out)

def ncboost_query(l_chr, l_pos, tb):
    query = f"{l_chr}:{l_pos}-{l_pos}"
    records = tb.querys(query)
    res = [x for x in records]
    return(res[0])

def add_ncboost_chr(l_chr, data, db_path):
    ncboost_path = f"{db_path}/WG_chr{l_chr}.tsv.gz"
    tb = tabix.open(ncboost_path)
    data = data.filter(pl.col('chr') == l_chr)
    data = data.with_columns(pl.struct(['chr', 'pos'])
                             .map_elements(lambda x: partial(ncboost_query, tb = tb)(x['chr'], x['pos']), 
                      return_dtype = pl.List(pl.String)
                      ).alias('ncboost_features'))
    return(data)

def add_ncboost_features(data, db_path):
    old_header = data.columns
    ncb_header = get_ncboost_header(db_path)
    result_dic = {}
    data_dict = data.partition_by('chr', as_dict = True)
    for l_chr in tqdm(data_dict.keys(), desc="Chromosomes"):
        result_dic[l_chr[0]] = add_ncboost_chr(l_chr = l_chr[0], data = data_dict[l_chr], db_path = db_path)
    data = pl.concat(result_dic.values())
    data = data.with_columns(struct = pl.col('ncboost_features').list.to_struct()).unnest('struct').drop('ncboost_features')
    data = data.drop([f'field_{x}' for x in [0, 1, 2]])
    data.columns = old_header + ncb_header[3:]
    return(data)

def ncboost_split_train_test(data, l_partition):
    l_training_set = data.filter(pl.col('partition') != l_partition)
    l_testing_set = data.filter(pl.col('partition') == l_partition)
    return(l_training_set, l_testing_set)

def get_model_feature_importance(model, l_partition):
    feature_importance_dict = model.get_score(fmap = '', importance_type = 'total_gain')
    total = sum(feature_importance_dict.values(), 0.0)
    feature_importance_dict = {k: np.round(v / total, 5) for k, v in feature_importance_dict.items()}
    l_feature_importance_df = pl.DataFrame(feature_importance_dict)
    l_feature_importance_df = l_feature_importance_df / l_feature_importance_df[:,:-1].sum_horizontal()
    l_feature_importance_df = l_feature_importance_df.with_columns(pl.lit(l_partition).cast(pl.Int64).alias('partition'))
    return(l_feature_importance_df)
    
def ncboost_train(data, features, save_path):
    # data: variant dataset as input, already annotated with NCBoost features and with 
    # 'label' information:  for negatives,  for positives.
    # features: list of features present in data to be used for training.
    # save_path: the path to save the trained models and the feature importance file.  
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

def plot_feature_importance(feature_importance_df, save_path):
    feature_names = feature_importance_df['features'].to_list()
    feature_importance_df = feature_importance_df.unpivot(index = 'features')
    fig, ax = plt.subplots(figsize=(5,8))
    sns.barplot(feature_importance_df.to_pandas(), y = 'features', x = 'value', order = feature_names, errorbar = 'sd', orient = 'h')
    plt.tight_layout()
    plt.savefig(f'{save_path}/Feature_importance.png', bbox_inches='tight', format='png', dpi=300)
    plt.savefig(f'{save_path}/Feature_importance.pdf', bbox_inches='tight', format='pdf')
    plt.show()

def get_fpr_tpr_auc(df, score):
    df = df.filter(pl.col(score).is_not_null())
    fpr, tpr, _ = roc_curve(df['label'], df[score])
    auroc = auc(fpr, tpr)
    return fpr, tpr, auroc

def get_pre_rec_auc(df, score):
    df = df.filter(pl.col(score).is_not_null())
    precision, recall, _ = precision_recall_curve(df['label'], df[score])
    auprc  = auc(recall, precision)
    return precision, recall, auprc


def plot_roc_prc(annotated_data, save_path, scores = ['NCBoost', 'CADD', 'ReMM'], figure_name='ROC_PRC'):
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
    plt.savefig(f'{save_path}/{figure_name}.pdf', format='pdf', bbox_inches='tight')
    plt.show()

def ncboost_score(data, model_name='ncboost_models'):
    geneDB_path = 'data/geneDB_ncboost2.tsv'
    partition_data = pl.read_csv(geneDB_path, 
                                    separator='\t', 
                                    has_header=True, 
                                    null_values = 'NA',
                                    schema_overrides={'chr' : str},
                                    columns=['ENSG', 'chr', 'partition']
                                    )
    data = data.join(partition_data, left_on=['closest_gene_ENSG', 'chr'], right_on = ['ENSG', 'chr'], how='left')
    model_path = f'models/{model_name}'
    models = load_models(model_path)
    data_model_tuples = [(data.filter(pl.col('partition') == l_id), models[l_id]) for l_id in range(1,11)]
    scored_list = []
    with Pool(processes = 10) as pool:
        for l_partition in pool.map(apply_model, data_model_tuples):
            scored_list.append(l_partition)
    scored_data = pl.concat(scored_list, how = 'vertical')
    scored_data = scored_data.sort(by=['chr', 'pos'])
    return(scored_data)

def load_models(path):
    model_dict = {}
    for l_id in range(1,11):
        model = xgb.XGBRegressor()
        model.load_model(f'{path}/model_{l_id}.json')
        model_dict[l_id] = model
    return(model_dict)

def apply_model(data_model_tuple):
    data, model = data_model_tuple
    y_pred = model.predict(np.asarray(data[model.feature_names_in_]))
    data = data.with_columns(pl.lit(y_pred).cast(float).alias('NCBoost'))
    return(data)