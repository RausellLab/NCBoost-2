
import polars as pl
import xgboost as xgb
import numpy as np
import pandas as pd

def ncboost_train(data, features, save_path):
    
    df_schema = {l_feature : float for l_feature in features}
    df_schema['partition'] = int
    feature_importance_df = pl.DataFrame(schema = df_schema)

    model_dict = {}

    annotated_data = pl.DataFrame()
    
    xgb.config_context(verbosity = 1)

    for l_partition in range(1,11):
        print(f'Training model {l_partition}')
        l_training_set = data.filter(pl.col('partition') != l_partition)
        l_testing_set = data.filter(pl.col('partition') == l_partition)
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
        model = xgb.train(xgb_params, dtrain, num_boost_round = 3000, verbose_eval = False, evals_result=evals_result, evals = watchlist)

        y_pred = model.predict(dtest)
        l_testing_set = l_testing_set.with_columns(pl.lit(y_pred).alias('NCBoost'))
        annotated_data = pl.concat([annotated_data, l_testing_set])

        feature_importance_dict = model.get_score(fmap = '', importance_type = 'total_gain')
        total = sum(feature_importance_dict.values(), 0.0)
        feature_importance_dict = {k: np.round(v / total, 5) for k, v in feature_importance_dict.items()}
        l_feature_importance_df = pl.DataFrame(feature_importance_dict)
        l_feature_importance_df = l_feature_importance_df / l_feature_importance_df[:,:-1].sum_horizontal()
        l_feature_importance_df = l_feature_importance_df.with_columns(pl.lit(l_partition).cast(pl.Int64).alias('partition'))
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

