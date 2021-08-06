import os
import random
import joblib
import numpy as np
import pandas as pd
import itertools
from operator import itemgetter
from sklearn.svm import SVC
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
import seaborn as sns
from keras.utils import to_categorical
from sklearn import preprocessing
from sklearn.manifold import TSNE
from sklearn import model_selection
from sklearn.cluster import AffinityPropagation, AgglomerativeClustering, SpectralClustering, OPTICS, MeanShift
from sklearn.mixture import GaussianMixture
from sklearn.cluster import MiniBatchKMeans, KMeans, DBSCAN, Birch
from sklearn.decomposition import PCA
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.metrics import confusion_matrix, f1_score, precision_score, recall_score, accuracy_score, auc, roc_curve


# TODO: convert this codes to class base code and find the relations between the functions


def amino_used(df):
    amino_dict = {}
    amino_appear = []
    amino_list = list('ACDEFGHIKLMNPQRSTVWY')
    for amino in amino_list:
        amino_dict[amino] = df.Sequence.str.count(amino).sum()
    for key, value in amino_dict.items():
        if value != 0:
            amino_appear.append(key)
    amino_appear.sort()
    return amino_appear


def amino_occurrence(df, amino):
    dict_count = {}
    for i in amino:
        dict_count[i] = []
        dict_count[i].append(df.Sequence.str.count(i).sum())
    df_counted = pd.DataFrame.from_dict(dict_count)
    df_counted.set_index('A', inplace=True)
    return df_counted


def random_sequence(df_m, amino_list, len_pep, sequence_number, file_name):
    cnv_list = list(df_m['Sequence'])
    amino_random = {'Sequence': []}
    bre_point = 0
    while bre_point < sequence_number:
        rnd_seq = ''.join(random.choices(amino_list, k=len_pep))
        if rnd_seq in cnv_list:
            pass
        else:
            amino_random['Sequence'].append(rnd_seq)
            bre_point = bre_point + 1
    df_r = pd.DataFrame.from_dict(amino_random)
    df_r['class'] = 0
    df_m['class'] = 1
    frames = [df_m.iloc[:sequence_number], df_r]
    df_rt = pd.concat(frames)
    df_rt.to_csv(file_name, index=False)
    df_rt.reset_index(inplace=True)
    df_rt.drop(['index'], axis=1, inplace=True)
    return df_rt


def random_sequence_generator(amino_list, len_pep, sequence_number, file_name):
    amino_random = {'Sequence': []}
    bre_point = 0
    while bre_point < sequence_number:
        rnd_seq = ''.join(random.choices(amino_list, k=len_pep))
        amino_random['Sequence'].append(rnd_seq)
        bre_point = bre_point + 1
    df = pd.DataFrame.from_dict(amino_random)
    df['class'] = [1] * (len(df) - int(len(df) / 2)) + [0] * int(len(df) / 2)
    df.reset_index(inplace=True)
    df.drop(['index'], axis=1, inplace=True)
    df.to_csv(file_name, index=False)
    return df


def random_select(df, subset='class', random_number=200):
    df_random = pd.concat([df[df[subset] == 1].sample(n=random_number), df[df[subset] == 0].sample(n=random_number)])
    df_random.reset_index(inplace=True, drop=True)
    df_random = df_random.reindex(np.random.permutation(df_random.index))
    df_random.reset_index(inplace=True, drop=True)
    return df_random


def encoding(df, amino_list):
    data = np.stack(df.loc[:, 'Sequence'].astype('str').apply(list).values, axis=0)
    e = LabelEncoder()
    e.fit(amino_list)
    encoded = e.transform(data.reshape(-1)).reshape(-1, len(data[0]))
    return encoded


def class_o_one(df):
    df['class'] = [1] * (len(df) - int(len(df) / 2)) + [0] * int(len(df) / 2)
    df.reset_index(inplace=True)
    df.drop(['index'], axis=1, inplace=True)
    df.to_csv('class01.csv', index=False)
    return df


# TODO: pca and tsne are needed to include the sequences


def pca_pep(df, encoded, n=2):
    pca = PCA(n_components=n)
    principalComponents = pca.fit_transform(encoded)
    df['pca_1'] = principalComponents[:, 0]
    df['pca_2'] = principalComponents[:, 1]
    return principalComponents, df


def t_sne_pep(df, encoded, n=2):
    t_sne = TSNE(n_components=n)
    t_distributed = t_sne.fit_transform(encoded)
    df['t_SNE_1'] = t_distributed[:, 0]
    df['t_SNE_2'] = t_distributed[:, 1]
    return t_distributed, df


def repeating_drop(df):
    position = []
    for i in range(len(df['Sequence'][0])):
        if len((df['Sequence'].astype(str).str[i]).unique()) != 1:
            position.append(i)
    df['Sequence'] = df['Sequence'].apply(lambda x: ''.join(itemgetter(*position)(x)))
    return df


def pattern_detection(df, save_name):
    ser1 = df.loc[:, 'Sequence']
    dot_list = ['.', '.']
    amino_list = list('ACDEFGHIKLMNPQRSTVWY')
    final = dot_list + amino_list
    keywords = [''.join(i) for i in itertools.product(final, repeat=4)]
    keywords_1 = [x for x in keywords if x.count('.') in [1, 2, 3]]
    res = []
    for i in keywords_1:
        if i not in res:
            res.append(i)
    final_dict = {'pattern': [], 'percentage': []}
    for i in res:
        fff = ser1.str.count(i).sum()
        if fff != 0:
            final_dict['pattern'].append(i)
            final_dict['percentage'].append(fff / len(df) * 100)
    df_final = pd.DataFrame.from_dict(final_dict)
    df_final.sort_values(by='percentage', ascending=False, inplace=True)
    df_final.to_csv(save_name, index=False)
    one_amino = []
    two_amino = []
    three_amino = []
    c1, c2, c3 = 0, 0, 0
    for i in df_final['pattern']:
        if c1 < 2:
            if i.count('.') == 1:
                one_amino.append(i)
                c1 += 1
        if c2 < 2:
            if i.count('.') == 2:
                two_amino.append(i)
                c2 += 1
        if c3 < 2:
            if i.count('.') == 3:
                three_amino.append(i)
                c3 += 1

    return one_amino, two_amino, three_amino


def locate_sequence(df, pattern, peptide_low_kd, peptide_high_kd):
    df_tmp = df.copy()
    for low_kd in peptide_low_kd:
        df_tmp.loc[(df_tmp['Sequence'].str.count(low_kd) > 0) & (df_tmp["class"] == 1), "class"] = 3

    for high_kd in peptide_high_kd:
        df_tmp.loc[(df_tmp['Sequence'].str.count(high_kd) > 0) & (df_tmp["class"] == 1), "class"] = 4

    df_tmp.loc[(df_tmp['Sequence'].str.count(pattern) > 0) & (df_tmp["class"] == 1), "class"] = 2

    return df_tmp


def sav_fig(df, label_list, title, x_label, y_label, save_name, size_list, color_column):
    plt.rcParams.update({'font.size': 65})
    fig, ax = plt.subplots(figsize=(30, 30))
    color_list = ('yellow', 'green', 'red', 'blue', 'brown', 'maroon', 'plum', 'cyan', 'teal', 'navy', 'black', 'orange'
                  , 'lime', 'purple', 'magenta', 'pink', 'olive', 'pink', 'darkblue', 'gold', 'fuchsia')
    # label_list = ('random', 'Keap1', '..G.', 'Kd sequence')
    enum_color = enumerate(color_list[:len(label_list)])
    colors = dict((i, j) for i, j in enum_color)
    labels_dict = dict(zip(label_list, color_list[:len(label_list)]))
    enum_size = enumerate(size_list[:len(label_list)])
    sizes = dict((i, j) for i, j in enum_size)
    ax.scatter(df[x_label], df[y_label],
               edgecolor='none', alpha=0.9,
               c=df[color_column].map(colors), s=df[color_column].map(sizes)
               )
    handles = [Line2D([0], [0], marker='o', color='w', markerfacecolor=v, label=k, markersize=50) for k, v in
               labels_dict.items()]
    ax.legend(title='', handles=handles, bbox_to_anchor=(1.05, 1), loc='upper left', frameon=False, prop={'size': 70})
    fig.suptitle(title + '\n', fontsize=60)
    ax.set_xlabel('\n' + x_label)
    ax.set_ylabel(y_label + '\n')
    plt.grid()
    plt.savefig(save_name, bbox_inches='tight')
    plt.figure().clear()
    plt.close()
    plt.cla()
    plt.clf()


def predict_sequence(saved_model, dataset, amino_acid, saving_path, sequence_number=200):
    df_initial = pd.read_csv(dataset)
    df = df_initial.iloc[sequence_number:]
    X = df['Sequence'].values
    X_ = np.array([list(_) for _ in X])
    le = preprocessing.LabelEncoder().fit(amino_acid)
    X_ = le.transform(X_.reshape(-1)).reshape(-1, len(df['Sequence'][sequence_number]))
    encoded = to_categorical(X_)
    X_predict = encoded.reshape(len(df['Sequence']), len(df['Sequence'][sequence_number]) * len(amino_acid))
    model = joblib.load(saved_model)

    dict_predict = {'test_Num': [], 'class': [], 'P_cl0': [], 'P_cl1': []}
    for test_sequences in range(len(X_predict)):
        y_pred = model.predict_proba(X_predict[test_sequences].reshape(1, -1))
        y_class = model.predict(X_predict[test_sequences].reshape(1, -1))
        dict_predict['test_Num'].append(test_sequences)
        dict_predict['P_cl0'].append(y_pred[0][0])
        dict_predict['P_cl1'].append(y_pred[0][1])
        dict_predict['class'].append(int(y_class))
    df_predict = pd.DataFrame.from_dict(dict_predict)
    df_predict.to_csv(saving_path + '/predicted_result.csv', index=False)
    return len(df_predict.loc[(df_predict['class'] == 1)]) / len(df_predict['class']) * 100


def pattern_evaluation(df, position, amino_acid, number_of_sequences='all'):
    sequence_occurrence = 0
    if number_of_sequences != 'all':
        df = df.head(number_of_sequences)
    for i in df['Sequence']:
        res = i[position - 1]
        if res == amino_acid:
            sequence_occurrence += 1
    return sequence_occurrence / len(df['Sequence']) * 100


def clustering(df, save_name):
    X = df[['t_SNE_1', 't_SNE_2']].to_numpy()
    df_1 = df.loc[df['class'] == 1]
    X_1 = df_1[['t_SNE_1', 't_SNE_2']].to_numpy()

    model = AffinityPropagation(damping=0.9)
    model.fit(X_1)
    yhat_1 = model.predict(X_1)
    cluster_count_propagate = len(np.unique(yhat_1))
    print(cluster_count_propagate)
    model_names = ["Birch"]
    classifiers = [Birch(threshold=0.01, n_clusters=cluster_count_propagate)]
    diff_models = zip(model_names, classifiers)

    for name, model in diff_models:
        if name == "AgglomerativeClustering" or name == "SpectralClustering":
            yhat = model.fit_predict(X)
        else:
            model.fit(X)
            yhat = model.predict(X)
        df['clusters'] = yhat

        dict_difference = {'cluster_number': [], 'c_class1': [], 'c_class0': [], 'difference': [], 'class1_ratio': []}
        for i in np.unique(yhat):
            dict_difference['cluster_number'].append(i)
            n1 = len(df.loc[(df['class'] == 1) & (df['clusters'] == i)])
            dict_difference['c_class1'].append(n1)
            n0 = len(df.loc[(df['class'] == 0) & (df['clusters'] == i)])
            dict_difference['c_class0'].append(n0)
            dict_difference['difference'].append(n1 - n0)
            dict_difference['class1_ratio'].append((n1 / (n1 + n0)) * 100)
            df_result_clustering = pd.DataFrame.from_dict(dict_difference)
            df_result_clustering.to_csv(save_name, index=False)

        df_reference = df_result_clustering.loc[(df_result_clustering['class1_ratio'] > 60)]
        cluster_list = list(df_reference['cluster_number'])
        df_select_sequence = df[df['clusters'].isin(cluster_list) & df['class'].isin([1])]
        sequence_list = list(df_select_sequence['Sequence'])

        return sequence_list, df, df_select_sequence


def svm_pep(df, amino_list, save_path, file_name):
    df = df.reindex(np.random.permutation(df.index))
    X = df['Sequence'].values
    X_ = np.array([list(_) for _ in X])
    le = preprocessing.LabelEncoder().fit(amino_list)
    X_ = le.transform(X_.reshape(-1)).reshape(-1, len(df['Sequence'][0]))
    encoded = to_categorical(X_)
    X__ = encoded.reshape(len(df['Sequence']), len(df['Sequence'][0]) * len(amino_list))
    y__ = df['class'].values
    model_names = ["SVM_Linear", "SVM_RBF", "SVM_Sigmoid"]
    classifiers = [
        SVC(kernel='linear', probability=True),
        SVC(kernel='rbf', probability=True),
        SVC(kernel='sigmoid', probability=True),
    ]

    diff_models = zip(model_names, classifiers)
    for name, model in diff_models:
        model_save_path = save_path + name
        if not os.path.exists(model_save_path):
            os.makedirs(model_save_path)
        for dr in range(1, 3):
            dict_para = {'model_name': [], 'random_state': [], 'c_cl0_train': [], 'c_cl1_train': [],
                         'c_cl0_test': [], 'c_cl1_test': [], 'accuracy': [], 'precision_0': [], 'precision_1': [],
                         'recall_0': [], 'recall_1': [], 'f1_0': [], 'f1_1': [], 'AUC_0': [], 'AUC_1': [], 'TN': [],
                         'FP': [],
                         'FN': [], 'TP': []}
            random_sequence_path = model_save_path + '/' + file_name
            if not os.path.exists(random_sequence_path):
                os.makedirs(random_sequence_path)
            random_state_path = random_sequence_path + '/random_state_' + str(dr).zfill(3)
            os.makedirs(random_state_path)
            # preparing the data
            X_train, X_test, y_train, y_test = model_selection.train_test_split(X__, y__
                                                                                , stratify=y__
                                                                                , random_state=dr
                                                                                )
            with open(random_state_path + '/X_train_' + '.npy', 'wb') as f:
                np.save(f, X_train)
            with open(random_state_path + '/X_test_' + '.npy', 'wb') as f:
                np.save(f, X_test)
            with open(random_state_path + '/y_train_' + '.npy', 'wb') as f:
                np.save(f, y_train)
            with open(random_state_path + '/y_test_' + '.npy', 'wb') as f:
                np.save(f, y_test)
            count_train_class_0 = 0
            count_train_class_1 = 0
            count_test_class_0 = 0
            count_test_class_1 = 0

            for ctr in range(len(y_train)):
                if int(y_train[ctr]) == 0:
                    count_train_class_0 = count_train_class_0 + 1
                else:
                    count_train_class_1 = count_train_class_1 + 1

            for cte in range(len(y_test)):
                if int(y_test[cte]) == 0:
                    count_test_class_0 = count_test_class_0 + 1
                else:
                    count_test_class_1 = count_test_class_1 + 1

            model.fit(X_train, y_train)
            dict_test = {'test_Num': [], 'Pred_R': [], 'Exp_R': [], 'PoF': [], 'P_cl0': [], 'P_cl1': []}
            for test_sequences in range(len(X_test)):
                y_pred = model.predict_proba(X_test[test_sequences].reshape(1, -1))
                y_class = model.predict(X_test[test_sequences].reshape(1, -1))
                dict_test['test_Num'].append(test_sequences)
                dict_test['P_cl0'].append(y_pred[0][0])
                dict_test['P_cl1'].append(y_pred[0][1])
                dict_test['Pred_R'].append(int(y_class))
                if int(y_test[test_sequences]) == 0:
                    dict_test['Exp_R'].append(0)
                elif int(y_test[test_sequences]) == 1:
                    dict_test['Exp_R'].append(1)
                if dict_test['Exp_R'][test_sequences] == dict_test['Pred_R'][test_sequences]:
                    dict_test['PoF'].append('Pass')
                else:
                    dict_test['PoF'].append('Fail')
            # save the predicted to the pandas dataframe
            df_predicted = pd.DataFrame.from_dict(dict_test)
            df_predicted.to_csv(random_state_path + '/pred_x_test.csv', index=False)
            # calculate parameters
            y_prob = model.predict_proba(X_test)
            y_test_decode = y_test
            # auc for class zero plot and csv save
            fpr0, tpr0, _ = roc_curve(y_test_decode, y_prob[:, 0], pos_label=0)
            roc_auc_0 = auc(fpr0, tpr0)
            plt.rcParams.update(plt.rcParamsDefault)
            fig = plt.figure()
            plt.plot([0, 1], [0, 1], '--')
            plt.plot(fpr0, tpr0)
            plt.xlabel("FPR-0")
            plt.ylabel("TPR-0")
            plt.title("AUC-0 : {:.2f}".format(roc_auc_0))
            fig.savefig(random_state_path + '/auc0.png')
            df_auc0 = pd.DataFrame(fpr0)
            df_auc0.columns = ["FPR-0"]
            df_auc0["TPR-0"] = tpr0.tolist()
            df_auc0.to_csv(random_state_path + '/class0.csv', index=False)
            # auc for class one plot and csv save
            fpr1, tpr1, _ = roc_curve(y_test_decode, y_prob[:, 1], pos_label=1)
            roc_auc_1 = auc(fpr1, tpr1)
            fig = plt.figure()
            plt.plot([0, 1], [0, 1], '--')
            plt.plot(fpr1, tpr1)
            plt.xlabel("FPR-1")
            plt.ylabel("TPR-1")
            plt.title("AUC-1 : {:.2f}".format(roc_auc_1))
            fig.savefig(random_state_path + '/auc1.png')
            df_auc1 = pd.DataFrame(fpr1)
            df_auc1.columns = ["FPR-1"]
            df_auc1["TPR-1"] = tpr1.tolist()
            df_auc1.to_csv(random_state_path + '/class1.csv', index=False)
            # calculate f1-score, recall
            y_predict = model.predict(X_test)
            f1_score_0 = f1_score(y_test_decode, y_predict, pos_label=0)
            f1_score_1 = f1_score(y_test_decode, y_predict, pos_label=1)
            precision_0 = precision_score(y_test_decode, y_predict, pos_label=0)
            precision_1 = precision_score(y_test_decode, y_predict, pos_label=1)
            recall_0 = recall_score(y_test_decode, y_predict, pos_label=0)
            recall_1 = recall_score(y_test_decode, y_predict, pos_label=1)
            acc_sco = accuracy_score(y_test_decode, y_predict)
            # True negative ...
            con_mat = confusion_matrix(y_test_decode, y_predict)
            tn = con_mat[0][0]
            fp = con_mat[0][1]
            fn = con_mat[1][0]
            tp = con_mat[1][1]
            fig = plt.figure()
            sns_plot = sns.heatmap(con_mat, annot=True, cmap='coolwarm')
            fig.savefig(random_state_path + '/TNFP.png')
            if f1_score_0 >= 0.1 and f1_score_1 >= 0.1:
                joblib.dump(model, random_state_path + '/trained_model' + '.sav')
            dict_para['model_name'].append(str(name))
            dict_para['random_state'].append(dr)
            dict_para['c_cl0_train'].append(count_train_class_0)
            dict_para['c_cl1_train'].append(count_train_class_1)
            dict_para['c_cl0_test'].append(count_test_class_0)
            dict_para['c_cl1_test'].append(count_test_class_1)
            dict_para['accuracy'].append(acc_sco)
            dict_para['precision_0'].append(precision_0)
            dict_para['precision_1'].append(precision_1)
            dict_para['recall_0'].append(recall_0)
            dict_para['recall_1'].append(recall_1)
            dict_para['f1_0'].append(f1_score_0)
            dict_para['f1_1'].append(f1_score_1)
            dict_para['AUC_0'].append(roc_auc_0)
            dict_para['AUC_1'].append(roc_auc_1)
            dict_para['TN'].append(tn)
            dict_para['FP'].append(fp)
            dict_para['FN'].append(fn)
            dict_para['TP'].append(tp)
            plt.close('all')
            # predicting the last sequences
            df_para = pd.DataFrame.from_dict(dict_para)
            df_para.to_csv(random_state_path + '/final_result.csv', index=False)


def gather_data(result_path, save_path):
    keap_csv_save = [str(path) for path in Path(result_path).rglob('final_result.csv')]
    appended_data = []
    random_seq = []
    for i in keap_csv_save:
        a = i.split('/')[:-1]
        df = pd.read_csv(i)
        df.drop(['c_cl0_train', 'c_cl1_train', 'c_cl0_test', 'c_cl1_test'], axis=1, inplace=True)
        appended_data.append(df)
        random_seq.append(str((a[2].split('_'))[-1]))
    appended_data = pd.concat(appended_data)

    appended_data.insert(loc=2, column='random_seq', value=random_seq)
    appended_data.sort_values(['model_name', 'random_seq'], ascending=[True, True], inplace=True)
    appended_data.to_csv(save_path, index=False)
