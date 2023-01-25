import time
import re
import operator
import os
import sys
import glob
import random
import joblib
import numpy as np
import pandas as pd
import itertools
import plotly.express as px
from operator import itemgetter
from sklearn.svm import SVC
from sklearn.manifold import TSNE
import matplotlib.pyplot as plt
from matplotlib.lines import Line2D
from pathlib import Path
import seaborn as sns
from tensorflow.keras.utils import to_categorical
from sklearn import preprocessing
from sklearn import model_selection
from sklearn.neighbors import NearestNeighbors
from sklearn.cluster import AffinityPropagation, AgglomerativeClustering, SpectralClustering, OPTICS, MeanShift
from sklearn.mixture import GaussianMixture
from sklearn.cluster import MiniBatchKMeans, KMeans, DBSCAN, Birch
from sklearn.decomposition import PCA
from sklearn.preprocessing import LabelEncoder, OneHotEncoder
from sklearn.preprocessing import StandardScaler
from sklearn.metrics import confusion_matrix, f1_score, precision_score, recall_score, accuracy_score, auc, roc_curve

def read_csv(prefix):
    path = os.getcwd()
    #csv_files = glob.glob(os.path.join(path, "Exp6R.motifs.gap_length*"))
    csv_files = glob.glob(os.path.join(path, "{}*".format(prefix)))
    list_df = []
    for f in csv_files:
        df = pd.read_csv(f)
        list_df.append(df)
        df_c = pd.concat(list_df, ignore_index=True)
    return df_c

def search_candidates(df1,df2):
    dict1 = {}
    for i in range(len(df1['Sequence'])):
        for j in range(len(df2['motifs'])):
            x = df2.iloc[j]['motifs']
            if re.search(x,df1.iloc[i]['Sequence']):
                dict1[df2.iloc[j]['motifs']] = [df2.iloc[j]['index_pos'],df2.iloc[j]['index_neg']]
    dict_new = {}
    for k in sorted(dict1, key=lambda k: len(k), reverse=True):
        dict_new[k] = dict1[k]
    return dict_new

def get_motif_index(dictionary,sequence):
    index_pos = dictionary[sequence][0].strip('[]').split(',')
    index_neg = dictionary[sequence][1].strip('[]').split(',')
    index_pos = [int(x) for x in index_pos]
    index_neg = [int(x) for x in index_neg]
    return index_pos,index_neg

def get_motif_df(df_pos,index_pos,n1,df_neg,index_neg,n2):
    df_sub1 = extract_subdf(df_pos,index_pos,n1)
    df_sub2 = extract_subdf(df_neg,index_neg,n2)
    df_sub = pd.concat([df_sub1,df_sub2])
    df_sub.reset_index(level=None, drop=True, inplace=True, col_level=0, col_fill='')
    return df_sub


def encoding(df, bases_list):
    data = np.stack(df.loc[:, 'Sequence'].astype('str').apply(list).values, axis=0)
    e = LabelEncoder()
    e.fit(bases_list)
    encoded = e.transform(data.reshape(-1)).reshape(-1, len(data[0]))
    return encoded

def pca_pep(df_pca, encoded, n=2):
    pca = PCA(n_components=n)
    principalComponents = pca.fit_transform(encoded)
    df_pca['pca_1'] = principalComponents[:, 0]
    df_pca['pca_2'] = principalComponents[:, 1]
    #df_pca['pca_3'] = principalComponents[:, 2]
    return principalComponents, df_pca

def t_sne_pep(df, encoded, n=2):
    t_sne = TSNE(n_components=n)
    t_distributed = t_sne.fit_transform(encoded)
    df['t_SNE_1'] = t_distributed[:, 0]
    df['t_SNE_2'] = t_distributed[:, 1]
    #df['t_SNE_3'] = t_distributed[:, 2]
    return t_distributed, df


bases_list = list('ACGT')
df1 = pd.read_csv('highDff.csv')
df = read_csv('Exp6R.motifs.length')
df_gap = read_csv('Exp6R.motifs.gap_length')
#dict1 = search_candidates(df1,df)
dict2 = search_candidates(df1,df_gap)
for x, y in dict2.items():
    print(x,len(y))

motif = 'AGCCCTTC'
orig_seq = 'AGCCCTTCACCACCAACT'
index_pos,index_neg = get_motif_index(dict1,motif)
df_pos = pd.read_csv('Exp6R_n.csv')
df_neg = pd.read_csv('Ctrl6R_n.csv')
df_sub = get_motif_df(df_pos,index_pos,1,df_neg,index_neg,2)
df_sub = mark_origial_seq(df_sub,orig_seq,3,10,20)

#encoded = encoding(df, bases_list)
t_distributed, df_sne =  t_sne_pep(df_sub, encoded_sub, n=2)
principalComponents, df_pca = pca_pep(df_sub, encoded_sub, n=2)

#plot 2d pca figure
fig = px.scatter(df_pca, x="pca_1", y="pca_2", hover_name="Class", color = "Class", size_max=60)
fig.update_layout(height=800)
fig.show()
#fig.write_image('{}_pca.png'.format(motif))

#plot 2d tsne figure
fig = px.scatter(df_sne, x="t_SNE_1", y="t_SNE_2", hover_name="Class", color = "Class", size_max=60)
fig.update_layout(height=800)
fig.show()
#fig.write_image('{}_sne.png'.format(motif))