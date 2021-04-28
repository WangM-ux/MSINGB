#!/usr/bin/env python 
# -*- coding:utf-8 -*-
import numpy as np
from ngboost import NGBClassifier
from sklearn.feature_selection import RFECV, RFE
from sklearn.model_selection import KFold
from sklearn.svm import SVC
from xgboost import XGBClassifier

from MSINGB import get_feature

ngbc_learning_rate = 0.01
ngbc_n_estimators = 450
ngbc_minibatch_frac = 1.0
ngbc_tol = 0.0016
MSI_feature=['SNP','INDEL','SNP_R','INDEL_R','t_mutation','t_mutation_R','SNP_R/SNP','INDEL_R/INDEL',
             'tm_R/tm','Frame_Shift_Del','Frame_Shift_Ins','In_Frame_Del','In_Frame_Ins',
             'Missense_Mutation','Nonsense_Mutation','Silent','Splice_Site','3\'UTR',
             '3\'Flank','5\'UTR','5\'Flank','Intron']

def rfe_select(feature_file, outPath):

    outfile = outPath + 'feature_sel.txt'
    outpos = outPath + '_pos.txt'
    # Load boston housing dataset as an example
    targetList, featureMatr = get_feature(feature_file)
    train_data = featureMatr
    train_label = np.array(targetList).astype(int)
    X = train_data
    Y = train_label
    print('the dimsion of features is: ' + str(len(X[0])))
    print('the sum of samples is: ' + str(len(X)))
    svc = XGBClassifier()
    # svc = SVC(kernel="linear")
    # svc = NGBClassifier(learning_rate=ngbc_learning_rate, n_estimators=ngbc_n_estimators, minibatch_frac=ngbc_minibatch_frac,tol=ngbc_tol)
    rfecv = RFECV(estimator=svc, cv=KFold(n_splits=5, shuffle=True, random_state=1), verbose=1,
                  scoring='neg_mean_squared_error')
    # svc = SVC(kernel="linear")
    # rfecv = RFE(estimator=svc, n_features_to_select=1, step=1)
    print('开始训练特征')
    rfecv.fit(X, Y)
    print('训练特征选择完成，输出特征：')
    importances = rfecv.support_
    posarr = []
    # print('***********************grid_scores_***********************')
    # print(rfecv.grid_scores_)
    print('***********************n_features_***********************')
    print(rfecv.n_features_)
    print('***********************ranking_***********************')
    print(rfecv.ranking_)
    print('***********************support_***********************')
    print(rfecv.support_)

    # print('***********************result***********************')
    # result = X.columns[rfecv.get_support()]
    # print(result)
    print('特征输出完毕')
    print('***********************保存***********************')
    for i in range(X.shape[1]):
        if importances[i] == True:
            posarr.append(i)


    feature_sel = X[:, posarr]
    feature_sel_Tag = np.c_[Y, feature_sel]
    np.savetxt(outfile, feature_sel_Tag, delimiter=',')
    print(posarr)
    feature = np.array(MSI_feature)
    opti_feature = feature[posarr]
    print(opti_feature)
    pos_feature = np.c_[posarr, opti_feature]
    print(pos_feature)
    np.savetxt(outpos, pos_feature, fmt="%s")
    return outpos, posarr

if __name__ == '__main__':
    rfe_select('../MSINGB/train/all_train.csv', '../output/10_RFE_')
