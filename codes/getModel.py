#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np

from ngboost import NGBClassifier
from mlxtend.feature_selection import ColumnSelector
from sklearn.pipeline import make_pipeline

def get_feature(fileName):

    SNP_Dic={}
    INDEL_Dic={}
    SNP_R_Dic={}
    INDEL_R_Dic={}
    t_mutation_Dic={}
    t_mutation_R_Dic={}
    SNP_R_SNP_Dic={}
    INDEL_R_INDEL_Dic={}
    tm_R_tm_Dic={}
    Frame_Shift_Del_Dic={}
    Frame_Shift_Ins_Dic={}
    In_Frame_Del_Dic={}
    In_Frame_Ins_Dic={}
    Missense_Mutation_Dic={}
    Nonsense_Mutation_Dic={}
    Silent_Dic={}
    Splice_Site_Dic={}
    UTR_3_Dic={}
    Flank_3_Dic={}
    UTR_5_Dic={}
    Flank_5_Dic={}
    Intron_Dic={}
    MSI_status_Dic={}

    with open(fileName,'r') as fp:
        fp.readline()
        lines = fp.readlines()
        for line in lines:
            line=line.strip()
            if(line != ''):
                arrs = line.split(',')
                Sample_id=arrs[0].strip()
                SNP_Dic[Sample_id]=arrs[2].strip()
                INDEL_Dic[Sample_id]=arrs[3].strip()
                SNP_R_Dic[Sample_id]=arrs[4].strip()
                INDEL_R_Dic[Sample_id]=arrs[5].strip()
                t_mutation_Dic[Sample_id]=arrs[6].strip()
                t_mutation_R_Dic[Sample_id]=arrs[7].strip()
                SNP_R_SNP_Dic[Sample_id]=arrs[8].strip()
                INDEL_R_INDEL_Dic[Sample_id]=arrs[9].strip()
                tm_R_tm_Dic[Sample_id]=arrs[10].strip()
                Frame_Shift_Del_Dic[Sample_id]=arrs[11].strip()
                Frame_Shift_Ins_Dic[Sample_id]=arrs[12].strip()
                In_Frame_Del_Dic[Sample_id]=arrs[13].strip()
                In_Frame_Ins_Dic[Sample_id]=arrs[14].strip()
                Missense_Mutation_Dic[Sample_id]=arrs[15].strip()
                Nonsense_Mutation_Dic[Sample_id]=arrs[16].strip()
                Silent_Dic[Sample_id]=arrs[17].strip()
                Splice_Site_Dic[Sample_id]=arrs[18].strip()
                UTR_3_Dic[Sample_id]=arrs[19].strip()
                Flank_3_Dic[Sample_id]=arrs[20].strip()
                UTR_5_Dic[Sample_id]=arrs[21].strip()
                Flank_5_Dic[Sample_id]=arrs[22].strip()
                Intron_Dic[Sample_id]=arrs[23].strip()
                MSI_status_Dic[Sample_id]=arrs[24].strip()

    probMatr = np.zeros((len(Intron_Dic),23))

    Sample_id_List = Intron_Dic.keys()
    sampleNo = 0
    for Sample_id in Sample_id_List:
        if(MSI_status_Dic[Sample_id]=='MSS' or MSI_status_Dic[Sample_id]=='MSI-L'):
            probMatr[sampleNo][0]=0
        elif(MSI_status_Dic[Sample_id]=='MSI-H'):
            probMatr[sampleNo][0]=1

        probMatr[sampleNo][1]=SNP_Dic[Sample_id]
        probMatr[sampleNo][2]=INDEL_Dic[Sample_id]
        probMatr[sampleNo][3]=SNP_R_Dic[Sample_id]
        probMatr[sampleNo][4]=INDEL_R_Dic[Sample_id]
        probMatr[sampleNo][5]=t_mutation_Dic[Sample_id]
        probMatr[sampleNo][6]=t_mutation_R_Dic[Sample_id]
        probMatr[sampleNo][7]=SNP_R_SNP_Dic[Sample_id]
        probMatr[sampleNo][8]=INDEL_R_INDEL_Dic[Sample_id]
        probMatr[sampleNo][9]=tm_R_tm_Dic[Sample_id]
        probMatr[sampleNo][10]=Frame_Shift_Del_Dic[Sample_id]
        probMatr[sampleNo][11]=Frame_Shift_Ins_Dic[Sample_id]
        probMatr[sampleNo][12]=In_Frame_Del_Dic[Sample_id]
        probMatr[sampleNo][13]=In_Frame_Ins_Dic[Sample_id]
        probMatr[sampleNo][14]=Missense_Mutation_Dic[Sample_id]
        probMatr[sampleNo][15]=Nonsense_Mutation_Dic[Sample_id]
        probMatr[sampleNo][16]=Silent_Dic[Sample_id]
        probMatr[sampleNo][17]=Splice_Site_Dic[Sample_id]
        probMatr[sampleNo][18]=UTR_3_Dic[Sample_id]
        probMatr[sampleNo][19]=Flank_3_Dic[Sample_id]
        probMatr[sampleNo][20]=UTR_5_Dic[Sample_id]
        probMatr[sampleNo][21]=Flank_5_Dic[Sample_id]
        probMatr[sampleNo][22]=Intron_Dic[Sample_id]
        sampleNo += 1

    targetList = probMatr[:,0]
    featureMatr = probMatr[:,1:]

    return targetList,featureMatr

def get_test_feature(feature):

    SNP_Dic={}
    INDEL_Dic={}
    SNP_R_Dic={}
    INDEL_R_Dic={}
    t_mutation_Dic={}
    t_mutation_R_Dic={}
    SNP_R_SNP_Dic={}
    INDEL_R_INDEL_Dic={}
    tm_R_tm_Dic={}
    Frame_Shift_Del_Dic={}
    Frame_Shift_Ins_Dic={}
    In_Frame_Del_Dic={}
    In_Frame_Ins_Dic={}
    Missense_Mutation_Dic={}
    Nonsense_Mutation_Dic={}
    Silent_Dic={}
    Splice_Site_Dic={}
    UTR_3_Dic={}
    Flank_3_Dic={}
    UTR_5_Dic={}
    Flank_5_Dic={}
    Intron_Dic={}
    for index, arrs in feature.iterrows():
        Sample_id=index
        SNP_Dic[Sample_id]=arrs[0]
        INDEL_Dic[Sample_id]=arrs[1]
        SNP_R_Dic[Sample_id]=arrs[2]
        INDEL_R_Dic[Sample_id]=arrs[3]
        t_mutation_Dic[Sample_id]=arrs[4]
        t_mutation_R_Dic[Sample_id]=arrs[5]
        SNP_R_SNP_Dic[Sample_id]=arrs[6]
        INDEL_R_INDEL_Dic[Sample_id]=arrs[7]
        tm_R_tm_Dic[Sample_id]=arrs[8]
        Frame_Shift_Del_Dic[Sample_id]=arrs[9]
        Frame_Shift_Ins_Dic[Sample_id]=arrs[10]
        In_Frame_Del_Dic[Sample_id]=arrs[11]
        In_Frame_Ins_Dic[Sample_id]=arrs[12]
        Missense_Mutation_Dic[Sample_id]=arrs[13]
        Nonsense_Mutation_Dic[Sample_id]=arrs[14]
        Silent_Dic[Sample_id]=arrs[15]
        Splice_Site_Dic[Sample_id]=arrs[16]
        UTR_3_Dic[Sample_id]=arrs[17]
        Flank_3_Dic[Sample_id]=arrs[18]
        UTR_5_Dic[Sample_id]=arrs[19]
        Flank_5_Dic[Sample_id]=arrs[20]
        Intron_Dic[Sample_id]=arrs[21]

    probMatr = np.zeros((len(Intron_Dic),23))
    Sample_id_List = Intron_Dic.keys()
    sampleNo = 0
    for Sample_id in Sample_id_List:
        probMatr[sampleNo][1] = SNP_Dic[Sample_id]
        probMatr[sampleNo][2] = INDEL_Dic[Sample_id]
        probMatr[sampleNo][3] = SNP_R_Dic[Sample_id]
        probMatr[sampleNo][4] = INDEL_R_Dic[Sample_id]
        probMatr[sampleNo][5] = t_mutation_Dic[Sample_id]
        probMatr[sampleNo][6] = t_mutation_R_Dic[Sample_id]
        if SNP_R_SNP_Dic[Sample_id] == '':
            probMatr[sampleNo][7] = 0
        else:
            probMatr[sampleNo][7] = SNP_R_SNP_Dic[Sample_id]
        if (INDEL_R_INDEL_Dic[Sample_id] != ''):
            probMatr[sampleNo][8] = INDEL_R_INDEL_Dic[Sample_id]
        probMatr[sampleNo][9] = tm_R_tm_Dic[Sample_id]
        probMatr[sampleNo][10] = Frame_Shift_Del_Dic[Sample_id]
        probMatr[sampleNo][11] = Frame_Shift_Ins_Dic[Sample_id]
        probMatr[sampleNo][12] = In_Frame_Del_Dic[Sample_id]
        probMatr[sampleNo][13] = In_Frame_Ins_Dic[Sample_id]
        probMatr[sampleNo][14] = Missense_Mutation_Dic[Sample_id]
        probMatr[sampleNo][15] = Nonsense_Mutation_Dic[Sample_id]
        probMatr[sampleNo][16] = Silent_Dic[Sample_id]
        probMatr[sampleNo][17] = Splice_Site_Dic[Sample_id]
        probMatr[sampleNo][18] = UTR_3_Dic[Sample_id]
        probMatr[sampleNo][19] = Flank_3_Dic[Sample_id]
        probMatr[sampleNo][20] = UTR_5_Dic[Sample_id]
        probMatr[sampleNo][21] = Flank_5_Dic[Sample_id]
        probMatr[sampleNo][22] = Intron_Dic[Sample_id]
        sampleNo += 1

    targetList = Sample_id_List
    featureMatr = probMatr[:,1:]

    return targetList,featureMatr

def calculate_MSI(results, test_targetList, cutoff):
    tumorNo = 0
    res_list = []
    tumorList = []
    for Tumor_Id in test_targetList:
        tumorList.append(Tumor_Id)
        prepro = results[tumorNo]
        if prepro > cutoff:
            MSI_stauts = 'MSI-H'
        else:
            MSI_stauts = 'MSS'
            prepro = 1-prepro
        prepro = "%.2f%%" % (prepro * 100)
        res_list.append([MSI_stauts, prepro])
        tumorNo += 1

    res_arr = np.array(res_list)
    df2 = pd.DataFrame(res_arr, index=tumorList, columns=["Pre_MSI_Status", "Pre_Probability"])
    return df2

def model_pre(trainSet, inputFeature, output):
    # Optimal feature
    ngbc_cols_t = [9, 15, 14, 19, 12, 18, 1, 16, 11, 8, 0, 10, 6, 17]
    # Optimal parameters of NGBoost
    ngbc_learning_rate = 0.01
    ngbc_n_estimators = 450
    ngbc_minibatch_frac = 1.0
    ngbc_tol = 0.0016
    cutoff = 0.3
    # model
    ngbc = make_pipeline(ColumnSelector(cols=ngbc_cols_t),
                         NGBClassifier(learning_rate=ngbc_learning_rate,
                                       n_estimators=ngbc_n_estimators,
                                       minibatch_frac=ngbc_minibatch_frac,
                                       tol=ngbc_tol
                                       ))

    rs_path = output + '/preres.csv'
    # Training dataset
    targetList, featureMatr = get_feature(trainSet)
    train_data = featureMatr
    train_label = np.array(targetList).astype(int)
    # Prediction dataset
    test_targetList,test_featureMatr = get_test_feature(inputFeature)
    predict_data = test_featureMatr
    # Fit the algorithm on the data
    MSImodel = ngbc.fit(train_data, train_label)
    y_predproba = MSImodel.predict_proba(predict_data)
    # Calculate the performance predicted by the model
    results = y_predproba[:, 1]
    res = calculate_MSI(results,test_targetList, cutoff)
    res.to_csv(rs_path)




