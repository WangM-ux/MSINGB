#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import numpy as np
import pandas as pd
from intervaltree import Interval, IntervalTree

from getFeature import make_feature_table

column_list =[	"bin","chrom","chromStart","chromEnd","name_tag",
        	    		"period_size","copyNUM","consensusSize","perMatch",
        	       		"perIndel","score","A","C","G","T","entropy","unit_sequence"]


candidate_chrom = [ 'chr1', 'chr2', 'chr3', 'chr4', 'chr5', 'chr6', 'chr7', 'chr8',
        	            	'chr9', 'chr10', 'chr11', 'chr12', 'chr13', 'chr14', 'chr15',
        	            	'chr16', 'chr17', 'chr18', 'chr19', 'chr20', 'chr21', 'chrX',
        	            	'chr22', 'chrY']


def getFile(fileUrl):
    chunksize = 10000
    chunks = []
    table_reader = pd.read_table(fileUrl, low_memory=False, comment='#', chunksize=chunksize)
    for chunk in table_reader:
        chunks.append(chunk)

    table_reader.close()
    table = pd.concat(chunks, axis=0)
    return table

def create_repeats_tree(chromosome,ref_repeats_df):
    target_ref_repeats = ref_repeats_df[ref_repeats_df['chrom']==chromosome]
    interval_tuples = target_ref_repeats.apply(lambda row: (row['chromStart'],row['chromEnd']+1),axis=1)
    target_interval_tree = IntervalTree.from_tuples(interval_tuples)
    return target_interval_tree


def tag_maf_row(maf_row,repeats_tree):
    query_tuple = (maf_row['Start_Position'],maf_row['End_Position']+1)
    if len(repeats_tree[query_tuple[0]:query_tuple[1]])>=1:
        maf_row['In_repeats'] = 1
        return maf_row
    elif len(repeats_tree[query_tuple[0]:query_tuple[1]])==0:
        maf_row['In_repeats'] = 0
        return maf_row


def tag_maf_table(maf_table, ref_repeats):

    tagged_group_frame = []
    grouped_maf_file = maf_table.groupby('Chromosome')

    for name, group_df in grouped_maf_file:
        ref_repeats_tree = create_repeats_tree(chromosome=name, ref_repeats_df=ref_repeats)
        tagged_group_df = group_df.apply(lambda row: tag_maf_row(row, ref_repeats_tree), axis=1)
        tagged_group_frame.append(tagged_group_df)
    return pd.concat(tagged_group_frame, ignore_index=True, axis=0)
    pass


def getMSI_tag(fileUrl, SRS_file):

    # 读取SRS文件
    ref_repeats = pd.read_table(SRS_file, names=column_list)
    ref_repeats = ref_repeats[(ref_repeats['chrom'].isin(candidate_chrom)) & (ref_repeats['period_size'] <= 5)][
        ['chrom', 'chromStart', 'chromEnd']]

    maf_table = getFile(fileUrl)

    MSI_annotated = tag_maf_table(maf_table, ref_repeats)
    return MSI_annotated

def setMSI_file(tumorMSI, param):
    pass




class Test_Input(object):
    def __init__(self, inputMaf, output):
        self.rawfile = inputMaf
        self.outDir = output

    def creat_tagged_feature(self,ref_repeats_file):
        rawfile = self.rawfile
        output = self.outDir
        tag_path = output + '/tagged.csv'
        feature_path = output + '/feature.csv'

        print('**************Add In_repeats is annotated**************')
        MSI_Tag = getMSI_tag(rawfile, ref_repeats_file)
        MSI_Tag.to_csv(tag_path.strip(), sep='\t', index=False)
        print('*****************Extracting 22 features*****************')
        MSI_Feature = make_feature_table(exome_size=44,tagged_maf_file = tag_path.strip())
        MSI_Feature.to_csv(feature_path)
        return MSI_Feature
