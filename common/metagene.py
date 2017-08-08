"""
Given a metagene file (a list of genes) and a cancer type, return a dataframe with
that metagene.

Data is read from rnaseq/cancer_type
"""
import pandas as pd
import numpy as np
import os
import glob

import utilities as util

def get_metagene_data(metagene_file, cancer_type):
  rnaseq_glob = os.path.join('rnaseq', cancer_type + '*.txt')
  rnaseq_file = glob.glob(rnaseq_glob)
  assert(len(rnaseq_file) == 1)
  rnaseq = pd.read_csv(rnaseq_file[0], sep='\t', low_memory=False, index_col=0)

  metagene_list = pd.read_csv(metagene_file)
  metagene_df = rnaseq.loc[metagene_list['RNASeq']].astype(float)
  metagene_df = metagene_df.transpose().reset_index()
  metagene_df = util.maybe_clear_non_01s(metagene_df, 'index', cancer_type)
  metagene_df = util.add_identifier_column(metagene_df, 'index').drop('index', axis=1)
  metagene_df = metagene_df.set_index('identifier')
  metagene_df = metagene_df.transpose()

  # now we normalize. take the mean of the base 2 log, then subtract that from each row.
  # then take the average across genes to get the metagene value
  metagene_df_log2 = metagene_df.apply(np.log2)
  metagene_clipped = np.clip(metagene_df_log2, 0, np.inf)
  metagene_means = metagene_clipped.mean(axis=1)

  metagene_normed = metagene_clipped.sub(metagene_means, axis=0)
  metagene = metagene_normed.mean()
  metagene.name = 'metagene'
  return metagene


