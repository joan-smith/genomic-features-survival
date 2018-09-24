#!/usr/bin/env python
# encoding: utf-8
'''
variant_allele_freq.py

Created by Joan Smith
on 2017-09-02.

Given the set of mutation files and the variant allele frequency key, calculate variante allel frequency distributions.

Copyright (c) 2018. All rights reserved.
'''

import pandas as pd
import numpy as np
import argparse
import sys
import os
import glob
import itertools
import pdb
from multiprocessing import Pool

from matplotlib import pyplot

sys.path.append('../common/')
import utilities as util
import analysis

MUTATION_PERCENT = .02

COMMONLY_MUTATED = ['\'TP53', '\'KRAS', '\'PIK3CA', '\'APC', '\'KMT2D', '\'ARID1A', '\'PTEN', '\'BRAF', '\'ATM',
                    '\'EGFR', '\'NF1', '\'RB1', '\'BRCA2', '\'ATRX', '\'NOTCH1', '\'CDKN2A', '\'SETD2', '\'CREBBP',
                    '\'SMAD4', '\'FBXW7', '\'ARID1B', '\'SMARCA4', '\'KMT2A', '\'EP300', '\'ERBB4', '\'IDH1',
                    '\'ARID2', '\'NRAS', '\'ROS1', '\'CTNNB1']

def get_options():
  parser = argparse.ArgumentParser(description='Get mutation and clinical dir')
  parser.add_argument('-i', action='store', dest='mutation_directory')
  parser.add_argument('-k', action='store', dest='key_file')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  namespace = parser.parse_args()

  return (namespace.mutation_directory, namespace.key_file,
          namespace.output_directory)

def prep_data(mutation, key):
  df = pd.read_csv(mutation, sep='\t', low_memory=False, dtype=str)
  cancer_type = util.get_cancer_type(mutation)

  # remove column headers from combined mutation sheet
  df = df[~df[u'Hugo_Symbol'].str.contains('Hugo_Symbol')]
  df['Hugo_Symbol'] = '\'' + df['Hugo_Symbol'].astype(str)
  df = df[df[u'Hugo_Symbol'].isin(COMMONLY_MUTATED)]

  df[u'Tumor_Sample_Barcode'] = df[u'Tumor_Sample_Barcode'].str.strip()

  number_barcodes_in_mutation_data = df[u'Tumor_Sample_Barcode'].unique().size
  print 'Number of total sequenced barcodes:   ', number_barcodes_in_mutation_data
  df = util.maybe_clear_non_01s(df, u'Tumor_Sample_Barcode', cancer_type)
  df = util.add_identifier_column(df, u'Tumor_Sample_Barcode')

  # include only nonsilent mutations
  non_silent = df.where(df[u'Variant_Classification'] != 'Silent')
  df = non_silent.dropna(subset=[u'Variant_Classification'])
  df = df.reset_index()

  df['VAF'] = calculate_vaf(df, key.loc[cancer_type])
  # use the largest VAF
  df = df.groupby(['Hugo_Symbol', 'identifier']).max()
  df = df.reset_index()
  pivoted = df.pivot(index='identifier', columns='Hugo_Symbol', values='VAF')

  minimum_vaf_count = MUTATION_PERCENT * number_barcodes_in_mutation_data
  enough_patients = pivoted.count() >= minimum_vaf_count
  too_few_patients = enough_patients[~enough_patients].index.values
  print 'Genes with too few patients:', too_few_patients
  pivoted = pivoted.drop(too_few_patients, axis=1)
  return pivoted

def calculate_vaf(df, key):
  if not pd.isnull(key['VAF?']):
    vaf = df[key['VAF?']].str[:-1]
    vaf = vaf.replace(r'\s+', np.nan, regex=True)
    vaf = vaf.replace(r'', np.nan, regex=True)
    return vaf.astype(float) / 100

  alt = df[key['Alt?']]
  if not pd.isnull(key['Total?']):
    total = df[key['Total?']].astype(float)
    print total
  else:
    total = df[key['Alt?']].astype(float) + df[key['Ref?']].astype(float)
  alt = alt.rename('alt')
  total = total.rename('total')
  var_allele_freq =  pd.DataFrame([alt, total], dtype=float).transpose()
  return var_allele_freq['alt'] / var_allele_freq['total']

def calculate_variant_allele_distribution(cancer_type, mutation, key, outdir):
  df = prep_data(mutation, key)
  boxplot_data = []
  boxplot_labels = []

  for gene in df:
    print gene
    boxplot_data.append(df[gene].dropna())
    boxplot_labels.append(gene)
  fig, ax = pyplot.subplots()
  pyplot.title(cancer_type)
  pyplot.boxplot(boxplot_data, labels=boxplot_labels)
  pyplot.setp(ax.get_xticklabels(), rotation=90, horizontalalignment='center')
  pyplot.savefig(cancer_type + '.png', pad_inches=1, bbox_inches='tight')

  df.to_csv(cancer_type + '_vaf_distribition.csv')


def main(argv=None):
  mutation_dir, key_file, outdir = get_options()
  mutation_files = glob.glob(mutation_dir + '*txt')
  key = pd.read_csv(key_file, na_values=['-'], index_col=0)
  key = key.dropna(how='all')
  print key

  p = Pool(1)

  args = []
  pancan = {}
  for mutation in mutation_files:
    cancer_type = util.get_cancer_type(mutation)
    if cancer_type in key.index:
      print cancer_type
      pancan[cancer_type] = calculate_variant_allele_distribution(cancer_type, mutation, key, outdir)

if __name__ == "__main__":
  main()
