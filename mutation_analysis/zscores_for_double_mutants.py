#!/usr/bin/env python
# encoding: utf-8
'''
zscores_for_double_mutants.py

Created by Joan Smith
on 2017-6-20.

Given a clinical file and a mutations file, calculate zscores for genes with mutations in >10% of patients

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
  parser.add_argument('-c', action='store', dest='clinical_directory')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  parser.add_argument('-u', action='store', dest='univariate_output', default='.')
  namespace = parser.parse_args()

  return (namespace.mutation_directory, namespace.clinical_directory,
          namespace.output_directory, namespace.univariate_output)

def prep_data(mutation, clinical_data):
  df = pd.read_csv(mutation, sep='\t', low_memory=False, dtype=str)
  cancer_type = util.get_cancer_type(mutation)

  # remove column headers from combined mutation sheet
  df = df[~df[u'Hugo_Symbol'].str.contains('Hugo_Symbol')]
  df[u'Tumor_Sample_Barcode'] = df[u'Tumor_Sample_Barcode'].str.strip()

  number_barcodes_in_mutation_data = df[u'Tumor_Sample_Barcode'].unique().size
  print 'Number of total sequenced barcodes:   ', number_barcodes_in_mutation_data
  df = util.maybe_clear_non_01s(df, u'Tumor_Sample_Barcode', cancer_type)
  df = util.add_identifier_column(df, u'Tumor_Sample_Barcode')

  # include only nonsilent mutations
  non_silent = df.where(df[u'Variant_Classification'] != 'Silent')
  df = non_silent.dropna(subset=[u'Variant_Classification'])

  df = df.reset_index()
  df['Hugo_Symbol'] = '\'' + df['Hugo_Symbol'].astype(str)

  gene_mutation_df = df.groupby(['Hugo_Symbol']).apply(mutations_for_gene)
  gene_mutation_df.index.set_names(['Hugo_Symbol', 'patient'], inplace=True)
  gene_mutation_df = gene_mutation_df.reset_index()
  gene_patient_mutations = gene_mutation_df.pivot(index='Hugo_Symbol', columns='patient', values='mutated')

  return gene_patient_mutations.transpose()

def mutations_for_gene(df):
  mutated_patients = df['identifier'].unique()
  return pd.DataFrame({'mutated': np.ones(len(mutated_patients))}, index=mutated_patients)

def calculate_cox(mutation, clinical_data, outdir, univariate_file=None):
  df = prep_data(mutation, clinical_data)
  df = df.join(clinical_data, how='inner')
  num_patients = len(df.index)

  gene_pairs = itertools.combinations(COMMONLY_MUTATED, 2)

  #prep output file
  cancer_type = os.path.basename(mutation).split('_')[0].split('.')[0]
  print cancer_type
  outfile = os.path.join(outdir, cancer_type + '_mutation-fraction-' + str(MUTATION_PERCENT) + '.zscores.out.csv')
  if univariate_file:
    univariate_data = pd.read_csv(univariate_file, index_col=0)
    outfile = os.path.join(outdir,
        cancer_type + '_mutation-fraction-' + str(MUTATION_PERCENT) + '.notalonesignificant.zscores.out.csv')
  formatstring = '{0}, {1}, {2}, {3}, {4}\n'

  with open(outfile, 'w') as out:
    out.write('gene,zscore,pvalue,num mutations,num patients\n')

    for gene_pair in gene_pairs:
      gene_pair = pd.Series(gene_pair)
      gene_pair_str = '-'.join(gene_pair)
      if gene_pair.isin(df.columns.values).sum() < 2:
        continue
      paired_mutations = df[list(gene_pair)]
      double_mutated_patients = paired_mutations[paired_mutations.sum(axis=1) == 2].index
      num_mutations = len(double_mutated_patients)
      print gene_pair_str, num_mutations

      if num_mutations >= MUTATION_PERCENT * num_patients:
        # if we have univariate data, check to see that neither gene is significant for survival independently
        #  before calculation
        if univariate_file:
          if univariate_data.loc[gene_pair[0]].zscore < -1.96 or univariate_data.loc[gene_pair[0]].zscore > 1.96:
            print 'Skipping pair', gene_pair_str, 'for gene 0'
            continue
          if univariate_data.loc[gene_pair[1]].zscore < -1.96 or univariate_data.loc[gene_pair[1]].zscore > 1.96:
            print 'Skipping pair', gene_pair_str, 'for gene 1'
            continue

        # analysis_data = pd.DataFrame({'mutated': np.ones(num_mutations)}, index=double_mutated_patients)
        analysis_data = pd.DataFrame()
        analysis_data['time'] = df['time']
        analysis_data['censor'] = df['censor']
        analysis_data['mutated'] = 0
        analysis_data.loc[double_mutated_patients,'mutated'] = 1

        #Do analysis!
        print 'Doing analysis for', gene_pair_str, 'with',  num_mutations, 'double mutations', 'of', num_patients

        name = cancer_type+ '_' + gene_pair_str.replace('\'', '')
        print name
        cox_dict = analysis.do_cox(analysis_data['time'], analysis_data['censor'], analysis_data['mutated'])
        out.write(formatstring.format(gene_pair_str,
                                      cox_dict['z'],
                                      cox_dict['p'],
                                      num_mutations,cox_dict['n']))
        analysis_data.to_csv(os.path.join(outdir, name + '_data.csv'),
                             columns=['time', 'censor', 'mutated'])

def multiprocess_zscores(args):
  mutation = args[0]
  clinical_data = args[1]
  cancer_type_outdir = args[2]
  univariate = args[3]
  calculate_cox(mutation, clinical_data, cancer_type_outdir, univariate_file=univariate)


def main(argv=None):
  mutation_dir, clinical_dir, outdir, univariate_output = get_options()

  clinical_files = os.listdir(clinical_dir)
  clinical_files = util.remove_extraneous_files(clinical_files)
  clinical_files = [os.path.join(clinical_dir, f) for f in clinical_files]
  p = Pool(16)

  args = []
  for clinical in clinical_files:
    cancer_type = util.get_cancer_type(clinical)
    print cancer_type
    mutation = glob.glob(os.path.join(mutation_dir, '*' + cancer_type + '*'))[0]

    univariate_file = None
    if univariate_output:
      univariate_file = glob.glob(os.path.join(univariate_output, cancer_type, cancer_type + '.zscores.out.csv'))[0]
      print univariate_file

    clinical_data = util.get_clinical_data(clinical)
    cancer_type_outdir = os.path.join(outdir, cancer_type)
    if not os.path.isdir(cancer_type_outdir):
      os.makedirs(cancer_type_outdir)
    args.append((mutation, clinical_data, cancer_type_outdir, univariate_file))
    # calculate_cox(mutation, clinical_data, cancer_type_outdir, univariate_file=univariate_file)

  print args
  p.map(multiprocess_zscores, args)


if __name__ == "__main__":
  main()
