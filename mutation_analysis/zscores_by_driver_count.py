#!/usr/bin/env python
# encoding: utf-8
'''
zscores_by_driver_count.py

Created by Joan Smith
on 2017-08-03.

Given a clinical file and a mutations file, calculate zscores for genes by driver mutation count.

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
  namespace = parser.parse_args()

  return (namespace.mutation_directory, namespace.clinical_directory,
          namespace.output_directory)

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

def calculate_cox(mutation, clinical_data, outdir):
  df = prep_data(mutation, clinical_data)
  df = df.join(clinical_data, how='inner')
  num_patients = len(df.index)

  #prep output file
  cancer_type = os.path.basename(mutation).split('_')[0].split('.')[0]
  print cancer_type
  outfile = os.path.join(outdir, cancer_type + '.driver_mutation_count.zscores.out.csv')

  print 'Missing driver genes:', set(COMMONLY_MUTATED) - set(df.columns)

  present_driver_genes = list(set(df.columns).intersection(set(COMMONLY_MUTATED)))
  print present_driver_genes
  driver_mutations = df[present_driver_genes]
  print driver_mutations
  driver_mutations['driver_mutation_count'] = driver_mutations.sum(axis=1, skipna=True)
  driver_mutations['time'] = df['time']
  driver_mutations['censor'] = df['censor']

  analysis_data = pd.DataFrame()
  analysis_data['time'] = driver_mutations['time']
  analysis_data['censor'] = driver_mutations['censor']
  analysis_data['driver_mutation_count'] = driver_mutations['driver_mutation_count']

  #Do analysis!
  cox_dict = analysis.do_cox(analysis_data['time'], analysis_data['censor'], analysis_data['driver_mutation_count'])
  with open(outfile, 'w') as out:
    out.write('Z: ' + str(cox_dict['z']) + ', P: ' + str(cox_dict['p']) + ', n: ' + str(cox_dict['n']) + '\n')
    driver_mutations.to_csv(out)
  return cox_dict

def multiprocess_zscores(args):
  mutation = args[0]
  clinical_data = args[1]
  outdir = args[2]
  calculate_cox(mutation, clinical_data, outdir)


def main(argv=None):
  mutation_dir, clinical_dir, outdir = get_options()

  clinical_files = os.listdir(clinical_dir)
  clinical_files = util.remove_extraneous_files(clinical_files)
  clinical_files = [os.path.join(clinical_dir, f) for f in clinical_files]
  p = Pool(1)

  args = []
  pancan = {}
  for clinical in clinical_files:
    cancer_type = util.get_cancer_type(clinical)
    print cancer_type
    mutation = glob.glob(os.path.join(mutation_dir, '*' + cancer_type + '*'))[0]

    clinical_data = util.get_clinical_data(clinical)
    #args.append((mutation, clinical_data, outdir))
    pancan[cancer_type] = calculate_cox(mutation, clinical_data, outdir)

  #print args
  #p.map(multiprocess_zscores, args)
  pancan_df = pd.DataFrame(pancan)
  pancan_df = pancan_df.transpose()
  pancan_df.to_csv(os.path.join(outdir, 'pancan.csv'))


if __name__ == "__main__":
  main()
