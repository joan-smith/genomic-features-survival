#!/usr/bin/env python
# encoding: utf-8
'''
Created by Joan Smith
on 2018-9-15.

Zscores for frequently mutated codons

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
import mutation_base

def get_options():
  parser = argparse.ArgumentParser(description='Get mutation and clinical dir')
  parser.add_argument('-m', action='store', dest='mutation_directory')
  parser.add_argument('-c', action='store', dest='clinical_directory')
  parser.add_argument('-d', action='store', dest='codon_count_directory')
  parser.add_argument('-f', action='store', dest='codon_count_file')
  parser.add_argument('-p', action='store', dest='cutoff_percent', type=float, default=0.05)
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  ns = parser.parse_args()

  return (ns.mutation_directory, ns.clinical_directory,
          ns.codon_count_directory, ns.codon_count_file,
          ns.cutoff_percent, ns.output_directory)

def prep_codon_list(codon_counts, codon_key):
  codon_counts[['Gene','Codon']] = codon_counts['Gene-Codon'].str.split('-', expand=True)
  codon_counts = codon_counts.set_index('Gene-Codon')

  codon_key = codon_counts.join(codon_key)
  codon_key[['chr', 'start_pos']] = codon_key['pos'].str.split(':', expand=True)
  codon_key['chr'] = codon_key['chr'].str[3:]
  codon_key['pos'] = codon_key['pos'].str[3:]
  codon_key = codon_key.reset_index()
  codon_key = codon_key.set_index('pos')

  return codon_key

def get_start_pos(df):
  upper_columns =  [i.upper() for i in df.columns]
  start_pos_index = upper_columns.index('START_POSITION')
  return df.columns[start_pos_index]

def mutations_for_gene(df):
  mutated_patients = df['identifier'].unique()
  return pd.DataFrame({'mutated': np.ones(len(mutated_patients))}, index=mutated_patients)

def make_cox(g, clinical_data, cancer_type, cutoff_percent, seq_patients, outdir):
  gene_codon = g['Gene-Codon'].iloc[0]

  by_patient = g.reset_index()
  by_patient = by_patient.pivot(index='level_3', columns='index', values='mutated')
  clinical_data = clinical_data.loc[seq_patients]
  num_seq_mutated = len(by_patient)
  if num_seq_mutated <= cutoff_percent * seq_patients.size:
    return None
  print gene_codon
  print num_seq_mutated, seq_patients.size, cutoff_percent * seq_patients.size
  print 'num patients w mut in more than 1 codon:', (by_patient.sum(axis=1) > 1).sum()
  by_patient['any_mut'] = by_patient.sum(axis=1) >= 1
  by_patient = by_patient.join(clinical_data, how='outer')
  by_patient['any_mut'] = by_patient[['any_mut']].fillna(0).astype(int)
  by_patient = by_patient.dropna(subset=['time', 'censor'], how='any')

  cox_dict = analysis.do_cox(by_patient.time, by_patient.censor, by_patient.any_mut)
  cox_dict['cancer_type'] = cancer_type
  cox_dict['num_mutated_w_clinical'] = by_patient['any_mut'].sum()
  cox_dict['num_sequence_mutated'] = num_seq_mutated
  by_patient.to_csv(os.path.join(outdir, cancer_type + '_' + gene_codon + '_' + str(cutoff_percent) + 'cutoff_clinical.csv'))
  return pd.Series(cox_dict)

def collect_mutation_data(mutation, clinical, genes):
  cancer_type = util.get_cancer_type(clinical)
  start_pos = ''
  chromosome = 'Chromosome'
  mutation_df = pd.DataFrame()
  if cancer_type in ['COADREAD', 'OV']:
    print 'Using translated NCBI build', cancer_type
    folder = os.path.dirname(mutation)
    new_path = os.path.join(folder, 'HG36_HG37', cancer_type + '_hg36_hg37.txt')
    mutation_df = mutation_base.prep_mutation_data_alone(new_path)
    start_pos = u'hg37_start'
    chromsome = u'hg37_chr'
  else:
    mutation_df = mutation_base.prep_mutation_data_alone(mutation)
    start_pos = get_start_pos(mutation_df)

  seq_patients = mutation_df['identifier'].unique()
  mutation_df = mutation_df[mutation_df['Hugo_Symbol'].isin(genes)]

  # include only missense mutations
  mutation_df = mutation_df[mutation_df[u'Variant_Classification'].str.contains('Missense')]
  mutation_df = mutation_df.reset_index()

  codon_mutation_df = mutation_df.groupby([chromosome, 'Hugo_Symbol', start_pos]).apply(mutations_for_gene)
  codon_mutation_df = codon_mutation_df.reset_index()
  codon_mutation_df['chr-pos'] = codon_mutation_df[chromosome] + ':' + codon_mutation_df[start_pos]
  codon_mutation_df.set_index('chr-pos', inplace=True)

  return codon_mutation_df, seq_patients


def main(argv=None):
  mutation_dir, clinical_dir, codon_count_dir, codon_count_file, cutoff_percent, outdir = get_options()

  codon_counts = pd.read_csv(codon_count_file)
  codon_key = pd.read_csv(os.path.join(codon_count_dir, 'codon key.csv'), index_col=0, header=None,
                                       names=['Gene-Codon', 'pos'])
  codon_key = prep_codon_list(codon_counts, codon_key)
  interesting_genes = codon_key['Gene'].values
  print codon_key[codon_key['chr'].isnull()]

  clinical_files = os.listdir(clinical_dir)
  clinical_files = util.remove_extraneous_files(clinical_files)
  clinical_files = [os.path.join(clinical_dir, f) for f in clinical_files]

  results = []
  my_codon_counts = []
  num_seq_patients = []
  for clinical in clinical_files:
    cancer_type = util.get_cancer_type(clinical)
    print cancer_type
    if '-' in cancer_type:
      continue
    mutation = glob.glob(os.path.join(mutation_dir, '*' + cancer_type + '*'))[0]
    clinical_data = util.get_clinical_data(clinical)

    codon_mutation_df, seq_patients = collect_mutation_data(mutation, clinical, interesting_genes)
    num_seq_patients.append({'cancer_type': cancer_type, 'num_seq': len(seq_patients)})

    # get mutated patients for interesting codons
    interesting_codons = codon_mutation_df.join(codon_key[['Gene', 'Codon', 'chr', 'Gene-Codon']], how='inner')
    interesting_codons = interesting_codons.reset_index()
    ctype_result = interesting_codons.groupby('Gene-Codon').apply(make_cox,
                                                             clinical_data=clinical_data,
                                                             cancer_type=cancer_type,
                                                             cutoff_percent=cutoff_percent,
                                                             seq_patients=seq_patients,
                                                             outdir=outdir)
    ctype_result = ctype_result.dropna(how='all', axis=0)
    results.append(ctype_result)

  num_sequenced_patients_df = pd.DataFrame(num_seq_patients)
  num_sequenced_patients_df = num_sequenced_patients_df.set_index('cancer_type')
  num_sequenced_patients_df.to_csv(os.path.join(outdir, 'num_sequenced_patients.csv'))

  results_df = pd.concat(results, axis=0, levels='cancer_type')
  results_df.to_csv(os.path.join(outdir, str(cutoff_percent) + 'cutoff_zscores.csv'))

  my_codon_counts_df = results_df.pivot(values='num_sequence_mutated', columns='cancer_type')
  my_codon_counts_df.to_csv(os.path.join(outdir, 'output_codon_counts.csv'))


if __name__ == "__main__":
  main()
