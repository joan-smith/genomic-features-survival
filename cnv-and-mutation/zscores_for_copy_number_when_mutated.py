#!/usr/bin/env python
# encoding: utf-8
'''

zscores_for_copy_number_when_mutated.py


Created by Joan Smith
on 2017-7-7.

Given a set of cnv files, mutation files and clinical files, all from TCGA
calculate cnv zscores partitioned by whether or not the gene in question is mutated

Copyright (c) 2017. All rights reserved.

'''

import pandas as pd
import numpy as np
import argparse
import sys
import os
import pdb
import collections
import glob
import rpy2
from multiprocessing import Pool

sys.path.append('../common/')
import utilities as util
import analysis

sys.path.append('../mutation-analysis')
import zscores_for_mutants as mutation_zscores


def get_options():
  parser = argparse.ArgumentParser(description='Get mutation, cnv, and clinical directories. optional output dir')
  parser.add_argument('-i', action='store', dest='cnv_directory')
  parser.add_argument('-m', action='store', dest='mutation_directory')
  parser.add_argument('-c', action='store', dest='clinical_directory')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  namespace = parser.parse_args()

  return (namespace.cnv_directory, namespace.mutation_directory, namespace.clinical_directory,
    namespace.output_directory)

def mutations_for_gene(df):
  mutated_patients = df['identifier'].unique()
  return pd.DataFrame({'mutated': np.ones(len(mutated_patients))}, index=mutated_patients)

def prep_mutation_data(mutation, clinical_data):
  mutation, clinical_data_w_seq_patients, num_patients = mutation_zscores.prep_data(mutation, clinical_data)

  # include only nonsilent mutations
  non_silent = mutation.where(mutation[u'Variant_Classification'] != 'Silent')
  mutation = non_silent.dropna(subset=[u'Variant_Classification'])

  mutation = mutation.reset_index()
  mutation['Hugo_Symbol'] = '\'' + mutation['Hugo_Symbol'].astype(str)

  gene_mutation_df = mutation.groupby(['Hugo_Symbol']).apply(mutations_for_gene)
  gene_mutation_df.index.set_names(['Hugo_Symbol', 'patient'], inplace=True)
  gene_mutation_df = gene_mutation_df.reset_index()
  gene_patient_mutations = gene_mutation_df.pivot(index='Hugo_Symbol', columns='patient', values='mutated')

  return gene_patient_mutations.transpose()

def calculate_cox(data, gene):
  data_cox_dict = collections.defaultdict(lambda: np.nan)
  if data[gene].count() > 10:
    try:
      data_cox_dict = analysis.do_cox(data.time,
                                 data.censor,
                                 data[gene])
    except rpy2.rinterface.RRuntimeError as e:
      print 'WARN: skipped', gene, 'due to R error'
  return data_cox_dict


def make_zscores(copy_number, mutation, clinical, outdir):
  clinical_data = util.get_clinical_data(clinical)
  cnv = pd.read_csv(copy_number, index_col=0)
  mutation = prep_mutation_data(mutation, clinical_data)

  cnv_by_patient = cnv.transpose()
  clinical_and_cnv = cnv_by_patient.join(clinical_data, how='inner')

  clinical_and_mutation_patients = list(set(mutation.index).intersection(set(clinical_and_cnv.index)))
  clinical_and_cnv_with_mutations = clinical_and_cnv.loc[clinical_and_mutation_patients]

  cancer_type = util.get_cancer_type(copy_number)
  outfile = os.path.join(outdir, cancer_type + '.cnv_with_mutation_zscores.csv')
  formatstring = '{0}, {1}, {2}, {3}, {4}, {5}, {6}\n'

  with open(outfile, 'w') as out:
    out.write('gene,mutated zscore,mutated pvalue,mutated patients,non-mutated zscore, non-mutated pvalue, non-mutated patients\n')
    for gene in clinical_and_cnv_with_mutations:
      if gene not in ('time', 'censor'): # skip metadata
        clinical_gene = clinical_and_cnv_with_mutations[[gene, 'time', 'censor']]

        if gene in mutation:
          mutations_for_gene = mutation[gene].rename('mutation')
          with_mutation = clinical_gene.join(mutations_for_gene.dropna(), how='inner')
          without_mutation = clinical_gene.join(mutations_for_gene[mutations_for_gene != 1], how='inner')
        else:
          with_mutation = pd.DataFrame({gene: []})
          without_mutation = clinical_gene

        without_mutation_cox_dict = calculate_cox(without_mutation, gene)
        with_mutation_cox_dict = calculate_cox(with_mutation, gene)
        out.write(formatstring.format(gene,
          with_mutation_cox_dict['z'], with_mutation_cox_dict['p'], with_mutation_cox_dict['n'],
          without_mutation_cox_dict['z'], without_mutation_cox_dict['p'], without_mutation_cox_dict['n']))

def multiprocess_zscores(args):
  copy_number = args[0]
  mutation = args[1]
  clinical = args[2]
  outdir = args[3]
  make_zscores(copy_number, mutation, clinical, outdir)

def main(argv=None):
  cnv_dir, mutation_dir, clinical_dir, outdir = get_options()
  cnv_files = os.listdir(cnv_dir)
  cnv_files = util.remove_extraneous_files(cnv_files)
  cnv_files = [os.path.join(cnv_dir, i) for i in cnv_files]

  zscore_inputs = []
  for cnv in cnv_files:
    cancer_type = util.get_cancer_type(cnv)
    mutation = glob.glob(os.path.join(mutation_dir, cancer_type + '*'))[0]
    clinical = glob.glob(os.path.join(clinical_dir, '*' + cancer_type + '*'))[0]
    zscore_inputs.append([cnv, mutation, clinical, outdir])

  p = Pool(4)
  p.map(multiprocess_zscores, zscore_inputs)

if __name__ == "__main__":
  main()
