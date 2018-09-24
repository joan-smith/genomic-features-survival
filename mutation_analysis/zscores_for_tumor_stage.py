#!/usr/bin/env python
# encoding: utf-8
'''
zscores_for_tumor_stage.py

Created by Joan Smith
on 2017-9-23.

Given a clinical file and a mutations file, and a tumor stage file calculate zscores for genes
with mutations in >x% of patients

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
import tumor_stage

MUTATION_PERCENT = .02
def get_options():
  parser = argparse.ArgumentParser(description='Get mutation and clinical dir')
  parser.add_argument('-i', action='store', dest='mutation_directory')
  parser.add_argument('-c', action='store', dest='clinical_directory')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  parser.add_argument('-t', action='store', dest='tumor_stage_dir', default='.')
  namespace = parser.parse_args()

  return (namespace.mutation_directory, namespace.clinical_directory,
          namespace.output_directory, namespace.tumor_stage_dir)

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

def calculate_cox(mutation, clinical_data, tumor_stage_file, outdir):
  df = prep_data(mutation, clinical_data)
  df = df.join(clinical_data, how='inner')

  tumor_stage_data = pd.read_csv(tumor_stage_file, index_col=0)
  tumor_stage_cols = [i for i in tumor_stage_data if 'group' in i]
  df = df.join(tumor_stage_data[tumor_stage_cols], how='inner')
  num_patients = len(df.index)

  #prep output file
  cancer_type = os.path.basename(mutation).split('_')[0].split('.')[0]
  outfile = os.path.join(outdir, cancer_type + '_mutation-fraction-' + str(MUTATION_PERCENT) + '.zscores.out.csv')
  header, formatstring = tumor_stage.tumor_stage_output_header_and_format(5, tumor_stage_cols)

  with open(outfile, 'w') as out:
    out.write('gene,zscore,pvalue,num mutations,num patients')
    out.write(header)
    out.write('\n')

    for gene in df:
      if gene in ['time', 'censor'] + tumor_stage_cols:
        continue
      num_mutations = df[gene].sum()
      if num_mutations >= MUTATION_PERCENT * num_patients:
        analysis_data = pd.DataFrame()
        analysis_data['time'] = df['time']
        analysis_data['censor'] = df['censor']
        analysis_data['mutated'] = df[gene].fillna(0)
        analysis_data[tumor_stage_cols] = df[tumor_stage_cols]

        #Do analysis!
        cox_dict = analysis.do_multivariate_cox(analysis_data['time'], analysis_data['censor'],
                                                analysis_data['mutated'], analysis_data[tumor_stage_cols])
        tumor_stage_zscores = tumor_stage.zscores_for_tumor_stage_cols(cox_dict, tumor_stage_cols)


        out.write(formatstring.format(gene,
                                      cox_dict['var-z'],
                                      cox_dict['var-p'], num_mutations, cox_dict['var-n'],
                                      *tumor_stage_zscores))
        analysis_data.to_csv(os.path.join(outdir, gene[1:] + '_data.csv'),
                             columns=['time', 'censor', 'mutated'] + tumor_stage_cols)

def multiprocess_zscores(args):
  mutation = args[0]
  clinical_data = args[1]
  tumor_stage = args[2]
  cancer_type_outdir = args[3]
  calculate_cox(mutation, clinical_data, tumor_stage, cancer_type_outdir)


def main(argv=None):
  mutation_dir, clinical_dir, outdir, tumor_stage_dir = get_options()

  clinical_files = os.listdir(clinical_dir)
  clinical_files = util.remove_extraneous_files(clinical_files)
  clinical_files = [os.path.join(clinical_dir, f) for f in clinical_files]
  p = Pool(16)

  args = []
  for clinical in clinical_files:
    cancer_type = util.get_cancer_type(clinical)
    print cancer_type
    mutation = glob.glob(os.path.join(mutation_dir, '*' + cancer_type + '*'))[0]

    tumor_stage = os.path.join(tumor_stage_dir, cancer_type + '_clinical.csv')
    if not os.path.isfile(tumor_stage):
      continue

    clinical_data = util.get_clinical_data(clinical)
    cancer_type_outdir = os.path.join(outdir, cancer_type)
    if not os.path.isdir(cancer_type_outdir):
      os.makedirs(cancer_type_outdir)
    args.append((mutation, clinical_data, tumor_stage, cancer_type_outdir))
    # calculate_cox(mutation, clinical_data, tumor_stage, cancer_type_outdir)

  p.map(multiprocess_zscores, args)


if __name__ == "__main__":
  main()
