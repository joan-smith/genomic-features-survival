#!/usr/bin/env python
# encoding: utf-8
'''
mskcc_utilities.py

Given an MSKCC clinical file, produce a dictionary of dataframes of clinical files for each cancer type.

Created by Joan Smith
on 2017-9-21


Copyright (c) 2018. All rights reserved.
'''

import argparse
import sys
import os
import glob
import pandas as pd
import numpy as np
import rpy2

sys.path.append('../common/')
import utilities as util
import analysis

INCLUDE_MUTATIONS = ['In_Frame_Ins', 'Nonstop_Mutation', 'Translation_Start_Site', 'In_Frame_Del',
                     'Splice_Region', 'Frame_Shift_Ins', 'Frame_Shift_Del', 'Splice_Site',
                     'Nonsense_Mutation', 'Missense_Mutation']

MUTATION_PERCENT = .02

def get_options():
  parser = argparse.ArgumentParser(description='Collect MSKCC clinical data')
  parser.add_argument('-c', action='store', dest='clinical_file')
  parser.add_argument('-a', action='store', dest='cna_file')
  parser.add_argument('-m', action='store', dest='mutation_file')
  parser.add_argument('-n', action='store', dest='name_conversions_file')
  parser.add_argument('-o', action='store', dest='outdir', default='.')

  namespace = parser.parse_args()
  return (namespace.clinical_file, namespace.cna_file, namespace.mutation_file,
          namespace.name_conversions_file, namespace.outdir)


def do_single_cancer_type_cna(cancer_type, cancer_type_clinical, cna_file, name_conversions, outdir):
  print 'Duplicate count:',  cancer_type_clinical.index.duplicated(keep='first').sum()
  cancer_type_cnas = cna_file[cancer_type_clinical.index].T
  print 'Patient count for CNAs:', cancer_type_cnas.shape
  cancer_type_cnas_and_clinical = cancer_type_cnas.join(cancer_type_clinical)
  print cancer_type_cnas_and_clinical.shape
  num_patients = cancer_type_cnas_and_clinical.shape[0]
  print 'num patients:', num_patients
  formatstring = '{0}, {1}, {2}, {3}\n'

  outfile = os.path.join(outdir, 'cnas', cancer_type.replace(' ', '-') + '.out.csv')
  print outfile
  with open(outfile, 'w') as out:
    out.write('gene,zscore,pvalue,num patients\n')
    for gene in cancer_type_cnas_and_clinical:
        if gene in ['Time', 'Censor']:
          continue
        number_non_zero = cancer_type_cnas_and_clinical[cancer_type_cnas_and_clinical[gene] != 0][gene].shape[0]
        try:
          # print cancer_type_cnas_and_clinical[['Time', 'Censor', gene]]
          cox_dict = analysis.do_cox(cancer_type_cnas_and_clinical.Time,
                                              cancer_type_cnas_and_clinical.Censor,
                                              cancer_type_cnas_and_clinical[gene])
          if gene in name_conversions.index:
            print 'Converting gene', gene, 'to', name_conversions['TCGA'].loc[gene]
            gene = name_conversions['TCGA'].loc[gene]
          out.write(formatstring.format(
                        gene, cox_dict['z'], cox_dict['p'], cox_dict['n']))
        except rpy2.rinterface.RRuntimeError as e:
          print 'Skipped ', gene, 'due to R error.'

def mutations_for_gene(df):
  mutated_patients = df['Tumor_Sample_Barcode'].unique()
  return pd.DataFrame({'mutated': np.ones(len(mutated_patients))}, index=mutated_patients)

def make_mutations_by_patient(mutations):
  df = mutations[mutations['Variant_Classification'].isin(INCLUDE_MUTATIONS)]
  number_barcodes_in_mutation_data = df[u'Tumor_Sample_Barcode'].unique().size
  print 'Number of total sequenced barcodes:', number_barcodes_in_mutation_data
  gene_mutation_df = df.groupby(['Hugo_Symbol']).apply(mutations_for_gene)
  gene_mutation_df.index.set_names(['Hugo_Symbol', 'Tumor_Sample_Barcode'], inplace=True)
  gene_mutation_df = gene_mutation_df.reset_index()
  gene_patient_mutations = gene_mutation_df.pivot(index='Hugo_Symbol', columns='Tumor_Sample_Barcode', values='mutated')
  return gene_patient_mutations.fillna(0)


def do_single_cancer_type_mutation(cancer_type, cancer_type_clinical, mutation_file, name_conversions, outdir):
  patients_in_both = list(set(mutation_file.columns).intersection(set(cancer_type_clinical.index)))
  cancer_type_mutations = mutation_file[patients_in_both].T
  print 'Patient count for Mutations:', cancer_type_mutations.shape
  cancer_type_mutations_and_clinical = cancer_type_mutations.join(cancer_type_clinical)
  formatstring = '{0}, {1}, {2}, {3}, {4}\n'
  num_patients = cancer_type_mutations_and_clinical.shape[0]
  outfile = os.path.join(outdir, 'mutations', cancer_type.replace(' ', '-') + '.out.csv')
  print outfile
  with open(outfile, 'w') as out:
    out.write('gene,zscore,pvalue,num mutations,num patients\n')
    for gene in cancer_type_mutations_and_clinical:
        if gene in ['Time', 'Censor']:
          continue
        if cancer_type_mutations_and_clinical[gene].sum() >= MUTATION_PERCENT*num_patients:
          try:
            cox_dict = analysis.do_cox(cancer_type_mutations_and_clinical.Time,
                                                cancer_type_mutations_and_clinical.Censor,
                                                cancer_type_mutations_and_clinical[gene])
            orig_gene = gene
            if gene in name_conversions.index:
              print 'Converting gene', gene, 'to', name_conversions['TCGA'].loc[gene]
              gene = name_conversions['TCGA'].loc[gene]
            out.write(formatstring.format(
                          gene, cox_dict['z'], cox_dict['p'], cancer_type_mutations_and_clinical[orig_gene].sum(), cox_dict['n']))
            cancer_type_mutations_and_clinical.to_csv(os.path.join(outdir, 'mutations/', cancer_type + '_' +  gene + '_mutations.csv'), columns=['Time', 'Censor', orig_gene])
          except rpy2.rinterface.RRuntimeError as e:
            print 'Skipped ', gene, 'due to R error.'


def do_work(clinical_file, cna_file, mutation_file, name_conversions_file, outdir):
  clinical = pd.read_csv(clinical_file)
  cna_file = pd.read_csv(cna_file, index_col=0)

  mutations = pd.read_excel(mutation_file, header=7)
  mutations_by_patient = make_mutations_by_patient(mutations)

  name_conversions = pd.read_csv(name_conversions_file, index_col=0)

  cancer_types = clinical.groupby('relevant_type')
  for cancer_type, g in cancer_types:
    print cancer_type, 'Patient Count:' , g['Sample ID'].unique().size
    cancer_type_clinical = g.set_index('Sample ID')
    cancer_type_clinical = cancer_type_clinical[['Time', 'Censor']]
    do_single_cancer_type_cna(cancer_type, cancer_type_clinical, cna_file, name_conversions, outdir,)
    do_single_cancer_type_mutation(cancer_type, cancer_type_clinical, mutations_by_patient, name_conversions, outdir)


def main():
  clinical_file, cna_file, mutation_file, name_conversions_file, outdir = get_options()
  do_work(clinical_file, cna_file, mutation_file, name_conversions_file, outdir)

if __name__ == "__main__":
  main()

