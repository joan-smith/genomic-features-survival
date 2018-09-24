#!/usr/bin/env python
# encoding: utf-8
'''
cbioportal_mutation_zscores.py

Created by Joan Smith
on 2017-9-20

Calculate zscores for each CBioPortal mutation file

Copyright (c) 2018. All rights reserved.
'''

import argparse
import sys
import os
import glob
import pandas as pd
import numpy as np
import cbioportal_util

sys.path.append('../common/')
import utilities as util
import analysis



MUTATION_PERCENT = .02

def get_options():
  parser = argparse.ArgumentParser(description='Get permutation and thread counts')
  parser.add_argument('-i', action='store', dest='input_directory')
  parser.add_argument('-c', action='store', dest='clinical_directory')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  namespace = parser.parse_args()

  return (namespace.input_directory, namespace.clinical_directory,
    namespace.output_directory)

def get_cbioportal_cancer_type(f):
  return '_'.join(f.split('_')[:2]).upper()


def prep_data(mutation_file, clinical):
  mutation_data = pd.read_csv(mutation_file, sep='\t', comment='#', low_memory=False)
  mutation_data = mutation_data[
      mutation_data[u'Variant_Classification'].isin(cbioportal_util.INCLUDED_MUTATIONS)
    ]

  number_patients_in_mutation_data = mutation_data[u'Tumor_Sample_Barcode'].unique().size
  print 'Number of total sequenced patients:   ', number_patients_in_mutation_data

  # Reduce mutation data to patients that also have clinical data
  df = mutation_data.join(clinical, on='Tumor_Sample_Barcode', how='inner')

  df.set_index([u'Hugo_Symbol', u'Tumor_Sample_Barcode'], inplace=True)

  # symmetrically filter clinical data down to patients that were also sequenced
  unique_patients =  df.index.get_level_values('Tumor_Sample_Barcode').unique()
  unique_patients_df = pd.DataFrame(unique_patients, index=unique_patients)
  clinical_data_with_sequenced_patients = clinical.join(unique_patients_df, how='inner')
  num_patients = clinical_data_with_sequenced_patients.shape[0]
  print 'Number of patients with sequence and clinical data: ', num_patients

  return df, clinical_data_with_sequenced_patients, num_patients


def calculate_zscores_for_file(mutation_file, clinical_data, outdir, cancer_type):
  df, clinical_data_with_sequenced_patients, num_patients = prep_data(mutation_file,
            clinical_data)

  formatstring = '{0}, {1}, {2}, {3}, {4}\n'
  outfile = os.path.join(outdir, cancer_type + '_mutation_percent_'+ str(MUTATION_PERCENT) +  '.cbioportal_zscores.out.csv')
  with open(outfile, 'w') as out:
    out.write('gene,zscore,pvalue,num mutations,num patients\n')

    #for every gene, collect the clinical data with the mutation data.
    patients_with_gene = df.groupby(level=u'Hugo_Symbol')
    for gene, gene_df in patients_with_gene:
      mutated_patient_list = gene_df.index.get_level_values('Tumor_Sample_Barcode').unique()
      num_mutations = len(mutated_patient_list)

      if num_mutations >= MUTATION_PERCENT * num_patients:
        # take the patients with mutations and without, and build an analysis dataframe with time and censor.
        analysis_data = pd.DataFrame(
            {'mutated': np.ones(num_mutations)}, index=mutated_patient_list)
        analysis_data = analysis_data.join(clinical_data_with_sequenced_patients, how='right')
        analysis_data['mutated'].fillna(0, inplace=True)

        #Do analysis!
        print 'Doing analysis for %s: mutated %d of %d' % (gene, num_mutations, num_patients)
        time = analysis_data['Time']
        censor = analysis_data['Censor']
        split = analysis_data['mutated']

        cox_dict = analysis.do_cox(time, censor, split)
        if cox_dict['n'] != len(analysis_data['Time']):
          print 'ERROR'
        if gene[0] != '\'':
          gene = '\'' + gene
        out.write(formatstring.format(gene, cox_dict['z'], cox_dict['p'], num_mutations,cox_dict['n']))
        analysis_data.to_csv(os.path.join(outdir, gene[1:] + '_data.csv'),
                             columns=['Time', 'Censor', 'mutated'], index_label='patient')


def get_cbioportal_clinical(clinical):
  clinical_data = pd.read_csv(clinical, sep=util.get_sep_from_filename(clinical))
  clinical_data = clinical_data.set_index('PATIENT_ID')
  relevant_clinical = clinical_data[[u'Time', u'Censor']].astype(float)
  relevant_clinical = relevant_clinical.dropna()

  return relevant_clinical


def main():
  indir, clinical_dir, outdir = get_options()
  files = os.listdir(indir)
  files = util.remove_extraneous_files(files)

  for f in files:
    cancer_type = get_cbioportal_cancer_type(f)
    print cancer_type
    clinical_file = os.path.join(clinical_dir, cancer_type + '_clinical.csv')
    cancer_type_outdir = os.path.join(outdir, cancer_type)
    if not os.path.isdir(cancer_type_outdir):
      os.makedirs(cancer_type_outdir)
      clinical = get_cbioportal_clinical(clinical_file)
      calculate_zscores_for_file(os.path.join(indir, f), clinical, cancer_type_outdir, cancer_type)


if __name__ == "__main__":
  main()

