#!/usr/bin/env python
# encoding: utf-8
'''
Created by Joan Smith
on 2018-9-16

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



def get_options():
  parser = argparse.ArgumentParser(description='Get permutation and thread counts')
  parser.add_argument('-i', action='store', dest='input_directory')
  parser.add_argument('-c', action='store', dest='clinical_directory')
  parser.add_argument('-g', action='store', dest='gene_file')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  ns = parser.parse_args()

  return (ns.input_directory, ns.clinical_directory, ns.gene_file,
    ns.output_directory)

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


def calculate_zscores_for_file(mutation_file, clinical_data, gene_list, cancer_type):
  df, clinical_data_with_sequenced_patients, num_patients = prep_data(mutation_file,
            clinical_data)
  df = df[df.index.get_level_values(0).isin(gene_list)]

  #for every gene, collect the clinical data with the mutation data.
  patients_with_gene = df.groupby(level=u'Hugo_Symbol')
  results = []
  for gene, gene_df in patients_with_gene:
    mutated_patient_list = gene_df.index.get_level_values('Tumor_Sample_Barcode').unique()
    num_mutations = len(mutated_patient_list)

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
    cox_dict['gene'] = gene
    cox_dict['num_mutations'] = num_mutations
    if cox_dict['n'] != len(analysis_data['Time']):
      print 'ERROR'
    if gene[0] != '\'':
      gene = '\'' + gene
    results.append(cox_dict)
  return results



def get_cbioportal_clinical(clinical):
  clinical_data = pd.read_csv(clinical, sep=util.get_sep_from_filename(clinical))
  clinical_data = clinical_data.set_index('PATIENT_ID')
  relevant_clinical = clinical_data[[u'Time', u'Censor']].astype(float)
  relevant_clinical = relevant_clinical.dropna()

  return relevant_clinical


def main():
  indir, clinical_dir, gene_file, outdir = get_options()
  files = os.listdir(indir)
  files = util.remove_extraneous_files(files)
  print files
  files = [f for f in files if 'metabric' in f]
  f = files[0]
  print files

  gene_list = pd.read_csv(gene_file,header=None)[0].values
  print gene_list

  cancer_type = get_cbioportal_cancer_type(f)
  print cancer_type
  clinical_file = os.path.join(clinical_dir, cancer_type + '_clinical.csv')
  clinical = get_cbioportal_clinical(clinical_file)
  results = calculate_zscores_for_file(os.path.join(indir, f),
                                        clinical, gene_list, cancer_type)
  results_df = pd.DataFrame(results)
  results_df = results_df.set_index('gene')
  results_df.to_csv(os.path.join(outdir, 'metabric_mutations.csv'))


if __name__ == "__main__":
  main()

