#!/usr/bin/env python
# encoding: utf-8
'''
zscores_for_mutants_with_vaf_median.py

Created by Joan Smith
on 2017-9-03.

Calculate zscores for mutants where patients are only counted as mutated if their vaf is >= median vaf for a gene

Copyright (c) 2018. All rights reserved.
'''

import pandas as pd
import numpy as np
import getopt
import sys
import os

sys.path.append('../common/')
import utilities as util
import analysis

import variant_allele_freq

MUTATION_PERCENT = .04
VARIANT_ALLELE_FREQ_CUTOFF = 'MEDIAN'


def prep_data(mutation, clinical_data, key):
  df = pd.read_csv(mutation, sep='\t', low_memory=False, dtype=str)
  cancer_type = util.get_cancer_type(mutation)

  # remove column headers from combined mutation sheet
  df = df[~df[u'Hugo_Symbol'].str.contains('Hugo_Symbol')]
  df[u'Tumor_Sample_Barcode'] = df[u'Tumor_Sample_Barcode'].str.strip()

  number_barcodes_in_mutation_data = df[u'Tumor_Sample_Barcode'].unique().size
  print 'Number of total sequenced barcodes:   ', number_barcodes_in_mutation_data
  df = util.maybe_clear_non_01s(df, u'Tumor_Sample_Barcode', cancer_type)

  df['VAF'] = variant_allele_freq.calculate_vaf(df, key.loc[cancer_type])

  # Reduce mutation data to patients that also have clinical data
  df = util.add_identifier_column(df, u'Tumor_Sample_Barcode')
  df = df.join(clinical_data, on='identifier', how='inner')
  df.set_index([u'Hugo_Symbol', 'identifier'], inplace=True)

  # symmetrically filter clinical data down to patients that were also sequenced
  unique_patients =  df.index.get_level_values('identifier').unique()
  unique_patients_df = pd.DataFrame(unique_patients, index=unique_patients)
  clinical_data_with_sequenced_patients = clinical_data.join(unique_patients_df, how='inner')
  num_patients = clinical_data_with_sequenced_patients.shape[0]
  print 'Number of patients with sequence and clinical data: ', num_patients
  return df, clinical_data_with_sequenced_patients, num_patients

def calculate_cox(mutation, clinical_data, key, outdir):
  df, clinical_data_with_sequenced_patients, num_patients = prep_data(mutation, clinical_data, key)

  #prep output file
  cancer_type = os.path.basename(mutation).split('_')[0].split('.')[0]
  print cancer_type
  outfile = os.path.join(outdir, (cancer_type + '_mutation-fraction-' + str(MUTATION_PERCENT) +
                                                '_vaf_cutoff-' + str(VARIANT_ALLELE_FREQ_CUTOFF) +'.zscores.out.csv'))
  formatstring = '\'{0}, {1}, {2}, {3}, {4}\n'

  with open(outfile, 'w') as out:
    out.write('gene,zscore,pvalue,num mutations,num patients\n')

    #for every gene, collect the clinical data with the mutation data.
    #  only for non-silent mutations
    patients_with_gene = df.groupby(level=u'Hugo_Symbol')
    for gene, gene_df in patients_with_gene:
      # Remove silent mutations
      non_silent = gene_df.where(gene_df[u'Variant_Classification'] != 'Silent')
      non_silent = non_silent.dropna(subset=[u'Variant_Classification'])
      mutated_patient_list = non_silent.index.get_level_values('identifier').unique()

      num_mutations = len(mutated_patient_list)

      if num_mutations >= MUTATION_PERCENT * num_patients:
        # Get "effectively mutated" patients: those who's VAF >= median
        median_vaf = non_silent['VAF'].median()
        greater_than_median = non_silent[non_silent['VAF'] >= median_vaf]
        effectively_mutated_patients = greater_than_median.index.get_level_values('identifier').unique()
        num_effective_mutations = len(effectively_mutated_patients)

        # take the patients with mutations and without, and build an analysis dataframe with time and censor.
        analysis_data = pd.DataFrame(
            {'mutated': np.ones(num_effective_mutations)},
            index=effectively_mutated_patients)
        analysis_data = analysis_data.join(clinical_data_with_sequenced_patients, how='right')
        analysis_data['mutated'].fillna(0, inplace=True)

        #Do analysis!
        print 'Doing analysis for ', gene, num_mutations
        time = analysis_data['time']
        censor = analysis_data['censor']
        split = analysis_data['mutated']

        name = cancer_type+ '_' + gene
        analysis.do_km(name, time, censor, split, outdir)
        cox_dict = analysis.do_cox(time, censor, split)
        if cox_dict['n'] != len(analysis_data['time']):
          print 'ERROR'
        out.write(formatstring.format(gene, cox_dict['z'], cox_dict['p'], num_mutations,cox_dict['n']))
        analysis_data.to_csv(os.path.join(outdir, name + '_data.csv'),
                             columns=['time', 'censor', 'mutated'])


def usage():
  print 'Provide a mutation file with -m and a clincial file with -c'
  print 'Output directory with -o'
  sys.exit(1)

def get_options(argv):
  try:
    opts, args = getopt.getopt(argv[1:], 'hm:c:o:g:k:',
      ['help', 'mutation=', 'clinical=', 'outdir=', 'key='])
  except getopt.error, msg:
    usage()

  mutation = None
  clinical = None
  outdir = '.'
  key = None

  for option, value in opts:
    if option in ('-o', '--outdir'):
      outdir = value
    if option in ('-h', '--help'):
      usage()
    if option in ('-m', '--mutation'):
      mutation = value
    if option in ('-c', '--clinical'):
      clinical = value
    if option in ('-k', '--key'):
      key = value

  return mutation, clinical, outdir, key


def main(argv=None):
  if argv is None:
    argv = sys.argv
    mutation, clinical, outdir, key_file = get_options(argv)
    key = pd.read_csv(key_file, index_col=0, na_values=['-'])
    key = key.dropna(how='all')

    cancer_type = util.get_cancer_type(mutation)
    if cancer_type in key.index:
      clinical_data = util.get_clinical_data(clinical)
      if not os.path.isdir(outdir):
        os.makedirs(outdir)
      calculate_cox(mutation, clinical_data, key, outdir)


if __name__ == "__main__":
  main()
