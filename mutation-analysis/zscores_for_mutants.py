#!/usr/bin/env python
# encoding: utf-8
'''
zscores_for_mutants.py

Created by Joan Smith
on 2016-2-20.

Given a clinical file and a mutations file, calculate zscores for genes with mutations in >10% of patients

Copyright (c) 2017 . All rights reserved.
'''

import pandas as pd
import numpy as np
import getopt
import sys
import os

sys.path.append('../common/')
import utilities as util
import analysis
import metagene as metagene_lib

MUTATION_PERCENT = .02


def prep_data(mutation, clinical_data):
  df = pd.read_csv(mutation, sep='\t', low_memory=False, dtype=str)
  cancer_type = util.get_cancer_type(mutation)

  # remove column headers from combined mutation sheet
  df = df[~df[u'Hugo_Symbol'].str.contains('Hugo_Symbol')]
  df[u'Tumor_Sample_Barcode'] = df[u'Tumor_Sample_Barcode'].str.strip()

  number_barcodes_in_mutation_data = df[u'Tumor_Sample_Barcode'].unique().size
  print 'Number of total sequenced barcodes:   ', number_barcodes_in_mutation_data
  df = util.maybe_clear_non_01s(df, u'Tumor_Sample_Barcode', cancer_type)

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

def calculate_cox(mutation, clinical_data, outdir, metagene_file=None):
  df, clinical_data_with_sequenced_patients, num_patients = prep_data(mutation, clinical_data)

  #prep output file
  cancer_type = os.path.basename(mutation).split('_')[0].split('.')[0]
  print cancer_type
  if metagene_file:
    formatstring = '{0}, {1}, {2}, {3}, {4}, {5}\n'
    outfile = os.path.join(outdir, cancer_type + '_mutation-fraction-' + str(MUTATION_PERCENT) + '_metagene_zscores.csv')

    print "Processing metagene..."
    metagene = metagene_lib.get_metagene_data(metagene_file, cancer_type)
    print "Complete"
  else:
    outfile = os.path.join(outdir, cancer_type + '_mutation-fraction-' + str(MUTATION_PERCENT) + '.zscores.out.csv')
    formatstring = '\'{0}, {1}, {2}, {3}, {4}\n'

  with open(outfile, 'w') as out:
    if metagene_file:
      out.write('gene,zscore,pvalue,metagene-zscore,metagene-pvalue,num patients\n')
    else:
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
        # take the patients with mutations and without, and build an analysis dataframe with time and censor.
        analysis_data = pd.DataFrame(
            {'mutated': np.ones(num_mutations)}, index=mutated_patient_list)
        analysis_data = analysis_data.join(clinical_data_with_sequenced_patients, how='right')
        analysis_data['mutated'].fillna(0, inplace=True)

        #Do analysis!
        print 'Doing analysis for ', gene, num_mutations
        time = analysis_data['time']
        censor = analysis_data['censor']
        split = analysis_data['mutated']

        if metagene_file:
            cox_dict = analysis.do_metagene_cox(time,
                                                censor,
                                                split,
                                                metagene)
            out.write(formatstring.format(gene, cox_dict['z'], cox_dict['p'], cox_dict['metagene-z'],
                                          cox_dict['metagene-p'],
                                          cox_dict['n']))
        else:
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
    opts, args = getopt.getopt(argv[1:], 'hm:c:o:g:',
      ['help', 'mutation=', 'clinical=', 'outdir=', 'metagene='])
  except getopt.error, msg:
    usage()

  mutation = None
  clinical = None
  outdir = '.'
  metagene_file = None

  for option, value in opts:
    if option in ('-o', '--outdir'):
      outdir = value
    if option in ('-h', '--help'):
      usage()
    if option in ('-m', '--mutation'):
      mutation = value
    if option in ('-c', '--clinical'):
      clinical = value
    if option in ('-g', '--metagene'):
      metagene_file = value

  return mutation, clinical, outdir, metagene_file


def main(argv=None):
  if argv is None:
    argv = sys.argv
    mutation, clinical, outdir, metagene_file = get_options(argv)
    clinical_data = util.get_clinical_data(clinical)
    if not os.path.isdir(outdir):
      os.makedirs(outdir)
    calculate_cox(mutation, clinical_data, outdir, metagene_file=metagene_file)


if __name__ == "__main__":
  main()
