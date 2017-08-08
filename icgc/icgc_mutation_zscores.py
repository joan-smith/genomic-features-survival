#!/usr/bin/env python
# encoding: utf-8
'''
icgc_mutation_zscores.py

Created by Joan Smith
on 2017-6-18

Calculate zscores for each ICGC mutation file

Copyright (c) 2017 . All rights reserved.
'''

import argparse
import sys
import os
import glob
import pandas as pd
import numpy as np

sys.path.append('../common/')
import utilities as util
import analysis

INCLUDED_MUTATIONS = ['disruptive_inframe_deletion',
          'disruptive_inframe_insertion',
          'frameshift_variant',
          'inframe_deletion',
          'missense_variant',
          'splice_acceptor_variant',
          'splice_donor_variant',
          'stop_gained',
          'stop_lost']
MUTATION_PERCENT = 0.02

def get_options():
  parser = argparse.ArgumentParser(description='Get permutation and thread counts')
  parser.add_argument('-i', action='store', dest='input_directory')
  parser.add_argument('-c', action='store', dest='clinical_directory')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  parser.add_argument('-a', action='store', dest='annotation_file')
  namespace = parser.parse_args()

  return (namespace.input_directory, namespace.clinical_directory,
    namespace.output_directory, namespace.annotation_file)

def get_icgc_cancer_type(f):
  return f.split('.')[2]


def prep_annotations(mutation_ensgs, hgnc):
  # get the ensgs we care about
  unique_ensgs = mutation_ensgs.to_frame().set_index('gene_affected').index.unique()
  unique_ensgs = pd.Series([True], name='mutation_ensgs', index=unique_ensgs)

  # only hold annotations for genes in our mutation file
  annotations = hgnc.join(unique_ensgs, how='right')

  # some ensgs don't have symbols, for these, carry through the ensg
  annotations = annotations.reset_index()
  annotations['Symbol'] = annotations['Symbol'].fillna(annotations['index'])
  annotations = annotations.set_index('index')
  annotations = annotations['Symbol']
  return annotations

def prep_data(mutation_file, clinical_file, hgnc):
  clinical = pd.read_csv(clinical_file, index_col=0)
  relevant_clinical = clinical[[u'Time', u'Censor']].astype(float)

  mutation_data = pd.read_csv(mutation_file, sep='\t')
  mutation_data = mutation_data[mutation_data[u'consequence_type'].isin(INCLUDED_MUTATIONS)]

  number_patients_in_mutation_data = mutation_data[u'icgc_donor_id'].unique().size
  print 'Number of total sequenced patients:   ', number_patients_in_mutation_data

  # Reduce mutation data to patients that also have clinical data
  df = mutation_data.join(relevant_clinical, on='icgc_donor_id', how='inner')

  annotations = prep_annotations(df['gene_affected'], hgnc)
  df['gene_affected'] = df['gene_affected'].map(annotations)

  df.set_index([u'gene_affected', u'icgc_donor_id'], inplace=True)

  # symmetrically filter clinical data down to patients that were also sequenced
  unique_patients =  df.index.get_level_values('icgc_donor_id').unique()
  unique_patients_df = pd.DataFrame(unique_patients, index=unique_patients)
  clinical_data_with_sequenced_patients = relevant_clinical.join(unique_patients_df, how='inner')
  num_patients = clinical_data_with_sequenced_patients.shape[0]
  print 'Number of patients with sequence and clinical data: ', num_patients

  return df, clinical_data_with_sequenced_patients, num_patients


def calculate_zscores_for_file(mutation_file, clinical_file, outdir, hgnc):
  df, clinical_data_with_sequenced_patients, num_patients = prep_data(mutation_file,
            clinical_file, hgnc)

  cancer_type = get_icgc_cancer_type(mutation_file)
  print cancer_type
  formatstring = '{0}, {1}, {2}, {3}, {4}\n'
  outfile = os.path.join(outdir, cancer_type + '_mutation_percent_'+ str(MUTATION_PERCENT) +  '.icgc_zscores.out.csv')
  with open(outfile, 'w') as out:
    out.write('gene,zscore,pvalue,num mutations,num patients\n')

    #for every gene, collect the clinical data with the mutation data.
    patients_with_gene = df.groupby(level=u'gene_affected')
    for gene, gene_df in patients_with_gene:
      mutated_patient_list = gene_df.index.get_level_values('icgc_donor_id').unique()
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



def main():
  indir, clinical_dir, outdir, hgnc_file = get_options()
  files = os.listdir(indir)
  files = util.remove_extraneous_files(files)

  hgnc = pd.read_csv(hgnc_file)
  hgnc = hgnc[['Approved Symbol', 'Ensembl ID(supplied by Ensembl)']]
  hgnc.columns = ['Symbol', 'Ensembl ID']
  hgnc.set_index('Ensembl ID', inplace=True)
  hgnc['Symbol'] = '\'' + hgnc['Symbol']

  for f in files:
    cancer_type = get_icgc_cancer_type(f)
    print cancer_type
    clinical_file = os.path.join(clinical_dir, cancer_type + '.csv')
    cancer_type_outdir = os.path.join(outdir, cancer_type)
    if not os.path.isdir(cancer_type_outdir):
      os.makedirs(cancer_type_outdir)
      calculate_zscores_for_file(os.path.join(indir, f), clinical_file, cancer_type_outdir, hgnc)


if __name__ == "__main__":
  main()

