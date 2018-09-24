#!/usr/bin/env python
# encoding: utf-8
'''

Created by Joan Smith
on 2018-9-16.

Copyright (c) 2018. All rights reserved.
'''

import pandas as pd
import numpy as np
import argparse
import sys
import os
import glob
import re

sys.path.append('../common/')
import utilities as util
import analysis

def get_options():
  parser = argparse.ArgumentParser(description='Copy numer for cbioportal data')
  parser.add_argument('-i', action='store', dest='input_directory')
  parser.add_argument('-c', action='store', dest='clinical_directory')
  parser.add_argument('-g', action='store', dest='genes_list')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  ns = parser.parse_args()

  return (ns.input_directory, ns.clinical_directory,
          ns.genes_list, ns.output_directory)

def make_zscores(copy_number, clinical, gene_list):
  clinical_data = pd.read_csv(clinical, sep=util.get_sep_from_filename(clinical))
  clinical_data = clinical_data.set_index('PATIENT_ID')
  relevant_clinical = clinical_data[[u'Time', u'Censor']].astype(float)
  relevant_clinical = relevant_clinical.dropna()

  df = pd.read_csv(copy_number, sep=util.get_sep_from_filename(copy_number))

  df = df.drop_duplicates(subset=['Hugo_Symbol'], keep='first')
  df = df.dropna(subset=['Hugo_Symbol'])

  df_by_patient = df.transpose()
  df_by_patient.columns = df_by_patient.loc['Hugo_Symbol']
  df_by_patient = df_by_patient[gene_list]
  print df_by_patient
  clinical_and_cnv = df_by_patient.join(relevant_clinical, how='inner')
  num_patients = clinical_and_cnv.shape[0]


  cancer_type = util.get_cancer_type(copy_number)

  results = []
  for gene in clinical_and_cnv:
    if gene in ('Time', 'Censor'): # skip metadata
      continue
    if clinical_and_cnv[gene].count() > 10:
      num_with_copy_number = (clinical_and_cnv[gene] != 0).sum()
      cox_dict = analysis.do_cox(clinical_and_cnv.Time,
                                 clinical_and_cnv.Censor,
                                 clinical_and_cnv[gene],
                                 float_time=True)
      cox_dict['gene'] = gene
      results.append(cox_dict)
  return results


def get_cbioportal_cancer_type(filename):
  return re.split('_CNV\.', filename)[0]

def main(argv=None):
  input_directory, clinical_directory, gene_file, outdir = get_options()

  gene_list = pd.read_csv(gene_file, header=None)[0].values

  cnv_files = os.listdir(input_directory)
  cnv_files = util.remove_extraneous_files(cnv_files)
  cnv_file = [f for f in cnv_files if 'METABRIC' in f]
  cnv = cnv_file[0]

  cancer_type = get_cbioportal_cancer_type(cnv)
  clinical_file = glob.glob(os.path.join(clinical_directory, '*' + cancer_type + '*'))[0]

  print cancer_type
  results = make_zscores(os.path.join(input_directory, cnv), clinical_file, gene_list)

  results_df = pd.DataFrame(results)
  results_df = results_df.set_index('gene')
  results_df.to_csv(os.path.join(outdir, 'metabric_copy_number.csv'))


if __name__ == "__main__":
  main()
