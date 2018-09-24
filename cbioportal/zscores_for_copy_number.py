#!/usr/bin/env python
# encoding: utf-8
'''
zscores_for_copy_number.py


Created by Joan Smith
on 2017-3-18.

Given a CBioportal clinical file and a cnv for genes file, calculate zscores for all genes.

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

COPY_NUMBER_PERCENT = 0.02 #unused

def get_options():
  parser = argparse.ArgumentParser(description='Copy numer for cbioportal data')
  parser.add_argument('-i', action='store', dest='input_directory')
  parser.add_argument('-c', action='store', dest='clinical_directory')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  namespace = parser.parse_args()

  return (namespace.input_directory, namespace.clinical_directory,
    namespace.output_directory)

def make_zscores(copy_number, clinical, outdir):
  clinical_data = pd.read_csv(clinical, sep=util.get_sep_from_filename(clinical))
  clinical_data = clinical_data.set_index('PATIENT_ID')
  relevant_clinical = clinical_data[[u'Time', u'Censor']].astype(float)
  relevant_clinical = relevant_clinical.dropna()

  df = pd.read_csv(copy_number, sep=util.get_sep_from_filename(copy_number))

  df = df.drop_duplicates(subset=['Hugo_Symbol'], keep='first')
  df = df.dropna(subset=['Hugo_Symbol'])

  df_by_patient = df.transpose()
  df_by_patient.columns = df_by_patient.loc['Hugo_Symbol']
  clinical_and_cnv = df_by_patient.join(relevant_clinical, how='inner')
  num_patients = clinical_and_cnv.shape[0]

  cancer_type = util.get_cancer_type(copy_number)
  outfile = os.path.join(outdir, cancer_type + '.cbioportal_zscores.csv')
  formatstring = '{0}, {1}, {2}, {3}\n'

  with open(outfile, 'w') as out:
    out.write('gene,zscore,pvalue,num patients\n')
    for gene in clinical_and_cnv:
      if gene not in ('Time', 'Censor'): # skip metadata
        if clinical_and_cnv[gene].count() > 10:

          num_with_copy_number = (clinical_and_cnv[gene] != 0).sum()
          cox_dict = analysis.do_cox(clinical_and_cnv.Time,
                                     clinical_and_cnv.Censor,
                                     clinical_and_cnv[gene],
                                     float_time=True)
          if gene[0] != '\'':
            gene = '\'' + gene
          out.write(formatstring.format(gene, cox_dict['z'], cox_dict['p'], cox_dict['n']))

def get_cbioportal_cancer_type(filename):
  return re.split('_CNV\.', filename)[0]

def main(argv=None):
  if argv is None:
    argv = sys.argv
    input_directory, clinical_directory, outdir = get_options()

    cnv_files = os.listdir(input_directory)
    cnv_files = util.remove_extraneous_files(cnv_files)
    for cnv in cnv_files:
      cancer_type = get_cbioportal_cancer_type(cnv)
      print cancer_type
      clinical_file = glob.glob(os.path.join(clinical_directory, '*' + cancer_type + '*'))[0]

      outglob = glob.glob(os.path.join(outdir, cancer_type + '*'))
      if len(outglob) == 0:
        print cancer_type
        make_zscores(os.path.join(input_directory, cnv), clinical_file, outdir)


if __name__ == "__main__":
  main()
