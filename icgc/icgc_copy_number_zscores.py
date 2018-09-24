#!/usr/bin/env python
# encoding: utf-8
'''
icgc_copy_number_zscores.py
Forked from copy-number/copy_number_zscores.py


Created by Joan Smith
on 2017-3-18.

Given a clinical file and a cnv for genes file, calculate zscores for all genes.

Copyright (c) 2018. All rights reserved.
'''

import pandas as pd
import numpy as np
import argparse
import sys
import os
import pdb

sys.path.append('../common/')
import utilities as util
import analysis

COPY_NUMBER_PERCENT = 0.02 #unused


def get_options():
  parser = argparse.ArgumentParser(description='Get permutation and thread counts')
  parser.add_argument('-i', action='store', dest='input_directory')
  parser.add_argument('-c', action='store', dest='clinical_directory')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  namespace = parser.parse_args()

  return (namespace.input_directory, namespace.clinical_directory,
    namespace.output_directory)

def make_zscores(copy_number, clinical_data, outdir):
  df = pd.read_csv(copy_number, sep=',')
  df_by_patient = df.transpose()
  df_by_patient.columns = df_by_patient.loc['Symbol']
  df_by_patient = df_by_patient.clip(upper=10)
  num_patients = df_by_patient.shape[0]
  clinical_and_cnv = df_by_patient.join(clinical_data, how='inner')

  cancer_type = util.get_cancer_type(copy_number)

  outfile = os.path.join(outdir, cancer_type + '_zscores.csv')
  formatstring = '{0}, {1}, {2}, {3}\n'

  with open(outfile, 'w') as out:
    out.write('gene,zscore,pvalue,num patients\n')
    for gene in clinical_and_cnv:
      if gene not in ('Time', 'Censor'): # skip metadata
        num_with_copy_number = (clinical_and_cnv[gene] != 0).sum()
        cox_dict = analysis.do_cox(clinical_and_cnv.Time,
                                   clinical_and_cnv.Censor,
                                   clinical_and_cnv[gene])
        out.write(formatstring.format(gene, cox_dict['z'], cox_dict['p'], cox_dict['n']))

def get_icgc_cancer_type(f):
  return f.split('.')[0]

def main(argv=None):
  indir, clinical_dir, outdir = get_options()

  files = os.listdir(indir)
  files = util.remove_extraneous_files(files)
  for copy_number in files:
    cancer_type = get_icgc_cancer_type(copy_number)
    print cancer_type
    clinical_file = os.path.join(clinical_dir, cancer_type + '.csv')

    relevant_clinical = pd.read_csv(clinical_file, index_col=0, low_memory=False)[['Time', 'Censor']].astype(float)
    make_zscores(os.path.join(indir, copy_number), relevant_clinical, outdir)

if __name__ == "__main__":
  main()
