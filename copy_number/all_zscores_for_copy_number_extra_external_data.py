#!/usr/bin/env python
# encoding: utf-8
'''
zscores_for_copy_number_extra_external_data.py


Created by Joan Smith
on 2017-9-5.
Given a clinical file and a cnv for genes file, calculate zscores for all genes, with multivariate extra data from another file

Copyright (c) 2018. All rights reserved.
'''

import pandas as pd
import numpy as np
import argparse
import sys
import os
import pdb
import glob

sys.path.append('../common/')
import utilities as util
import analysis


def get_options():
  parser = argparse.ArgumentParser(description='Run all cancer type zscores for platform')
  parser.add_argument('-i', action='store', dest='input_directory', default='.')
  parser.add_argument('-c', action='store', dest='clinical_directory', default='.')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  parser.add_argument('-e', action='store', dest='extra_data_dir')

  namespace = parser.parse_args()

  return (namespace.input_directory, namespace.clinical_directory, namespace.output_directory, namespace.extra_data_dir)

def make_zscores(copy_number, clinical_data, outdir, extra_data, extra_data_col):
  df = pd.read_csv(copy_number)
  df_by_patient = df.transpose()
  df_by_patient.columns = df_by_patient.loc['Symbol']
  clinical_and_cnv = df_by_patient.join(clinical_data, how='inner')
  clinical_and_cnv_and_extra = clinical_and_cnv.join(extra_data, how='inner')

  cancer_type = util.get_cancer_type(copy_number)
  formatstring = '{0}, {1}, {2}, {3}, {4}, {5}\n'
  outfile = os.path.join(outdir, cancer_type + '_extra_clinical_zscores.csv')

  with open(outfile, 'w') as out:
    out.write('gene,zscore,pvalue,'+extra_data_col+'-zscore,'+extra_data_col+'-pvalue,num patients\n')
    for gene in clinical_and_cnv_and_extra:
      if gene not in ('time', 'censor'): # skip metadata
        if clinical_and_cnv_and_extra[gene].count() > 10:
          cox_dict = analysis.do_metagene_cox(clinical_and_cnv_and_extra.time,
                                              clinical_and_cnv_and_extra.censor,
                                              clinical_and_cnv_and_extra[gene],
                                              clinical_and_cnv_and_extra[extra_data_col].rename('metagene'))
          out.write(formatstring.format(
                        gene, cox_dict['z'], cox_dict['p'],
                        cox_dict['metagene-z'], cox_dict['metagene-p'],
                        cox_dict['n']))

def prep_extra_data(extra_data_directory, cancer_type):
  extra_data_path = os.path.join(extra_data_directory, cancer_type + '.txt')
  extra_data = pd.read_csv(extra_data_path, sep='\t', na_values=['-'])
  extra_data = util.maybe_clear_non_01s(extra_data, 'SampleName', cancer_type)
  extra_data = util.add_identifier_column(extra_data, 'SampleName')
  extra_data = extra_data.drop('SampleName', axis=1)
  extra_data = extra_data.set_index('identifier')
  return extra_data

def main(argv=None):
  if argv is None:
    argv = sys.argv
    input_directory, clinical, outdir, extra_data_dir = get_options()
    clinical_files = os.listdir(clinical)
    clinical_files = util.remove_extraneous_files(clinical_files)
    extra_data_col = 'Purity_InfiniumPurify'

    for c in clinical_files[3:]:
      cancer_type = util.get_cancer_type(c)
      print cancer_type

      if cancer_type == 'COADREAD':
        extra_data = prep_extra_data(extra_data_dir, 'COAD')
      else:
        extra_data = prep_extra_data(extra_data_dir, cancer_type)
      clinical_data = util.get_clinical_data(os.path.join(clinical, c))

      copy_number = glob.glob(os.path.join(input_directory, cancer_type + '*.csv'))[0]

      make_zscores(copy_number, clinical_data, outdir, extra_data, extra_data_col)


if __name__ == "__main__":
  main()
