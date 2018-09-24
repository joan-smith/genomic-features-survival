#!/usr/bin/env python
# encoding: utf-8
'''
zscores_for_copy_number.py


Created by Joan Smith
on 2017-3-18.

Given a clinical file and a cnv for genes file, calculate zscores for all genes, with multivariate extra clinical row.

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
  parser.add_argument('-e', action='store', dest='extra_clinical_rows', help='File providing extra clinical rows for each cancer type')

  namespace = parser.parse_args()

  return (namespace.input_directory, namespace.clinical_directory, namespace.output_directory, namespace.extra_clinical_rows)

def make_zscores(copy_number, clinical_data, outdir, extra_clinical_rows=None):
  df = pd.read_csv(copy_number)
  df_by_patient = df.transpose()
  df_by_patient.columns = df_by_patient.loc['Symbol']
  clinical_and_cnv = df_by_patient.join(clinical_data, how='inner')

  cancer_type = util.get_cancer_type(copy_number)
  formatstring = '{0}, {1}, {2}, {3}, {4}, {5}\n'
  outfile = os.path.join(outdir, cancer_type + '_extra_clinical_zscores.csv')

  with open(outfile, 'w') as out:
    out.write('gene,zscore,pvalue,clinical-row-zscore,clinical-row-pvalue,num patients\n')
    for gene in clinical_and_cnv:
      if gene not in ('time', 'censor'): # skip metadata
        if clinical_and_cnv[gene].count() > 10:
          cox_dict = analysis.do_metagene_cox(clinical_and_cnv.time,
                                              clinical_and_cnv.censor,
                                              clinical_and_cnv[gene],
                                              clinical_and_cnv[extra_clinical_rows[0]].rename('metagene'))
          out.write(formatstring.format(
                        gene, cox_dict['z'], cox_dict['p'],
                        cox_dict['metagene-z'], cox_dict['metagene-p'],
                        cox_dict['n']))

def main(argv=None):
  if argv is None:
    argv = sys.argv
    input_directory, clinical, outdir, extra_clinical_rows_file = get_options()
    clinical_files = os.listdir(clinical)
    clinical_files = util.remove_extraneous_files(clinical_files)

    all_extra_clinical_rows = pd.read_csv(extra_clinical_rows_file, index_col=0, header=None)

    for c in clinical_files:
      cancer_type = util.get_cancer_type(c)
      extra_rows = [all_extra_clinical_rows.loc[cancer_type][1]]
      print cancer_type
      clinical_data = util.get_clinical_data(os.path.join(clinical, c),
                                             extra_rows=extra_rows)
      print clinical_data

      copy_number = glob.glob(os.path.join(input_directory, cancer_type + '*.csv'))[0]
      print copy_number

      make_zscores(copy_number, clinical_data, outdir, extra_rows)


if __name__ == "__main__":
  main()
