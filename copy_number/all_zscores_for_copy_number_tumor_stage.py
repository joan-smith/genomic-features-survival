#!/usr/bin/env python
# encoding: utf-8
'''
zscores_for_copy_number_extra_external_data.py


Created by Joan Smith
on 2017-9-5.
Given a clinical file and a cnv for genes file, calculate zscores for all genes, with multivariate extra data from another file

Copyright (c) 2018. All rights reserved.  '''

import pandas as pd
import numpy as np
import argparse
import sys
import os
import pdb
import glob
from multiprocessing import Pool

sys.path.append('../common/')
import utilities as util
import analysis
import tumor_stage


def get_options():
  parser = argparse.ArgumentParser(description='Run all cancer type zscores for platform')
  parser.add_argument('-i', action='store', dest='input_directory', default='.')
  parser.add_argument('-c', action='store', dest='clinical_directory', default='.')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  parser.add_argument('-e', action='store', dest='extra_data_dir')

  namespace = parser.parse_args()

  return (namespace.input_directory, namespace.clinical_directory, namespace.output_directory, namespace.extra_data_dir)


def make_zscores(copy_number, clinical_data, tumor_stage_data_dir, outdir):
  cancer_type = util.get_cancer_type(copy_number)

  df = pd.read_csv(copy_number)
  df_by_patient = df.transpose()
  df_by_patient.columns = df_by_patient.loc['Symbol']
  clinical_and_cnv = df_by_patient.join(clinical_data, how='inner')

  tumor_stage_data, tumor_stage_cols = tumor_stage.prep_tumor_stage_data(tumor_stage_data_dir, cancer_type)
  if tumor_stage_data is None:
    return

  clinical_and_cnv_and_extra = clinical_and_cnv.join(tumor_stage_data[tumor_stage_cols], how='inner')

  outfile = os.path.join(outdir, cancer_type + '_extra_clinical_zscores.csv')
  header, formatstring = tumor_stage.tumor_stage_output_header_and_format(
            4, tumor_stage_cols)

  with open(outfile, 'w') as out:
    out.write('gene,zscore,pvalue,num patients')
    out.write(header)
    out.write('\n')
    for gene in clinical_and_cnv_and_extra:
      if gene in ['time', 'censor'] + tumor_stage_cols: # skip metadata
        continue
      if clinical_and_cnv_and_extra[gene].count() > 10:
        cox_dict = analysis.do_multivariate_cox(clinical_and_cnv_and_extra.time,
                                            clinical_and_cnv_and_extra.censor,
                                            clinical_and_cnv_and_extra[gene],
                                            clinical_and_cnv_and_extra[tumor_stage_cols])
        group_zscores = tumor_stage.zscores_for_tumor_stage_cols(
                                      cox_dict, tumor_stage_cols)
        out.write(formatstring.format(
                      gene, cox_dict['var-z'], cox_dict['var-p'], cox_dict['var-n'],
                      *group_zscores))


def multiprocess_zscores(args):
  copy_number = args[0]
  clinical_data = args[1]
  cancer_type_outdir = args[3]
  tumor_stage_data = args[2]
  make_zscores(copy_number, clinical_data, tumor_stage_data, cancer_type_outdir)



def main(argv=None):
  if argv is None:
    argv = sys.argv
    input_directory, clinical, outdir, extra_data_dir = get_options()
    clinical_files = os.listdir(clinical)
    clinical_files = util.remove_extraneous_files(clinical_files)

    args = []
    for c in clinical_files:
      cancer_type = util.get_cancer_type(c)
      print cancer_type

      clinical_data = util.get_clinical_data(os.path.join(clinical, c))
      copy_number = glob.glob(os.path.join(input_directory, cancer_type + '*.csv'))[0]

      args.append((copy_number, clinical_data, extra_data_dir, outdir))
      # make_zscores(copy_number, clinical_data, extra_data_dir, outdir)
    p = Pool(4)
    p.map(multiprocess_zscores, args)



if __name__ == "__main__":
  main()
