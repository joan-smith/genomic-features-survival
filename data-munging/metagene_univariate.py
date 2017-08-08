#!/usr/bin/env python
# encoding: utf-8
'''
pcna25_univariate.py

Created by Joan Smith
on 2017-7-22.

Calculate univariate zscores for pcna25 for rnaseq

Copyright (c) 2017 . All rights reserved.
'''

import pandas as pd
import numpy as np
import argparse
import sys
import os
import pdb
import rpy2

sys.path.append('../common/')
import utilities as util
import analysis
import metagene as metagene_lib

def get_options():
  parser = argparse.ArgumentParser(description='Output univariate cox from pcna25')
  parser.add_argument('-r', action='store', dest='rnaseq_directory')
  parser.add_argument('-c', action='store', dest='clinical_directory')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  parser.add_argument('-m', action='store', dest='metagene_file')
  namespace = parser.parse_args()

  return (namespace.rnaseq_directory, namespace.clinical_directory,
    namespace.output_directory, namespace.metagene_file)


def make_zscore(clinical_file, metagene):
  clinical_data = util.get_clinical_data(clinical_file)
  metagene_and_clinical = clinical_data.join(metagene, how='inner')
  cox_dict = analysis.do_cox(metagene_and_clinical.time,
                             metagene_and_clinical.censor,
                             metagene_and_clinical.metagene)
  print cox_dict
  return cox_dict


def main():
  rnaseq_directory, clinical_dir, output_dir, metagene_file = get_options()
  files = os.listdir(clinical_dir)
  files = util.remove_extraneous_files(files)

  outdict = {}
  for f in files:
    cancer_type = util.get_cancer_type(f)
    print cancer_type
    clinical_file = os.path.join(clinical_dir, cancer_type + '.clin.merged.txt')
    metagene_data = metagene_lib.get_metagene_data(metagene_file, cancer_type)
    cox_dict = make_zscore(clinical_file, metagene_data)
    outdict[cancer_type] = cox_dict
  outdf = pd.DataFrame(outdict)
  outdf = outdf.transpose()
  print outdf
  outfile = os.path.join(output_dir, 'univariate_pcna25.csv')
  outdf.to_csv(outfile)


if __name__ == "__main__":
  main()

