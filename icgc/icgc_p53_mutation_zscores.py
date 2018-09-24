#!/usr/bin/env python
# encoding: utf-8
'''
icgc_mutation_zscores.py

Created by Joan Smith
on 2017-6-18

Calculate P53 zscores for each ICGC mutation file

Copyright (c) 2018. All rights reserved.
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

MUTATION_PERCENT = .02

def get_options():
  parser = argparse.ArgumentParser(description='Get univariate cox for P53')
  parser.add_argument('-i', action='store', dest='input_directory')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  ns = parser.parse_args()

  return (ns.input_directory, ns.output_directory)


def main():
  indir, outdir = get_options()

  print os.path.join(indir, 'TP53_data.csv')
  files = glob.glob(os.path.join(indir, '*', 'TP53_data.csv'))

  results = []
  for f in files:
    cancer_type = os.path.basename(os.path.dirname(f))
    print cancer_type
    df = pd.read_csv(f, index_col=0)
    cox_dict = analysis.do_cox(df.Time, df.Censor, df.mutated)
    cox_dict['cancer_type'] = cancer_type
    results.append(cox_dict)

  results_df = pd.DataFrame(results)
  print results_df
  results = results_df.set_index('cancer_type')
  results.to_csv(os.path.join(outdir, 'icgc_p53_mutation_zscores.csv'))


if __name__ == "__main__":
  main()

