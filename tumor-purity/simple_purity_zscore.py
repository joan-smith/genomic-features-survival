#!/usr/bin/env python
# encoding: utf-8
'''
quick_zscores.py

Created by Joan Smith
on 2017-7-29

Calculate simple tumor purity zscores for all cancer types.

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

def get_options():
  parser = argparse.ArgumentParser(description='Calculate simple tumor purity zscore for all cancer types')
  parser.add_argument('-c', action='store', dest='clinical_dir', required=True)
  parser.add_argument('-o', action='store', dest='output_dir', default='.')
  parser.add_argument('-d', action='store', dest='header_file')
  namespace = parser.parse_args()

  return (namespace.clinical_dir, namespace.output_dir, namespace.header_file)


def main():
  clinical_dir, output_dir, header_file = get_options()
  headers = pd.read_csv(header_file, index_col=0, header=None)
  clinical_files = os.listdir(clinical_dir)
  clinical_files = util.remove_extraneous_files(clinical_files)

  zscore_data = {}
  for f in clinical_files:
    clinical_path = os.path.join(clinical_dir, f)
    cancer_type = util.get_cancer_type(f)
    purity_header = headers.get_value(cancer_type, 1)

    clinical = util.get_clinical_data(clinical_path, extra_rows=[purity_header])

    cox_dict = analysis.do_cox(clinical.time, clinical.censor, clinical[purity_header])
    zscore_data[cancer_type] = cox_dict

    purity_mean = clinical[purity_header].mean()
    print purity_mean
    clinical['km_split'] = np.where(clinical[purity_header] <= purity_mean, 0, 1)
    analysis.do_km(cancer_type, clinical.time, clinical.censor, clinical.km_split, output_dir)

  out_df = pd.DataFrame(zscore_data).transpose()
  out_df.to_csv(os.path.join(output_dir, 'simple_purity_zscores.csv'))





if __name__ == "__main__":
  main()

