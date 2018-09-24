#!/usr/bin/env python
# encoding: utf-8
'''
tumor_stage_group_counts.py

Created by Joan Smith
on 2017-9-21.

Given a folder with TCGA clinical files, and a folder with tumor stage groups, calculate counts

Copyright (c) 2017. All rights reserved.
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
import tumor_stage_util

def get_options():
  parser = argparse.ArgumentParser(description='Tumor stage group counts')
  parser.add_argument('-c', action='store', dest='clinical_directory')
  parser.add_argument('-i', action='store', dest='tumor_groups')
  namespace = parser.parse_args()

  return (namespace.clinical_directory, namespace.tumor_groups)

def count_tumor_groups(clinical_file, tumor_group_file):
  cancer_type = util.get_cancer_type(clinical_file)
  stage_row = tumor_stage_util.TUMOR_STAGE[cancer_type]
  if stage_row:
    tumor_groups = pd.read_csv(tumor_group_file)
    clinical = util.get_clinical_data(clinical_file, extra_rows=[stage_row], extra_rows_numeric=False)
    clinical[stage_row] = clinical[stage_row].str.strip()

    included_stages = []
    for i, group in tumor_groups.iterrows():
      tg = group.dropna().values
      if len(tg) > 0:
        print ', '.join(tg) +  ': ', \
              clinical[clinical[stage_row].isin(tg)][stage_row].count()
        included_stages.extend(tg)
      excluded_patients = clinical[~clinical[stage_row].isin(included_stages)]
    print 'Excluded:'
    print excluded_patients[stage_row].value_counts()


def main():
  clinical_dir, tumor_groups = get_options()
  files = os.listdir(clinical_dir)
  files = util.remove_extraneous_files(files)

  for f in files:
    cancer_type = util.get_cancer_type(f)
    clinical_file = os.path.join(clinical_dir, f)
    tumor_group_file = os.path.join(tumor_groups, cancer_type + '.csv')
    count_tumor_groups(clinical_file, tumor_group_file)


if __name__ == "__main__":
  main()

