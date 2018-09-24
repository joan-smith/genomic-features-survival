#!/usr/bin/env python
# encoding: utf-8
'''
tumor_stage_group_counts.py

Created by Joan Smith
on 2017-9-21.

Given a folder with TCGA clinical files, and a folder with tumor stage groups, calculate new clinical data

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
  parser.add_argument('-o', action='store', dest='outdir')
  parser.add_argument('-g', action='store_true', dest='grade')
  namespace = parser.parse_args()

  return (namespace.clinical_directory, namespace.tumor_groups,
          namespace.outdir, namespace.grade)

def make_clinical_data(clinical_file, tumor_group_file, outdir, grade):
  cancer_type = util.get_cancer_type(clinical_file)
  if grade:
    row = tumor_stage_util.TUMOR_GRADE
    if not cancer_type in tumor_stage_util.TUMOR_GRADE_TYPES:
      return
  else:
    row = tumor_stage_util.TUMOR_STAGE[cancer_type]
  if row:
    tumor_groups = pd.read_csv(tumor_group_file)
    tumor_groups = tumor_groups.dropna(how='all')

    clinical = util.get_clinical_data(clinical_file, extra_rows=[row], extra_rows_numeric=False)
    clinical[row] = clinical[row].str.strip()

    included_stages = []
    for i, group in tumor_groups.iterrows():
      tg = group.dropna().values
      if len(tg) > 0:
        print ', '.join(tg) +  ': ', \
              clinical[clinical[row].isin(tg)][row].count()
        included_stages.extend(tg)
        clinical['group_' + str(i)] = np.where(
                            clinical[row].isin(included_stages),
                            0, 1)

    clinical = clinical.drop('group_' + str(i), axis=1)
    clinical = clinical[clinical[row].isin(included_stages)]
    clinical.to_csv(os.path.join(outdir, cancer_type + '_clinical.csv'),
                                 index_label='patient_id')


def main():
  clinical_dir, tumor_groups, outdir, grade = get_options()
  files = os.listdir(clinical_dir)
  files = util.remove_extraneous_files(files)

  for f in files:
    cancer_type = util.get_cancer_type(f)
    clinical_file = os.path.join(clinical_dir, f)
    tumor_group_file = os.path.join(tumor_groups, cancer_type + '.csv')
    make_clinical_data(clinical_file, tumor_group_file, outdir, grade)


if __name__ == "__main__":
  main()

