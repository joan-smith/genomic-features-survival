#!/usr/bin/env python
# encoding: utf-8
'''
stage_counts.py

Created by Joan Smith
on 2017-9-20

Count patients with staged cancer for each type

Copyright (c) 2018. All rights reserved.
'''

import argparse
import sys
import os
import glob
import pandas as pd

sys.path.append('../common/')
import utilities as util
import analysis
import tumor_stage_util


def get_options():
  parser = argparse.ArgumentParser(description='Produce intermediate file for analysis.')
  parser.add_argument('-i', action='store', dest='input')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  namespace = parser.parse_args()
  return (namespace.input, namespace.output_directory)


def main():
 indir, outdir = get_options()
 clinical_files = os.listdir(indir)
 clinical_files = util.remove_extraneous_files(clinical_files)
 stage_row = 'patient.stage_event.pathologic_stage'

 for clinical_f in clinical_files:
   f = os.path.join(indir, clinical_f)
   cancer_type = util.get_cancer_type(clinical_f)
   stage_row = tumor_stage_util.TUMOR_STAGE[cancer_type]
   if stage_row:
     clinical = util.get_clinical_data(f, extra_rows=[stage_row], extra_rows_numeric=False)
     clinical[stage_row] = clinical[stage_row].str.strip()
     print cancer_type
     print clinical[stage_row].value_counts()

if __name__ == "__main__":
  main()

