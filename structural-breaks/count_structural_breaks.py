#!/usr/bin/env python
# encoding: utf-8
'''
Created by Joan Smith
on 2018-09-08

Copyright (c) 2018. All rights reserved.
'''

import sys
import os
import argparse
import glob

import pandas as pd
import numpy as np

sys.path.append('../common/')
import utilities as util
import analysis

def get_options():
  parser = argparse.ArgumentParser(description='Count strucutral breaks for every cancer type')
  parser.add_argument('-c', action='store', dest='clinical', default='.')
  parser.add_argument('-d', action='store', dest='copy_number_loc', default='.')
  parser.add_argument('-o', action='store', dest='outdir', default='.')
  ns = parser.parse_args()

  return (ns.copy_number_loc, ns.clinical,
          ns.outdir)

def do_count(g):
  return g.groupby('Chromosome')['Segment_Mean'].nunique()

def count_breaks(cn_file):
  cn = pd.read_csv(cn_file, sep='\t')
  breaks = cn.groupby('Sample').apply(do_count).sum(axis=1)
  breaks.name = 'breaks'
  return pd.DataFrame(breaks)


def main():
  copy_number_loc, clinical, outdir = get_options()
  cnas = os.listdir(copy_number_loc)
  cnas = util.remove_extraneous_files(cnas)

  results = pd.DataFrame()
  for c in cnas:
    cancer_type = util.get_cancer_type(c)
    print cancer_type

    clinical_file = glob.glob(os.path.join(clinical, '*' + cancer_type + '*.txt'))[0]
    clin = util.get_clinical_data(clinical_file)

    patient_breaks = count_breaks(os.path.join(copy_number_loc, c))
    patient_breaks = patient_breaks.reset_index()
    patient_breaks = util.maybe_clear_non_01s(patient_breaks, 'Sample', cancer_type)
    patient_breaks = util.add_identifier_column(patient_breaks, 'Sample')
    patient_breaks = patient_breaks.set_index('identifier')
    patient_breaks = patient_breaks.drop('Sample', axis=1)

    breaks_and_clin = patient_breaks.join(clin, how='inner')
    breaks_and_clin.to_csv(os.path.join(outdir, cancer_type + '_breaks.csv'))
    cox = analysis.do_cox(breaks_and_clin.time, breaks_and_clin.censor, breaks_and_clin.breaks)
    cox['cancer_type'] = cancer_type
    results = results.append(cox, ignore_index=True)

  results.to_csv(os.path.join(outdir, 'cox_results.csv'))

if __name__ == '__main__':
  main()



