#!/usr/bin/env python
# encoding: utf-8
'''
zscores.py

Created by Joan Smith
on 2018-5-3.

Given a clinical file and a data for genes file, calculate zscores for all genes.

Copyright (c) 2018. All rights reserved.
'''

import pandas as pd
import numpy as np
import argparse
import sys
import os
import pdb
import rpy2
from  multiprocessing import Pool

sys.path.append('../common/')
import utilities as util
import analysis

HEADERS_BY_FILETYPE = {
    'RPPA-data': 'Composite.Element.REF',
    'microRNA': 'miRseqMature',
    }


def get_options():
  parser = argparse.ArgumentParser(description='Make histological type zscores')
  parser.add_argument('-b', action='store', dest='base_dir')
  parser.add_argument('-c', action='store', dest='clinical_dir')
  parser.add_argument('-o', action='store', dest='outdir', default='.')

  ns = parser.parse_args()

  return (ns.base_dir, ns.clinical_dir, ns.outdir)

def prep_data(data):
  df = pd.read_csv(data)
  df_by_patient = df.transpose()
  df_by_patient.columns = df_by_patient.loc['Symbol']
  return df_by_patient

def make_zscores(data, clinical, outdir):
  subtype = clinical.split('.')[1]
  print clinical
  clinical_data = pd.read_csv(clinical, index_col=0, header=0)
  print clinical_data
  clinical_data = clinical_data.dropna(subset=['time', 'censor'], how='any')
  subtype_col = clinical_data.columns[-1]
  print subtype_col

  cancer_type = util.get_cancer_type(data)
  df = prep_data(data)
  print df

  print cancer_type
  print 'Number of patients present in both:', len(set(clinical_data.index) & set(df.index))

  clinical_and_data = df.join(clinical_data, how='inner')

  outfile = os.path.join(outdir, cancer_type + '_' + subtype + '_zscores.csv')
  formatstring = '{0}, {1}, {2}, {3}\n'

  zscore_count = 0
  zscore_skipped = 0
  with open(outfile, 'w') as out:
    out.write('gene,zscore,pvalue,num patients\n')
    for gene in clinical_and_data:
      if gene not in ('time', 'censor', 'index', subtype_col): # skip metadata
        if clinical_and_data[gene].count() <= 10:
          zscore_skipped += 1
          continue
        try:
          cox_dict = analysis.do_cox(clinical_and_data.time,
                                     clinical_and_data.censor,
                                     clinical_and_data[gene])
          out.write(formatstring.format(gene, cox_dict['z'], cox_dict['p'], cox_dict['n']))
          zscore_count += 1
        except rpy2.rinterface.RRuntimeError as e:
          print 'WARN: skipped ', gene, ' due to R error'
          zscore_skipped += 1
          continue

    print 'Total:', clinical_and_data.shape[1] - 3 # minus time, censor, index
    print 'Output length:', zscore_count
    print 'Skipped:', zscore_skipped


def multiprocess_zscores(input):
  make_zscores(*input)

def main():
  basedir, clinical_dir, outdir = get_options()

  data_files = os.listdir(basedir)
  data_files = util.remove_extraneous_files(data_files)
  data_files_by_cancer_type = {util.get_cancer_type(f): f for f in data_files}

  clinical_files = os.listdir(clinical_dir)
  clinical_files = util.remove_extraneous_files(clinical_files)

  inputs = []
  for clinical in clinical_files:
    cancer_type = clinical.split('.')[0]
    data_file = data_files_by_cancer_type[cancer_type]

    inputs.append((os.path.join(basedir, data_file),
                os.path.join(clinical_dir, clinical),
                outdir))
     # make_zscores(os.path.join(basedir, data_file),
     #            os.path.join(clinical_dir, clinical),
     #            outdir)

  p = Pool(4)

  p.map(multiprocess_zscores, inputs)

if __name__ == "__main__":
  main()
