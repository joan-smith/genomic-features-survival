#!/usr/bin/env python
# encoding: utf-8
'''
Created by Joan Smith
on 2018-5-4.

Produce mutation zscores for only those patients that are NOT hypermutated.

Copyright (c) 2018. All rights reserved.
'''

import pandas as pd
import numpy as np
import argparse
import sys
import os
import pdb
import rpy2
from multiprocessing import Pool

sys.path.append('../common/')
import utilities as util
import analysis
import mutation_base as mb

HEADERS_BY_FILETYPE = {
    'RPPA-data': 'Composite.Element.REF',
    'microRNA': 'miRseqMature',
    }


MUTATION_PERCENT = 0.02

def get_options():
  parser = argparse.ArgumentParser(description='Mutation zscores for patients that are not hypermutated')
  parser.add_argument('-b', action='store', dest='base_dir')
  parser.add_argument('-c', action='store', dest='clinical_dir')
  parser.add_argument('-y', action='store', dest='hypermutated_patients')
  parser.add_argument('-o', action='store', dest='outdir', default='.')

  ns = parser.parse_args()

  return (ns.base_dir, ns.clinical_dir, ns.hypermutated_patients, ns.outdir)

def make_zscores(data, clinical, hypermutated_patients, outdir):
  clinical_data = util.get_clinical_data(clinical)
  hypermutated = set(clinical_data.index).intersection(hypermutated_patients['patients'])
  print 'Hypermutated in clinical file:', len(hypermutated)
  clinical_data = clinical_data.drop(hypermutated)

  cancer_type = util.get_cancer_type(data)
  df = mb.prep_mutation_data(data, clinical_data)

  print 'Remaining hypermutated:', set(df.index).intersection(hypermutated)
  num_patients = len(set(clinical_data.index) & set(df.index))
  print 'Number of patients present in both:', num_patients

  clinical_and_data = df.join(clinical_data, how='inner')
  print 'Num patients, other count:', len(df.index)

  outfile = os.path.join(outdir, cancer_type + '_non-hypermutated_zscores.csv')
  formatstring = '{0}, {1}, {2}, {3}, {4}\n'

  zscore_count = 0
  zscore_skipped = 0
  with open(outfile, 'w') as out:
    out.write('gene,zscore,pvalue,num patients,num mutations\n')
    for gene in clinical_and_data:
      if gene not in ('time', 'censor', 'index'): # skip metadata
        num_mutations = clinical_and_data[gene].sum()
        # print gene, num_mutations
        if num_mutations >= MUTATION_PERCENT * num_patients:
          try:
            cox_dict = analysis.do_cox(clinical_and_data.time,
                                       clinical_and_data.censor,
                                       clinical_and_data[gene])
            out.write(formatstring.format(gene, cox_dict['z'], cox_dict['p'], cox_dict['n'], num_mutations))
            zscore_count += 1
          except rpy2.rinterface.RRuntimeError as e:
            print 'WARN: skipped ', gene, ' due to R error'
            zscore_skipped += 1
            continue
        else:
          zscore_skipped += 1
          continue

def multiprocess_zscores(ins):
  make_zscores(*ins)

def main():
  basedir, clinical_dir, hypermutated_patients, outdir = get_options()

  hypermutated = pd.read_csv(hypermutated_patients, header=None, names=['patients'])

  data_files = os.listdir(basedir)
  data_files = util.remove_extraneous_files(data_files)
  data_files_by_cancer_type = {util.get_cancer_type(f): f for f in data_files}


  clinical_files = os.listdir(clinical_dir)
  clinical_files = util.remove_extraneous_files(clinical_files)
  inputs = []
  for clinical in clinical_files:
    cancer_type = clinical.split('.')[0]
    data_file = data_files_by_cancer_type[cancer_type]

    make_zscores(os.path.join(basedir, data_file),
               os.path.join(clinical_dir, clinical),
               hypermutated,
               outdir)
    # inputs.append((os.path.join(basedir, data_file),
    #             os.path.join(clinical_dir, clinical),
    #             outdir))
    #
  # p = Pool(10)
  # p.map(multiprocess_zscores, inputs)

if __name__ == "__main__":
  main()
