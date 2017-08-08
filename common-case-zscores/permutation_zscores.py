#!/usr/bin/env python
# encoding: utf-8
'''
permutation_zscores.py

Created by Joan Smith
on 2017-5-28.

Given a clinical file and a  data for genes file, calculate zscores for all genes, but permute the data randomly n times.

Copyright (c) 2017 . All rights reserved.
'''

import sys
import argparse
import pandas as pd
import numpy as np
import rpy2
import multiprocessing
import os

import zscores
import all_zscores

sys.path.append('../common/')
import utilities as util
import analysis
import metagene as metagene_lib


def get_permutation_options(argv):
  parser = argparse.ArgumentParser(description='Get permutation and thread counts')
  parser.add_argument('-t', action='store', dest='threads', type=int, default=4)
  parser.add_argument('-p', action='store', dest='num_permutations', default=4, type=int)
  parser.add_argument('-c', action='store_true', dest='cnv')
  namespace, leftovers = parser.parse_known_args()

  return namespace.num_permutations, namespace.threads, namespace.cnv, leftovers

def prep_cnv_data(copy_number):
  df = pd.read_csv(copy_number)
  df_by_patient = df.transpose()
  df_by_patient.columns = df_by_patient.loc['Symbol']
  return df_by_patient

def calculate(args):
  np.random.seed() # multiprocessing spawns each worker with the same seed. override.

  df = args[0]
  clinical = args[1]
  shuffled = df.apply(np.random.permutation, axis=1)

  clinical =  clinical_data = util.get_clinical_data(clinical)
  clinical_and_data = shuffled.join(clinical_data, how = 'inner')

  calculated_zscores = []
  zscore_skipped = 0
  for gene in clinical_and_data:
    if gene not in ('time', 'censor'): # skip metadata
      if clinical_and_data[gene].count() <= 10:
        zscore_skipped += 1
        continue
      try:
        cox_dict = analysis.do_cox(clinical_and_data.time,
                                   clinical_and_data.censor,
                                   clinical_and_data[gene])
        if not np.isnan(cox_dict['z']):
          calculated_zscores.append(cox_dict['z'])
        else:
          zscore_skipped += 1
      except rpy2.rinterface.RRuntimeError as e:
        print 'WARN: skipped ', gene, ' due to R error'
        zscore_skipped += 1
        continue
  return calculated_zscores

def permutations(data, clinical, outdir, delete_rows, header=u'Hybridization REF',
    pool=None, num_permutations=4, cnv=False):
  cancer_type = util.get_cancer_type(data)
  print cancer_type
  if not cnv:
    df = zscores.prep_data(data, delete_rows, header, cancer_type)
    df = df.drop('index', 1)
  else:
    df = prep_cnv_data(data)

  histogram_data = pool.map(calculate, [(df, clinical)]*num_permutations)
  histogram_data = [i for iteration in histogram_data for i in iteration]
  cancer_type = util.get_cancer_type(data)
  outfile = os.path.join(outdir, cancer_type + '_permuted_zscores.csv')
  np.savetxt(outfile, histogram_data)

def single_permutation(argv=None):
  if argv is None:
    argv = sys.argv
    num_permutations, threads, leftovers = get_permutation_options(argv)
    data, clinical, outdir, delete_rows, header, metagene_file = zscores.get_options(
        ['binary'] + leftovers)
    if metagene_file:
      print "Error: not implemented"
      sys.exit(1)

    print num_permutations
    permutations(data, clinical, outdir, delete_rows, header=unicode(header),
        threads=threads, num_permutations=num_permutations)

def make_pancan_file(outdir, num_permutations, show=True):
  print 'Assembling pancan file'
  files = os.listdir(outdir)
  pancan = pd.DataFrame()
  files = util.remove_extraneous_files(files)
  for f in files:
    cancer_type = util.get_cancer_type(f).split('_')[0]
    print cancer_type
    df = pd.read_csv(os.path.join(outdir, f), names=[cancer_type], header=None, sep=',')
    pancan[cancer_type] = df

  print pancan.columns
  df['stouffer_unweighted'] = analysis.stouffer_unweighted(pancan)
  pancan.to_csv(os.path.join(outdir, 'pancan_permutations.csv'))
  util.make_histogram(df['stouffer_unweighted'].dropna().values, outdir,
                      permutation_count=num_permutations, show=show)


def main(argv=None):
  if argv is None:
    argv = sys.argv
    num_permutations, threads, cnv, leftovers = get_permutation_options(argv)
    filetype, delete_rows, outdir, metagene_file = all_zscores.get_options(['binary'] + leftovers)
    if metagene_file is not None:
      print "Metagene is not currently supported"
      sys.exit(1)

    input_directory = os.path.join('.', filetype)
    files = os.listdir(input_directory)
    files = util.remove_extraneous_files(files)

    header = 'Hybridization REF'
    if filetype in all_zscores.HEADERS_BY_FILETYPE:
      header = all_zscores.HEADERS_BY_FILETYPE[filetype]
    print header
    pool = multiprocessing.Pool(threads)

    data_files = []
    for f in files:
      cancer_type = util.get_cancer_type(f)
      clinical_file = os.path.join('.', 'clinical', cancer_type + '.clin.merged.txt')
      data_file = os.path.join('.', filetype, f)
      outfile =  cancer_type + '_permuted_zscores.csv'
      if outfile  in os.listdir(outdir):
        print 'Alredady done', cancer_type, '. Skipping'
        continue
      permutations(data_file, clinical_file, outdir, delete_rows, header,
          cnv=cnv, num_permutations=num_permutations, pool=pool)

    make_pancan_file(outdir, num_permutations, show=False)


if __name__ == "__main__":
  main()
