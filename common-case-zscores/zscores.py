#!/usr/bin/env python
# encoding: utf-8
'''
zscores.py

Created by Joan Smith
on 2017-4-1.

Given a clinical file and a  data for genes file, calculate zscores for all genes.

Copyright (c) 2017 . All rights reserved.
'''

import pandas as pd
import numpy as np
import getopt
import sys
import os
import pdb
import rpy2

sys.path.append('../common/')
import utilities as util
import analysis
import metagene as metagene_lib

def get_options(argv):
  try:
    opts, args = getopt.getopt(argv[1:], 'hc:d:o:e:r:m:',
    ['help', 'data=', 'clinical=', 'outdir=', 'delete_rows=', 'header=', 'metagene='])
  except getopt.error, msg:
    print msg
    usage()

  data = None
  clinical = None
  outdir = '.'
  delete_rows = []
  header = 'Hybridization REF'
  metagene_file = None

  for option, value in opts:
    if option in ('-o', '--outdir'):
      outdir = value
    if option in ('-h', '--help'):
      usage()
    if option in ('-e', '--delete_rows'):
      delete_rows = [int(i) for i in value.split(',')]
    if option in ('-d', '--data'):
      data = value
    if option in ('-c', '--clinical'):
      clinical = value
    if option in ('-r', '--header'):
      header = value
    if option in ('-m', '--metagene'):
      metagene_file = value

  return data, clinical, outdir, delete_rows, header, metagene_file

def prep_data(data, delete_rows, header, cancer_type):
  df = pd.read_csv(data, low_memory=False, sep='\t')
  if delete_rows:
    print 'Deleting rows: ', delete_rows
    df = df.drop(delete_rows)
  df = df.set_index(header)

  df = df.transpose().reset_index()
  df = util.maybe_clear_non_01s(df, 'index', cancer_type)
  df = util.add_identifier_column(df, 'index')
  if df['identifier'].duplicated().sum() > 0:
    print 'ERROR: duplicate patients in input data'
    print df[df['identifier'].duplicated()]['identifier']
    sys.exit(1)

  df = df.set_index('identifier')
  return df

def make_zscores(data, clinical, outdir, delete_rows=None, header=u'Hybridization REF', metagene_file=None):
  clinical_data = util.get_clinical_data(clinical)
  cancer_type = util.get_cancer_type(data)
  df = prep_data(data, delete_rows, header, cancer_type)

  print cancer_type
  print 'Number of patients present in both:', len(set(clinical_data.index) & set(df.index))

  clinical_and_data = df.join(clinical_data, how='inner')

  if metagene_file:
    formatstring = '{0}, {1}, {2}, {3}, {4}, {5}\n'
    outfile = os.path.join(outdir, cancer_type + '_metagene_zscores.csv')

    print "Processing metagene..."
    metagene = metagene_lib.get_metagene_data(metagene_file, cancer_type)
    print "Complete"
  else:
    outfile = os.path.join(outdir, cancer_type + '_zscores.csv')
    formatstring = '{0}, {1}, {2}, {3}\n'

  zscore_count = 0
  zscore_skipped = 0
  with open(outfile, 'w') as out:
    if metagene_file:
      out.write('gene,zscore,pvalue,metagene-zscore,metagene-pvalue,num patients\n')
    else:
      out.write('gene,zscore,pvalue,num patients\n')
    for gene in clinical_and_data:
      if gene not in ('time', 'censor', 'index'): # skip metadata
        if clinical_and_data[gene].count() <= 10:
          zscore_skipped += 1
          continue
        try:
          if metagene_file:
            cox_dict = analysis.do_metagene_cox(clinical_and_data.time,
                                                clinical_and_data.censor,
                                                clinical_and_data[gene],
                                                metagene)
            out.write(formatstring.format(gene, cox_dict['z'], cox_dict['p'], cox_dict['metagene-z'], cox_dict['metagene-p'],
              cox_dict['n']))
          else:
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

def main(argv=None):
  if argv is None:
    argv = sys.argv
    data, clinical, outdir, delete_rows, header, metagene_file = get_options(argv)

    make_zscores(data, clinical, outdir, delete_rows, header=unicode(header), metagene_file=metagene_file)

if __name__ == "__main__":
  main()
