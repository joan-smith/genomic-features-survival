#!/usr/bin/env python
# encoding: utf-8
'''
zscores_for_copy_number.py

Created by Joan Smith
on 2017-3-18.

Given a clinical file and a cnv for genes file, calculate zscores for all genes.

Copyright (c) 2017. All rights reserved.
'''

import pandas as pd
import numpy as np
import argparse
import sys
import os
import glob
import multiprocessing

sys.path.append('../common/')
import utilities as util
import analysis
import metagene as metagene_lib

def get_options(argv):
  parser = argparse.ArgumentParser(description='Copy number zscores')
  parser.add_argument('-n', action='store', dest='copy_number')
  parser.add_argument('-c', action='store', dest='clinical_dir')
  parser.add_argument('-m', action='store', dest='metagene', default=None)
  parser.add_argument('-o', action='store', dest='outdir', default='.')
  ns = parser.parse_args()

  return (ns.copy_number, ns.clinical_dir, ns.metagene, ns.outdir)

def make_zscores(copy_number, clinical, outdir, metagene_file=None):
  clinical_data = util.get_clinical_data(clinical)

  df = pd.read_csv(copy_number)
  df = df.drop(['Chromosome', 'Location'], axis=1)
  df_by_patient = df.transpose()
  df_by_patient.columns = df_by_patient.loc['Symbol']
  clinical_and_cnv = df_by_patient.join(clinical_data, how='inner')

  cancer_type = util.get_cancer_type(copy_number)
  if metagene_file:
    formatstring = '{0}, {1}, {2}, {3}, {4}, {5}\n'
    outfile = os.path.join(outdir, cancer_type + '_metagene_zscores.csv')

    print "Processing metagene..."
    metagene = metagene_lib.get_metagene_data(metagene_file, cancer_type)
    print "Complete"
  else:
    outfile = os.path.join(outdir, cancer_type + '_zscores.csv')
    formatstring = '{0}, {1}, {2}, {3}\n'

  with open(outfile, 'w') as out:
    if metagene_file:
      out.write('gene,zscore,pvalue,metagene-zscore,metagene-pvalue,num patients\n')
    else:
      out.write('gene,zscore,pvalue,num patients\n')
    for gene in clinical_and_cnv:
      if gene not in ('time', 'censor'): # skip metadata
        if clinical_and_cnv[gene].count() > 10:
          if metagene_file:
            cox_dict = analysis.do_metagene_cox(clinical_and_cnv.time,
                                                clinical_and_cnv.censor,
                                                clinical_and_cnv[gene],
                                                metagene)
            out.write(formatstring.format(gene, cox_dict['z'], cox_dict['p'], cox_dict['metagene-z'], cox_dict['metagene-p'],
              cox_dict['n']))
          else:
            cox_dict = analysis.do_cox(clinical_and_cnv.time,
                                       clinical_and_cnv.censor,
                                       clinical_and_cnv[gene])
            out.write(formatstring.format(gene, cox_dict['z'], cox_dict['p'], cox_dict['n']))

def multiprocess(args):
  infile, clinical, outdir = args
  make_zscores(infile, clinical, outdir)

def all_cancer_types(copy_number_dir, clinical_dir, outdir, parallel_workers=0):
  copy_number_files = os.listdir(copy_number_dir)
  copy_number_files = util.remove_extraneous_files(copy_number_files)

  args = []
  for c in copy_number_files:
    infile = os.path.join(copy_number_dir, c)
    cancer_type = util.get_cancer_type(infile)
    clinical_file = glob.glob(os.path.join(clinical_dir, '*' + cancer_type + '*'))[0]

    if parallel_workers == 0:
      make_zscores(infile, clinical_file, outdir)
    else:
      args.append((infile, clinical_file, outdir))

  p = multiprocessing.Pool(parallel_workers)
  p.map(multiprocess, args)


def main(argv=None):
  if argv is None:
    argv = sys.argv
    copy_number, clinical,metagene_file, outdir = get_options(argv)
    make_zscores(copy_number, clinical, outdir, metagene_file=metagene_file)


if __name__ == "__main__":
  main()
