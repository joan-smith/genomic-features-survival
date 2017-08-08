#!/usr/bin/env python
# encoding: utf-8
'''
zscores_for_copy_number.py


Created by Joan Smith
on 2017-3-18.

Given a clinical file and a cnv for genes file, calculate zscores for all genes.

Copyright (c) 2017 . All rights reserved.
'''

import pandas as pd
import numpy as np
import getopt
import sys
import os
import pdb

sys.path.append('../common/')
import utilities as util
import analysis
import metagene as metagene_lib


def get_options(argv):
  try:
      opts, args = getopt.getopt(argv[1:], 'hn:c:o:m:',
    ['help', 'copy_number=', 'clinical=', 'outdir=', 'metagene='])
  except getopt.error, msg:
    usage()

  copy_number = None
  clinical = None
  outdir = '.'
  metagene_file = None

  for option, value in opts:
    if option in ('-o', '--outdir'):
      outdir = value
    if option in ('-h', '--help'):
      usage()
    if option in ('-n', '--copy_number'):
      copy_number = value
    if option in ('-c', '--clinical'):
      clinical = value
    if option in ('-m', '--metagene'):
      metagene_file = value

  return copy_number, clinical, outdir, metagene_file

def make_zscores(copy_number, clinical_data, outdir, metagene_file=None):
  df = pd.read_csv(copy_number)
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

def main(argv=None):
  if argv is None:
    argv = sys.argv
    copy_number, clinical, outdir, metagene_file = get_options(argv)
    clinical_data = util.get_clinical_data(clinical)
    make_zscores(copy_number, clinical_data, outdir, metagene_file=metagene_file)


if __name__ == "__main__":
  main()
