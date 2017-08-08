#!/usr/bin/env python
# encoding: utf-8
'''
all_zscores.py

Created by Joan Smith
on 2017-4-1.

Run zscores.py for all cancer types, for a single data type

Copyright (c) 2017 . All rights reserved.
'''

import getopt
import sys
import os
import glob

import zscores

from multiprocessing import Pool

sys.path.append('../common/')
import utilities as util

HEADERS_BY_FILETYPE = {
    'RPPA-data': 'Composite.Element.REF',
    'microRNA': 'miRseqMature',
    }

def get_options(argv):
  try:
    opts, args = getopt.getopt(argv[1:], 'ho:f:d:m:',
    ['help', 'outdir=', 'delete_rows=', 'filetype=', 'metagene='])
  except getopt.error, msg:
    usage()

  delete_rows = []
  filetype = None
  outdir = '.'
  metagene = None

  for option, value in opts:
    if option in ('-o', '--outdir'):
      outdir = value
    if option in ('-h', '--help'):
      usage()
    if option in ('-d', '--delete_rows'):
      delete_rows = [int(i) for i in value.split(',')]
    if option in ('-f', '--filetype'):
      filetype = value
    if option in ('-m', '--metagene'):
      metagene = value

  return filetype, delete_rows, outdir, metagene

def multiprocess_zscores(args):
  data_file = args[0]
  clinical_file = args[1]
  outdir = args[2]
  delete_rows = args[3]
  header = args[4]
  metagene_file = args[5]
  p = zscores.make_zscores(data_file, clinical_file, outdir, delete_rows, header=header, metagene_file=metagene_file)


def main(argv=None):
  if argv is None:
    argv = sys.argv
    filetype, delete_rows, outdir, metagene_file = get_options(argv)
    input_directory = os.path.join('.', filetype)
    files = os.listdir(input_directory)
    files = util.remove_extraneous_files(files)

    header = 'Hybridization REF'
    if filetype in HEADERS_BY_FILETYPE:
      header = HEADERS_BY_FILETYPE[filetype]
    print header

    data_files = []
    for f in files:
      cancer_type = util.get_cancer_type(f)
      clinical_file = os.path.join('.', 'clinical', cancer_type + '.clin.merged.txt')
      data_file = os.path.join('.', filetype, f)
      data_files.append((data_file, clinical_file, outdir, delete_rows, header, metagene_file))

    p = Pool(16)
    p.map(multiprocess_zscores, data_files)

if __name__ == "__main__":
  main()
