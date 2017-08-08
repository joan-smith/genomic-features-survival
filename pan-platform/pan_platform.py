#!/usr/bin/env python
# encoding: utf-8
'''
pan_platform.py

Created by Joan Smith
on 2017-4-28.

Given a directory that has the specified folders, combine pancan data into a single file, with
the row prefixed by the platform (filetype)

Copyright (c) 2017 . All rights reserved.
'''

import pandas as pd
import numpy as np
import getopt
import sys
import os
import glob

PLATFORMS = ['mutation', 'methylation-mean', 'CNV', 'rnaseq', 'rppa', 'microRNA']

def usage():
  print 'o: output directory, i: input directory. uses input directory to look for pancan.csv files'
  sys.exit(1)

def get_options(argv):
  try:
    opts, args = getopt.getopt(argv[1:], 'ho:i:ms:',
    ['help', 'outdir=', 'indir=', 'metagene', 'suffix'])
  except getopt.error, msg:
    usage()

  indir = None
  outdir = None
  metagene = False
  suffix = None

  for option, value in opts:
    print option, value
    if option in ('-o', '--outdir'):
      outdir = value
    if option in ('-h', '--help'):
      usage()
    if option in ('-i', '--indir'):
      indir = value
    if option in ('-m', '--metagene'):
      metagene = True
    if option in ('-s', '--suffix'):
      suffix = value


  if not outdir:
    outdir = indir
  return outdir, indir, metagene, suffix

def make_panplatform(outdir, indir, metagene, suffix):
  if metagene:
    filename = 'metagene_pancan.csv'
    outfile = 'metagene_panplatform.csv'
  else:
    filename = 'pancan.csv'
    outfile = 'panplatform.csv'
  files = []
  for platform in PLATFORMS:
    if platform == 'mutation':
      f = glob.glob(os.path.join(indir, platform + suffix + '-2percent', filename))
    else:
      f = glob.glob(os.path.join(indir, platform + suffix, filename))
    assert(len(f) == 1)
    files.append(f[0])

  panplatform_df = pd.DataFrame()
  for f in files:
    prefix = os.path.split(os.path.split(f)[0])[1]

    df = pd.read_csv(f)
    df['gene'] = prefix + ':' + df['gene']
    df = df.set_index('gene')
    panplatform_df = panplatform_df.append(df, verify_integrity=True)

  panplatform_df.to_csv(os.path.join(outdir, outfile))



def main(argv=None):
  if argv is None:
    argv = sys.argv
    outdir, indir, metagene, suffix = get_options(argv)
    make_panplatform(outdir, indir, metagene, suffix)

if __name__ == "__main__":
  main()

