#!/usr/bin/env python
# encoding: utf-8
'''
quick_zscores.py

Created by Joan Smith
on 2017-7-29

Given a file with patient data (time, censor) and some variable, make zscores

Copyright (c) 2017 . All rights reserved.
'''

import argparse
import sys
import os
import glob
import pandas as pd
import numpy as np

sys.path.append('../common/')
import utilities as util
import analysis

def get_options():
  parser = argparse.ArgumentParser(description='Produce intermediate file for analysis.')
  parser.add_argument('-i', action='store', dest='input_file')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  parser.add_argument('-d', action='store', dest='header')
  namespace = parser.parse_args()

  return (namespace.input_file, namespace.output_directory, namespace.header)


def main():
  infile, outdir, header = get_options()
  data = pd.read_csv(infile, index_col=0)
  cox_dict = analysis.do_cox(data.time, data.censor, data[header])
  print cox_dict
  output = pd.DataFrame(cox_dict, index=[header])
  output.to_csv(os.path.join(outdir, os.path.basename(infile).split('.')[0] + '.single_zscore_output.csv'))


if __name__ == "__main__":
  main()

