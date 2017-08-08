#!/usr/bin/env python
# encoding: utf-8
'''
find-gene-copy-number.py

Created by Joan Smith
on 2016-01-02

Copyright (c) 2016 . All rights reserved.
'''

import sys
import os
import getopt

import pandas as pd
import numpy as np
from intervaltree import Interval, IntervalTree

sys.path.append('../common/')
import utilities as util

def process_input_file(input_file):
  df = pd.read_csv(input_file, sep='\t')
  df = util.maybe_clear_non_01s(df, u'Sample', input_file)
  df = util.add_identifier_column(df, u'Sample')

  patient_data = {}
  for index, row in df.iterrows():
    identifier = row['identifier']
    chromosome = row['Chromosome']
    start = row['Start']
    end = row['End']
    copy_number = row['Segment_Mean']

    # Add the new patient to the data dict, initializing the chromosome list
    if not identifier in patient_data:
      # note the length of the list is 24, so we can use chromosome number (and not worry about index 0)
      patient_data[identifier] = [0] * 24

    # Initialize the interval tree for the new chromosome for this patient
    if patient_data[identifier][chromosome] == 0:
      patient_data[identifier][chromosome] = IntervalTree()

    # Add the range and copy number in this row to the correct patient_data/chromosome location
    # Note the interval tree implementation uses half close intervals, but copy number data
    # uses closed intervals, so we add 1 to the end to ensure our intervaltree matches the data.
    patient_data[identifier][chromosome][start:end+1] = copy_number
  return patient_data

def process_annotation_row(g):
  chromosome = g['chromosome']
  if chromosome == 'X':
    chromosome = '23'
  try:
    chromosome = int(chromosome)
  except ValueError as e:
    # print g
    return None

  return pd.Series({'gene': g['Approved Symbol'], 'start': g.txStart,
    'chromosome': chromosome})

def process_annotation_file(annotation_file):
  # columns: 'Approved Symbol', 'Chromosome' 'txStart' (could use 'txEnd' too)
  # plan of record, use txStart as gene location
  annotation = pd.read_csv(annotation_file, low_memory=False)
  processed_annotation = annotation.apply(process_annotation_row, axis=1)
  return processed_annotation.dropna(how="any")

def get_copy_number_for_gene_start(chromosome_copy_numbers, start, patient, gene):
  intervals = chromosome_copy_numbers[start]
  if len(intervals) > 1:
    print 'Err: more than one interval found for patient:', patient, 'for gene:', gene
    print 'Got:', intervals
    return np.nan
  if len(intervals) == 0:
    # print 'Err: no interval found for patient:', patient, 'for gene:', gene
    return np.nan
  else:
    interval = intervals.pop()
    return interval.data

def process_and_write_data(outfile, annotation_file, patient_data):
  patients = sorted(patient_data.keys())
  # The best way of writing this data is the old fashioned way: line by line.
  # Building up an in-memory data structure and then writing to disk is far too slow for the
  # 100s of MB files we hae here. So: classic it is.
  with open(outfile, 'w') as f:
    f.write('Symbol,Chromosome,Location,'+ ','.join(patients) + '\n')
    for index, row in annotation_file.iterrows():
      if not np.isnan(row.chromosome):
        f.write('\'' + str(row.gene) + ',' + str(row.chromosome) + ',' + str(row.start) + ',')
        for i, patient in enumerate(patients):
          all_chromosome_copy_numbers = patient_data[patient]
          this_chromosome_copy_numbers = all_chromosome_copy_numbers[int(row.chromosome)]
          copy_number = get_copy_number_for_gene_start(
              this_chromosome_copy_numbers, row.start, patient, row.gene)
          f.write(str(copy_number))
          #  we do need to make sure there's no trailing comma though:
          if i != len(patients) - 1:
            f.write(',')
        f.write('\n')

def get_options(argv):
  try:
    opts, args = getopt.getopt(argv[1:], 'hi:a:o:',
    ['help', 'input=', 'annotation=', 'outdir='])
  except getopt.error, msg:
    help_message.usage()

  infile = None
  annotation_file = 'annotation-file'
  outdir = '.'
  for option, value in opts:
    if option == '-h':
      print 'use -i to provide an input file, TCGA copy number file, and -a to provide an annotation file.'
      print 'use -o to provide an output directory'
    if option in ('-i', '--input'):
      infile = value
    if option in ('-o', '--outdir'):
      outdir = value
    if option in ('-a', '--annotation'):
      annotation_file = value
  return infile, annotation_file, outdir

def main(argv=None):
  if argv is None:
    argv = sys.argv

  infile, annotation_file, outdir = get_options(argv)
  type_name = os.path.basename(infile).split('.')[0]
  outfile = os.path.join(outdir, type_name + '.cnv.csv')

  # returns a dataframe indexed by gene name, with chr number and txstart
  annotation_data = process_annotation_file(annotation_file)

  # returns a dict of patient_ids => lists of interval trees containing range data for each chromosome
  patient_data = process_input_file(infile)

  process_and_write_data(outfile, annotation_data, patient_data)

if __name__ == '__main__':
  main()



