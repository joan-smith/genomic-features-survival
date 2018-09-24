#!/usr/bin/env python
# encoding: utf-8
'''
process_copy_numbers_to_genes.py

Forked from copy-number/process_copy_numbers_to_genes.py for ICGC (instead of TCGA)

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

SEGMENT_VALUE = {
    'ESAD-UK': 'segment_mean',
    'OV-AU': 'segment_median',
    'PACA-AU': 'segment_median',
    'SKCA-BR': 'segment_mean'
}

def process_input_file(input_file, cancer_type):
  df = pd.read_csv(input_file, sep='\t')
  df['chromosome'] = df['chromosome'].astype(str)

  patient_data = {}
  for index, row in df.iterrows():
    identifier = row['icgc_donor_id']
    chromosome = row['chromosome']
    start = row['chromosome_start']
    end = row['chromosome_end']
    copy_number = row[SEGMENT_VALUE[cancer_type]]

    # Add the new patient to the data dict, initializing the chromosome dict
    if not identifier in patient_data:
      patient_data[identifier] = {}

    # Initialize the interval tree for the new chromosome for this patient
    if chromosome not in patient_data[identifier].keys():
      patient_data[identifier][chromosome] = IntervalTree()

    # Add the range and copy number in this row to the correct patient_data/chromosome location
    # Note the interval tree implementation uses half close intervals, but copy number data
    # uses closed intervals, so we add 1 to the end to ensure our intervaltree matches the data.
    patient_data[identifier][chromosome][start:end+1] = copy_number
  return patient_data

def process_annotation_row(g):
  return pd.Series({'gene': g['Approved Symbol'], 'start': g.txStart,
    'chromosome': g['chromosome']})

def process_annotation_file(annotation_file):
  # columns: 'Approved Symbol', 'Chromosome' 'txStart' (could use 'txEnd' too)
  # plan of record, use txStart as gene location
  annotation = pd.read_csv(annotation_file, low_memory=False)
  processed_annotation = annotation.apply(process_annotation_row, axis=1)
  return processed_annotation.dropna(how="any")

def get_shortest_conaining_interval_val(intervals):
  shortest_interval = None
  shortest_interval_size = sys.maxint
  for i in intervals:
    i_size = i.end - i.begin
    if i_size < shortest_interval_size:
      shortest_interval_size = i_size
      shortest_interval = i
  return shortest_interval.data

def get_copy_number_for_gene_start(chromosome_copy_numbers, start, patient, gene):
  intervals = chromosome_copy_numbers[start]
  if len(intervals) > 1:
    return get_shortest_conaining_interval_val(intervals)
  if len(intervals) == 0:
    # icgc cna files don't include chromosomes that have no copy number alteration,
    # thus if the location is missing, it means that the gene has an unaltered copy number,
    # so we return 0 (since CNAs are measured as log2(copy number)
    return 0
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
      f.write('\'' + str(row.gene) + ',' + str(row.chromosome) + ',' + str(row.start) + ',')
      for i, patient in enumerate(patients):
        all_chromosome_copy_numbers = patient_data[patient]
        # sometimes the whole chromosome has no copy number alterations, which means it will be absent from the dict
        #  in that case, the gene has 0 copy number alteration.
        if row.chromosome in all_chromosome_copy_numbers.keys():
          this_chromosome_copy_numbers = all_chromosome_copy_numbers[row.chromosome]
          copy_number = get_copy_number_for_gene_start(
              this_chromosome_copy_numbers, row.start, patient, row.gene)
        else:
          # print 'Patient ', patient, 'has not chromsome data for ', row.chromosome
          copy_number = 0
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
  type_name = os.path.basename(infile).split('.')[-2]
  outfile = os.path.join(outdir, type_name + '.cnv.csv')

  # returns a dataframe indexed by gene name, with chr number and txstart
  annotation_data = process_annotation_file(annotation_file)

  # returns a dict of patient_ids => lists of interval trees containing range data for each chromosome
  patient_data = process_input_file(infile, type_name)

  process_and_write_data(outfile, annotation_data, patient_data)

if __name__ == '__main__':
  main()



