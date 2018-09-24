#!/usr/bin/env python
# encoding: utf-8
'''
find_mutated_codons.py

Created by Joan Smith
on 2017-5-29.

Given a set of mutation files, find count the codon mutations per cancer type, and in total.
Used for finding mutation hotspots.

Copyright (c) 2018. All rights reserved.
'''

import pdb

import pandas as pd
import numpy as np
import argparse
import sys
import os

sys.path.append('../common/')
import utilities as util

def get_options():
  parser = argparse.ArgumentParser(description='Count codon mutations across cancer types')
  parser.add_argument('-d', action='store', dest='data', default='.', type=str)
  parser.add_argument('-o', action='store', dest='outdir', default='.', type=str)
  namespace = parser.parse_args()

  return namespace.data, namespace.outdir

def count_codons_in_file(f):
  cancer_type = util.get_cancer_type(f)
  print cancer_type

  df = pd.read_csv(f, sep='\t', low_memory=False)
  # Some of the columns are named Start_position. Others are Start_Position. some are start_position. :|
  upper_columns =  [i.upper() for i in df.columns]
  start_pos_index = upper_columns.index('START_POSITION')
  start_pos = df.columns[start_pos_index]
  chromosome = u'Chromosome'

  ncbi_builds = df[u'NCBI_Build'].value_counts()
  if '36' in ncbi_builds.index:
    print 'Using translated NCBI build', cancer_type
    folder = os.path.dirname(f)
    new_path = os.path.join(folder, 'HG36_HG37', cancer_type + '_hg36_hg37.txt')
    print new_path
    df = pd.read_csv(new_path, sep='\t', dtype=str)
    start_pos = u'hg37_start'
    chromsome = u'hg37_chr'
  wild_type_allele_col = u'Reference_Allele'

  df[u'Tumor_Sample_Barcode'] = df[u'Tumor_Sample_Barcode'].str.strip()
  df = df[~df[u'Hugo_Symbol'].str.contains('Hugo_Symbol')]
  df = df[df[u'Variant_Classification'].str.contains('Missense')] # only include missense
  df[u'Hugo_Symbol'] = '\'' + df[u'Hugo_Symbol'].astype(str)
  df = util.add_identifier_column(df, u'Tumor_Sample_Barcode')

  # Some files have the same mutation from different samples listed under one patient,
  # we only care about the number of patients with a given mutation, so drop duplicates
  df = df.drop_duplicates(
        subset=[u'Hugo_Symbol', chromosome, start_pos, u'identifier'],
        keep='last')
  counts = df.groupby([u'Hugo_Symbol', chromosome, start_pos, wild_type_allele_col]).size()
  count_df = pd.DataFrame(counts)
  count_df.columns = [cancer_type]
  count_df.index.rename(['Gene', 'Chromosome', 'Start Position', 'Wild Type Allele'], inplace=True)
  return count_df


def count_codons(data, outdir):
  files = os.listdir(data)
  files = util.remove_extraneous_files(files)
  files.remove('HG36_HG37')

  outdata = []
  ncbi_outdata = []
  for f in files:
    file_name = os.path.join(data, f)
    cancer_type = util.get_cancer_type(file_name)
    codon_counts = count_codons_in_file(file_name)
    outdata.append(codon_counts)
  df = pd.concat(outdata, axis=1, verify_integrity=True)
  df['sum'] = df.sum(axis=1)
  df.to_csv('codon_counts.csv', index_label=['Gene','Chromosome', 'Start Position', 'Wild Type Allele'])

def main(argv=None):
  argv = sys.argv
  data, outdir = get_options()
  count_codons(data, outdir)


if __name__ == "__main__":
  main()

