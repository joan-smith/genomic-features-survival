#!/usr/bin/env python
# encoding: utf-8
"""
probes_to_gene_ids_pandas.py

Created by Joan Smith
on 2017-8-6.

Copyright (c) 2017 . All rights reserved.
"""

import sys
import argparse
import os
import glob
import re
import numpy as np
import pandas as pd
from multiprocessing import Pool


from probe_to_gene_ids_helpers import *
COLUMN_FUNCTION_FOR_DB = {
    'GB': get_probe_and_gb_col,
    'ORF': get_probe_and_orf_col,
    'ENTREZ': lambda gpl, sheet: get_probe_and_col_strmatch(sheet, 'entrez'),
    'ENSEMBL': lambda gpl, sheet: get_probe_and_col_strmatch(sheet, 'ensembl'),
}


GENE_DB_GPLS = {
  'GB':      ['GPL1390', 'GPL1479', 'GPL3877', 'GPL5049', 'GPL5186', 'GPL6650', 'GPL7015', 'GPL7759', 'GPL13425'],
  'ENSEMBL': ['GPL5345'],
  'ENTREZ':  ['GPL80', 'GPL96', 'GPL570', 'GPL571', 'GPL887', 'GPL1708', 'GPL2986', 'GPL3883',
              'GPL3921', 'GPL4133', 'GPL6102', 'GPL6480', 'GPL6884', 'GPL6947', 'GPL7264', 'GPL8300',
              'GPL8542', 'GPL10558', 'GPL15048', 'GPL6883', 'GPL4372'],
  'ORF':     ['GPL20769', 'GPL9053', 'GPL10295', 'GPL6098', 'GPL15973']
}


GENE_DB_SPLIT_PATTERNS = {
  'GB': ',',
  'ORF': ',|;',
  'ENTREZ': ',| | /// |;|///|\|',
  'ENSEMBL': ',| |;',
}

def UNUSED_process_nonspecific_probes():
  '''used to save the code snipped to take a split gene_id and turn it into multiple rows with the same probe'''
  split_annotation_sheet = annotation_sheet.assign(
    **{gene_id_col:annotation_sheet[gene_id_col].str.split(pattern)}
  )

  annotation_sheet = pd.DataFrame({
      probe_col:np.repeat(
        split_annotation_sheet[probe_col].values, split_annotation_sheet[gene_id_col].str.len()
      )
    }).assign(
        **{gene_id_col:np.concatenate(split_annotation_sheet[gene_id_col].values)
          })[split_annotation_sheet.columns.tolist()]
  annotation_sheet[gene_id_col] = annotation_sheet[gene_id_col].str.strip()

def get_annotation_sheet(datasets, gpl):
  annotation_path = os.path.join(datasets, gpl, gpl + '.xlsx')
  annotation = pd.read_excel(annotation_path, sheet='Sheet1', na_values='null')
  return annotation

def clean_annotation_sheet(db, probe_col, gene_id_col, annotation_sheet):
  ''' Removes nans, splits gene ids, returns new dataframe with one row per probe, gene id pair.'''
  annotation_sheet = annotation_sheet.dropna(subset=[gene_id_col])
  annotation_sheet = annotation_sheet[[probe_col, gene_id_col]]
  print 'TYPE', annotation_sheet[gene_id_col].dtype
  if db in ['ENTREZ', 'ENSEMBL'] and annotation_sheet[gene_id_col].dtype == float:
    # entrez has some files with integer-like gene ids. Let's not munge them to have ".0" at the end
    annotation_sheet[gene_id_col] = annotation_sheet[gene_id_col].map(lambda x: '{:.0f}'.format(x))
  annotation_sheet[gene_id_col] = annotation_sheet[gene_id_col].astype(str)
  annotation_sheet[probe_col] = annotation_sheet[probe_col].astype(str)

  # filter out all the nonspecific probes
  pattern = GENE_DB_SPLIT_PATTERNS[db]
  annotation_sheet['split'] = annotation_sheet[gene_id_col].str.split(pattern)
  # small mess in order to rule out empty strings as cause for the nonspecificity
  annotation_sheet['nonspecific'] = annotation_sheet['split'].apply(lambda x: len([i for i in x if len(i)])) > 1
  annotation_sheet = annotation_sheet[~annotation_sheet['nonspecific']]

  if db == 'ORF' and not '\'' in annotation_sheet[gene_id_col][0]:
    annotation_sheet[gene_id_col] = '\'' + annotation_sheet[gene_id_col]

  return annotation_sheet


def get_gses_in_gpl(datasets, gpl):
  gse_path = os.path.join(datasets, gpl, 'GSE*.csv')
  gses = glob.glob(gse_path)

  # the first group is 1, not 0 since 0 holds the whole match
  return [re.match('.*(GSE\d+).*', gse).group(1) for gse in gses]


def get_zscore_files_from_gses(infolder, gses, gpl):
  survival_files = []
  for gse in gses:
    gse_path = os.path.join(infolder, gpl, gse + '*.csv')
    gse_matches = glob.glob(gse_path)
    if len(gse_matches) == 1:
      survival_files.append(gse_matches[0])
    elif len(gse_matches) > 1:
      no_match = True
      for gse_match in gse_matches:
        if gpl in gse_match:
          survival_files.append(gse_match)
          no_match = False
      if no_match:
        print 'Warn: no matching univariate cox outputs found for', gse, gpl
    else:
      print 'Warn: no univariate cox outputs found for', gse, gpl
  return survival_files


def get_genes_and_zscores(zscore_file, annotations, probe_col):
  zscores = pd.read_csv(zscore_file, header=2)
  zscores['Gene/Probe'] = zscores['Gene/Probe'].astype(str)
  zscores.set_index('Gene/Probe', inplace=True)
  return annotations.join(zscores, on=probe_col, how='inner')

def get_gse_from_filename(filename):
  return re.match('.*(GSE\d+).*', filename).group(1)

def make_output_column(genes_and_zscores, gene_order_dict , gene_set_size):
  output = [np.nan] * (gene_set_size)
  for gene, zscore in genes_and_zscores:
    ind = gene_order_dict[gene]
    while not np.isnan(output[ind]):
      ind += 1
    output[ind] = zscore
  return output


def get_options():
  parser = argparse.ArgumentParser(description=('Use per-gpl annotation file to convert'
  'from probe to gene database ids (e.g. gb, ensembl)'))

  parser.add_argument('-i', action='store', dest='input_dir')
  parser.add_argument('-o', action='store', dest='output_dir')
  parser.add_argument('-d', action='store', dest='raw_datasets_dir')
  parser.add_argument('-n', action='store', dest='notes_file')
  namespace = parser.parse_args()

  return (namespace.input_dir, namespace.output_dir,
          namespace.raw_datasets_dir, namespace.notes_file)

def process_zscore_file(args):
  zscore_file, gpl, annotation_sheet, probe_col, gene_id_col = args
  gse = get_gse_from_filename(zscore_file)

  zscores_with_gene_ids = get_genes_and_zscores(zscore_file, annotation_sheet, probe_col)
  gene_duplicate_count = zscores_with_gene_ids.groupby(gene_id_col)[gene_id_col].count()
  zscores_with_gene_ids['gene_id_count'] =  zscores_with_gene_ids.groupby(gene_id_col)[gene_id_col].transform(
      lambda x: pd.Series(range(x.count())))
  zscores_with_gene_ids.set_index([gene_id_col, 'gene_id_count'], inplace=True)
  gse_zscore_series = zscores_with_gene_ids[' gene Z Score']

  gse_zscore_series.rename(gpl + ' ' + gse, inplace=True)
  return gse_zscore_series


def main(argv=None):
  input_dir, output_dir, raw_datasets_dir, notes_file = get_options()
  p = Pool(4)

  for db in GENE_DB_GPLS.keys():
    print db
    zscore_series = []
    for gpl in GENE_DB_GPLS[db]:
      print gpl
      annotation_sheet = get_annotation_sheet(raw_datasets_dir, gpl)
      probe_col, gene_id_col = COLUMN_FUNCTION_FOR_DB[db](gpl, annotation_sheet)
      annotation_sheet = clean_annotation_sheet(db, probe_col, gene_id_col, annotation_sheet)
      gses = get_gses_in_gpl(raw_datasets_dir, gpl)
      zscore_files = get_zscore_files_from_gses(input_dir, gses, gpl)
      print 'GSE count for GPL:', len(zscore_files)
      args = [(f, gpl, annotation_sheet, probe_col, gene_id_col) for f in zscore_files]
      gpl_zscores = p.map(process_zscore_file, args)
      zscore_series.extend(gpl_zscores)
    print 'NUMBER ZSCORE FILES', len(zscore_series)
    df = pd.concat(zscore_series, axis=1)
    print df
    df = df.reset_index()
    df = df.drop(df.columns[1], axis=1)
    df.set_index(df.columns[0], inplace=True)
    df.to_csv(db + '.csv', index_label='Gene')

if __name__ == '__main__':
  main()
