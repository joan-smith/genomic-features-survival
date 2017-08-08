#!/usr/bin/env python
# encoding: utf-8
"""
combine-types-to-pancan.py

Given an entrez.csv, gb.csv, orf.csv, ensembl.csv, the survival notes,
and the annotation data, combine into one pancan file
"""

import pandas as pd
import numpy as np
import sys
import argparse
import os


def join_with_genes(df, conversion_path, conversion_file):
  path = os.path.join(conversion_path, conversion_file)
  conversion = pd.read_csv(path)
  conversion = conversion.dropna(how='any')
  conversion = conversion.set_index(u'Annotation')

  joined = df.join(conversion[u'Symbol'], how='outer')
  joined = joined.dropna(subset=[u'Symbol'])

  gene_means = joined.groupby(u'Symbol').mean()
  return gene_means

def join_with_survival_notes(df, notes):
  notes = pd.read_excel(notes)
  notes = notes.dropna(subset=[u'Platform', u'GSE', u'Patient number', u'Cancer type'])
  notes[u'gpl gse'] = notes[u'Platform'].str.strip() + ' ' + notes[u'GSE'].str.strip()
  notes = notes.set_index(u'gpl gse')
  transposed_pancan = df.transpose()
  metadata = notes[[u'Patient number', u'Cancer type']]
  metadata.columns = ['patient count', 'type']
  pancan_with_metadata_tnsp = metadata.join(transposed_pancan)
  with_gse_gpl = make_gpl_gse_rows(pancan_with_metadata_tnsp)
  return with_gse_gpl.transpose()

def make_gpl_gse_rows(pancan_tnsp):
  splitter = lambda x: pd.Series([i for i in x.split(' ')])
  keys = pd.Series(pancan_tnsp.index).apply(splitter)
  keys = keys.set_index(pancan_tnsp.index)
  keys.columns = ['GPL', 'GSE']
  return keys.join(pancan_tnsp)

def do_work(infolder, outfolder, notes, conversions):
  gb_path = os.path.join(infolder, 'GB.csv')
  gb = pd.read_csv(gb_path, index_col=0, low_memory=False)
  gb = gb.astype(float)
  gb_with_gene_names = join_with_genes(gb, conversions, 'genbank to HGNC.csv')

  entrez_path = os.path.join(infolder, 'ENTREZ.csv')
  entrez = pd.read_csv(entrez_path, index_col=0, low_memory=False)
  entrez = entrez.astype(float)
  entrez_with_gene_names = join_with_genes(entrez, conversions, 'Entrez to HGNC.csv')

  orf_path = os.path.join(infolder, 'ORF.csv')
  orf = pd.read_csv(orf_path, index_col=0, low_memory=False)
  orf = orf.astype(float)
  orf = orf.reset_index()
  orf_by_gene = orf.groupby(u'Gene').agg(np.mean)

  ensembl_path = os.path.join(infolder, 'ENSEMBL.csv')
  ensembl = pd.read_csv(ensembl_path, index_col=0, low_memory=False)
  ensembl = ensembl.astype(float)
  ensembl_with_gene_names = join_with_genes(ensembl, conversions, 'ensembl to HGNC.csv')

  df = pd.concat([orf_by_gene,
                 gb_with_gene_names,
                 ensembl_with_gene_names,
                 entrez_with_gene_names], axis=1)

  final = join_with_survival_notes(df, notes)
  out = os.path.join(outfolder, 'pancan_with_metadata.csv')
  final.to_csv(out, header=False)

def get_options():
  parser = argparse.ArgumentParser(description=('Given per-cancer type gene-id files, build',
    'to a pancan file')


  parser.add_argument('-i', action='store', dest='input_dir')
  parser.add_argument('-o', action='store', dest='output_dir')
  parser.add_argument('-c', action='store', dest='conversions')
  parser.add_argument('-n', action='store', dest='notes_file')
  namespace = parser.parse_args()

  return (namespace.input_dir, namespace.output_dir,
          namespace.conversions, namespace.notes_file)

def main(argv=None):
  if argv == None:
    argv = sys.argv
  infolder, outfolder, notes, conversions = get_options(argv)

  do_work(infolder, outfolder, notes, conversions)

if __name__ == '__main__':
  main()
