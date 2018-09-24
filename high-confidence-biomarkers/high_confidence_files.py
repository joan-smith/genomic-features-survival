#!/usr/bin/env python
# encoding: utf-8
'''
high_confidence_files.py

Created by Joan Smith
on 2017-11-04.

Collect zscores across TCGA, MSKCC, and ICGC for finding high confidence biomarkers
This script produces output files for mutation and cna for each cancer type included in
the input file

Copyright (c) 2018. All rights reserved.
'''

import pandas as pd
import numpy as np
import argparse
import sys
import os
import pdb
import rpy2

sys.path.append('../common/')
import utilities as util
import analysis


SOURCE_FOLDERS = {
  'cBioportal': '~/Dropbox/genomic-features-survival/CBioPortal',
  'MSKCC': '~/Dropbox/genomic-features-survival/MSKCC_met/no-dupes',
  'TCGA': '~/Dropbox/genomic-features-survival/gene-only-time-since-surgery',
  'ICGC': '~/Dropbox/genomic-features-survival/ICGC'
}

CNA_FILES = {
  'cBioportal': 'cnv-0p-cutoff/{0}_CNV.cbioportal_zscores.csv',
  'ICGC': 'copy-number-0p/{0}_zscores.csv',
  'MSKCC': 'cnas-0p-cutoff/{0}',
  'TCGA': 'cna-time-since-surgery/{0}_zscores.csv'
}

MUTATION_FILES = {
  'cBioportal': 'mutation_zscores/{0}/{0}_mutation_percent_0.02.cbioportal_zscores.out.csv',
  'ICGC': 'mutation-zscores/{0}/{0}_mutation_percent_0.02.icgc_zscores.out.csv',
  'MSKCC': 'mutations-2p-cutoff/{0}',
  'TCGA': 'mutation-time-since-surgery-2percent/{0}/{0}.zscores.out.csv'
}


def get_options():
  parser = argparse.ArgumentParser(description='Output univariate cox from pcna25')
  parser.add_argument('-i', action='store', dest='cancer_type_file')
  parser.add_argument('-a', action='store', dest='annotation_file', default='.')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  namespace = parser.parse_args()

  return (namespace.cancer_type_file, namespace.annotation_file, namespace.output_directory)

def get_cna_zscores(source, source_cancer_type):
  source_folder = SOURCE_FOLDERS[source]
  cna_file = os.path.join(source_folder, CNA_FILES[source].format(source_cancer_type))
  cna_file = os.path.expanduser(cna_file)
  if os.path.exists(cna_file):
    df = pd.read_csv(cna_file, index_col=0)
    if source == 'MSKCC':
      df.index = '\'' + df.index
    return df['zscore'].astype(float)
  else:
    print 'missing file for ', source, source_cancer_type

def get_mutation_zscores(source, source_cancer_type):
  source_folder = SOURCE_FOLDERS[source]
  mutations_file = os.path.join(source_folder, MUTATION_FILES[source].format(source_cancer_type))
  mutations_file = os.path.expanduser(mutations_file)
  if os.path.exists(mutations_file):
    df = pd.read_csv(mutations_file, index_col=0)
    if source == 'MSKCC':
      df.index = '\'' + df.index
    return df['zscore']

def get_annotation_data(annotation_file):
  annotation = pd.read_csv(annotation_file, low_memory=False)
  annotation = annotation[['Approved Symbol', 'chromosome', 'Chromosome band', 'txStart']]
  annotation['Gene'] = '\'' + annotation['Approved Symbol']
  annotation = annotation.set_index('Gene')
  annotation = annotation.dropna(how='any')
  annotation = annotation.drop(['Approved Symbol'], axis=1)
  annotation.columns = ['Chromosome', 'Chromosome Band', 'Location']
  return annotation

def main():
  cancer_type_file, annotation_file, outdir = get_options()
  cancer_types = pd.read_csv(cancer_type_file, header=None, index_col=0)
  annotation = get_annotation_data(annotation_file)
  for ctype, row in cancer_types.iterrows():
    print ctype
    row = row.dropna()
    ctype_cnas = pd.DataFrame()
    ctype_mutations = pd.DataFrame()
    for i in row:
      source, source_cancer_type = i.split('_', 1)
      ctype_cna_zscores = get_cna_zscores(source, source_cancer_type)
      if ctype_cna_zscores is not None:
        ctype_cnas[i] = ctype_cna_zscores

      ctype_mutation_zscores = get_mutation_zscores(source, source_cancer_type)
      if ctype_mutation_zscores is not None:
        ctype_mutations[i] = ctype_mutation_zscores


    safe_ctype = ctype.replace(' ', '-')
    ctype_cnas['Stouffers'] = analysis.stouffer_unweighted(ctype_cnas)
    ctype_cnas = annotation.join(ctype_cnas, how='right')
    ctype_cnas.to_csv(os.path.join(outdir, 'cnas', safe_ctype + '.out.csv'))

    ctype_mutations['Stouffers'] = analysis.stouffer_unweighted(ctype_mutations)
    ctype_mutations.to_csv(os.path.join(outdir, 'mutations', safe_ctype + '.out.csv'))


if __name__ == "__main__":
  main()

