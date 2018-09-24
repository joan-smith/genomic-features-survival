#!/usr/bin/env python
# encoding: utf-8
'''
interesting_genes_raw_subtype_data.py

Created by Joan Smith
on 2018-5-6.

Copyright (c) 2018. All rights reserved.
'''

import pandas as pd
import numpy as np
import argparse
import sys
import os
import pdb
import rpy2
import glob

sys.path.append('../common/')
import utilities as util
import analysis

def get_options():
  parser = argparse.ArgumentParser(description='Histologic subtype Zscores')
  parser.add_argument('-c', action='store', dest='clinical_directory')
  parser.add_argument('-i', action='store', dest='histologic_subtype_rows')
  parser.add_argument('-b', action='store', dest='basedir')
  parser.add_argument('-f', action='store', dest='interesting_genes')
  parser.add_argument('-o', action='store', dest='outdir')
  ns = parser.parse_args()

  return (ns.clinical_directory, ns.histologic_subtype_rows,
          ns.basedir, ns.interesting_genes, ns.outdir)

def prep_BRCA_data(extra_data_file, cancer_type):
  extra_data = pd.read_csv(extra_data_file, sep=',', skiprows=[0], index_col=0)
  extra_data['subtype'] = np.nan
  data = extra_data[['ER Status', 'PR Status', 'HER2 Final Status', 'subtype']].copy()

  data.loc[
        (data['HER2 Final Status'] == 'Positive')
      , 'subtype'] = 'HER2+'
  data['subtype'][
        (data['HER2 Final Status'] == 'Negative') &
        ((data['ER Status'] == 'Positive') | (data['PR Status'] == 'Positive'))
      ] = 'ER|PR+ HER2-'
  data['subtype'][
        (data['HER2 Final Status'] == 'Negative') &
        (data['ER Status'] == 'Negative') &
        (data['PR Status'] == 'Negative')
      ] = 'Triple Negative'
  return data

def save_subtype_files(clinical, subtype_col, cancer_type, outdir):
  subtype_counts = clinical[subtype_col].value_counts()
  # print subtype_counts

  subtype_data = []
  for subtype, c in subtype_counts.iteritems():
    subtype_clinical_data = clinical[clinical[subtype_col] == subtype]
    if c >= 100 or (cancer_type == 'SARC' and c >= 50):
      subtype_data.append(subtype_clinical_data)
  return pd.concat(subtype_data)


def make_clinical_data(clinical_file, histologic_subtype_col, outdir):
  cancer_type = util.get_cancer_type(clinical_file)
  clinical = util.get_clinical_data(clinical_file, extra_rows=[histologic_subtype_col],
                  extra_rows_numeric=False)
  return save_subtype_files(clinical, histologic_subtype_col, cancer_type, outdir)



def main():
  clinical_dir, row_names_file, basedir, interesting_genes, outdir = get_options()
  files = os.listdir(clinical_dir)
  files = util.remove_extraneous_files(files)
  clinical_by_cancer_type = {util.get_cancer_type(f): f for f in files}

  row_names = pd.read_csv(row_names_file, header=0)

  interesting_genes = pd.read_csv(interesting_genes, header=0, index_col=1)

  for i, row  in row_names.iterrows():
    cancer_type = row['cancer_type']
    cancer_type_fname = cancer_type
    print cancer_type
    clinical_file = clinical_by_cancer_type[cancer_type]
    clinical_file = os.path.join(clinical_dir, clinical_file)

    if row['histological_subtype_row'] != 'EXTERNAL':
      clinical_data = make_clinical_data(clinical_file, row['histological_subtype_row'], outdir)
    else:
      subtype_data = prep_BRCA_data(row['external_file'], cancer_type)
      # subtype_data.to_csv(os.path.join(outdir, 'BRCA_annotation_subtype_data.csv'))
      clinical = util.get_clinical_data(clinical_file)
      subtype_clinical = clinical.join(subtype_data['subtype'], how='outer')
      clinical_data = save_subtype_files(subtype_clinical, 'subtype', cancer_type, outdir)
      cancer_type_fname = 'BRCA_HER2'


    cna_file = glob.glob(os.path.join(basedir, cancer_type + '*.csv'))[0]
    cna = pd.read_csv(cna_file, header=0, index_col=0).T
    genes = '\'' + interesting_genes['Gene']
    genes = genes.loc[cancer_type]
    print genes
    if type(genes) == str:
      print cna[[genes]]
      joined = cna[[genes]].join(clinical_data, how='outer')
    else:
      joined = cna[genes].join(clinical_data, how='outer')
    joined.to_csv(os.path.join(outdir, cancer_type_fname + '.csv'))


if __name__ == "__main__":
  main()

