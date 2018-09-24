#!/usr/bin/env python
# encoding: utf-8
'''
Created by Joan Smith
on 2018-9-17.

Copyright (c) 2018. All rights reserved.
'''

import argparse
import urllib
import os
import sys
import pandas as pd

sys.path.append('common')
from common import utilities as util
from mutation_analysis import zscores_for_mutants
from copy_number import process_copy_numbers_to_genes
from copy_number import zscores_for_copy_number

CANCER_TYPES = ['BLCA', 'BRCA', 'COADREAD', 'GBMLGG', 'HNSC', 'KIPAN', 'LIHC', 'LUAD', 'LUSC',
                'OV', 'PAAD', 'PRAD', 'SARC', 'SKCM', 'STES', 'UCEC']

BUCKET_URL = 'http://storage.googleapis.com/public-smith-sheltzer-cancer-analysis/'
COPY_NUMBER_SUFFIX = '.snp__genome_wide_snp_6__broad_mit_edu__Level_3__segmented_scna_minus_germline_cnv_hg19__seg.seg.txt'

def get_options():
  parser = argparse.ArgumentParser(description='Download data and run analysis')
  parser.add_argument('-p', action='store', dest='parallel_workers', type=int, default=0)
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  ns = parser.parse_args()

  return (ns.parallel_workers, ns.output_directory)

def maybe_create_dir(d):
  try:
    os.makedirs(d)
  except OSError:
    pass
  return

def download_data(output_directory):
  maybe_create_dir(os.path.join(output_directory, 'input-data/TCGA-clinical'))
  maybe_create_dir(os.path.join(output_directory, 'input-data/mutation-data'))
  maybe_create_dir(os.path.join(output_directory, 'input-data/copy-number-data'))

  for c in CANCER_TYPES:
    print 'Downloading', c, 'data...'
    clin_file = c + '.clin.merged.txt'
    urllib.urlretrieve(BUCKET_URL + 'TCGA-clinical-files/' + clin_file,
                       os.path.join(output_directory, 'input-data', 'TCGA-clinical', clin_file))

    copy_number_file = c + COPY_NUMBER_SUFFIX
    urllib.urlretrieve(BUCKET_URL + 'copy-number-data/' + copy_number_file,
                      os.path.join(output_directory, 'input-data', 'copy-number-data', copy_number_file))

    mutation_file = c + '.txt'
    urllib.urlretrieve(BUCKET_URL + 'mutation-data/' + mutation_file,
                      os.path.join(output_directory, 'input-data', 'mutation-data', mutation_file))

  urllib.urlretrieve(BUCKET_URL + 'HUGO_Gene_Nomenclature_Committee_Annotations.csv',
              os.path.join(output_directory, 'input-data', 'HUGO_Gene_Nomenclature_Committee_Annotations.csv'))


def single_file_errors(cancer_type, file_path, fkind):
  try:
    clin = pd.read_csv(file_path, nrows=5, sep='\t')
    if len(clin) != 5:
      print cancer_type, fkind, 'file not downloaded correctly at', file_path
      return True
  except IOError as e:
    print cancer_type, fkind, 'file was not found at', file_path
    return True
  return False

def check_data(input_directory):
  print 'Checking Downloads...'
  errors = False
  for c in CANCER_TYPES:
    clin_file_path = os.path.join(input_directory, 'TCGA-clinical', c + '.clin.merged.txt')
    errors |= single_file_errors(c, clin_file_path, 'clinical')

    mutation_file_path = os.path.join(input_directory, 'mutation-data', c + '.txt')
    errors |= single_file_errors(c, mutation_file_path, 'mutation data')

    copy_number_file_path = os.path.join(input_directory, 'copy-number-data', c + COPY_NUMBER_SUFFIX)
    errors |= single_file_errors(c, copy_number_file_path, 'copy number data')

  if errors:
    print 'File download errors, please fix and rerun'
    sys.exit(1)
  print 'Downloads OK...'


def run_univariate_mutation(parallel_workers, output_directory):
  maybe_create_dir(os.path.join(output_directory, 'mutation-zscores'))
  zscores_for_mutants.all_cancer_types(os.path.join(output_directory, 'input-data', 'mutation-data'),
                                    os.path.join(output_directory, 'input-data', 'TCGA-clinical'),
                                    os.path.join(output_directory, 'mutation-zscores'),
                                    parallel_workers=parallel_workers)

def run_univariate_copy_number(parallel_workers, output_directory):
  maybe_create_dir(os.path.join(output_directory, 'copy-number-by-gene'))
  maybe_create_dir(os.path.join(output_directory, 'copy-number-zscores'))

  process_copy_numbers_to_genes.all_cancer_types(
        os.path.join(output_directory, 'input-data', 'copy-number-data'),
        os.path.join(output_directory, 'input-data', 'HUGO_Gene_Nomenclature_Committee_Annotations.csv'),
        os.path.join(output_directory, 'copy-number-by-gene'),
        parallel_workers=parallel_workers)
  zscores_for_copy_number.all_cancer_types(os.path.join(output_directory, 'copy-number-by-gene'),
        os.path.join(output_directory, 'input-data', 'TCGA-clinical'),
        os.path.join(output_directory, 'copy-number-zscores'),
        parallel_workers=parallel_workers)


def main():
  parallel_workers, output_directory = get_options()
  download_data(output_directory)
  check_data(os.path.join(output_directory, 'input-data'))

  run_univariate_mutation(parallel_workers, output_directory)
  run_univariate_copy_number(parallel_workers, output_directory)


if __name__ == "__main__":
  main()
