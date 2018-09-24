#!/usr/bin/env python
# encoding: utf-8
'''

zscores_for_copy_number_when_mutated.py


Created by Joan Smith
on 2017-7-7.

Given a set of cnv files, and a list of interesting genes with amplifications/deletions, find out how much of the relevant chromosome has a similar aleteration.

Copyright (c) 2018. All rights reserved.

'''

import pandas as pd
import sys
import argparse
import os
import glob

from IPython import embed

sys.path.append('../common/')
import utilities as util
import analysis
import mutation_base

def get_options():
  parser = argparse.ArgumentParser(description='Get cnv directory. csv file, and optional output dir')
  parser.add_argument('-i', action='store', dest='cnv_directory')
  parser.add_argument('-c', action='store', dest='clinical_dir')
  parser.add_argument('-f', action='store', dest='interesting_file')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  namespace = parser.parse_args()

  return (namespace.cnv_directory, namespace.clinical_dir,
          namespace.output_directory, namespace.interesting_file)


def find_continuous_region(patient_data, starting_at=None, alteration_type=None):
  chromosome_length = patient_data.index[-1]

  threshold = patient_data.loc[starting_at]*0.8

  if alteration_type == 'Amplification':
    passes = patient_data > threshold
  else:
    passes = patient_data < threshold


  streaks = passes.mul(passes.cumsum()).diff().where(lambda x: x < 0).ffill().add(
                                                                            passes.cumsum(), fill_value=0)

  # look at the relevant gene's streak value, and count backward that number of locations to get
  # the first element of the streak.
  start_of_streak_pos = streaks.index.get_loc(starting_at) - int(streaks[starting_at]) + 1

  #  get_loc will return a slice if there's a duplicate position. We don't care, so just take the first
  end_of_streak_pos = streaks.index.get_loc(streaks[starting_at:].where(lambda x: x == 0).dropna().index[0])
  if type(end_of_streak_pos) == slice:
    end_of_streak_pos = end_of_streak_pos.start

  relevant_streak = streaks.iloc[start_of_streak_pos:end_of_streak_pos]

  distance_passing_threshold = relevant_streak.index[-1] - relevant_streak.index[0]

  return distance_passing_threshold,  chromosome_length

def copy_number_changes(cnv, clinical,  outdir, cancer_type_genes):
  cancer_type = util.get_cancer_type(cnv)
  print cancer_type
  clinical = util.get_clinical_data(clinical)
  copy_numbers = pd.read_csv(cnv, index_col=0)


  for i, gene in cancer_type_genes.iterrows():
    results = pd.DataFrame()

    gene_name = gene['Gene']
    print gene_name
    gene_cnas = copy_numbers.loc[gene_name]
    chrom = gene_cnas['Chromosome']
    gene_location = copy_numbers.loc[gene_name]['Location']


    if gene['Type'] == 'Amplification':
      threshold_passed = gene_cnas > 0.3
    else:
      threshold_passed = gene_cnas < -0.3
    threshold_passed = threshold_passed.drop(['Chromosome', 'Location'])
    threshold_passed = threshold_passed[threshold_passed]

    copy_numbers_on_same_chrom = copy_numbers[copy_numbers['Chromosome'] == chrom]
    for patient in copy_numbers_on_same_chrom:
      if patient not in clinical.index:
        continue
      if patient in ['Chromosome', 'Location']:
        continue
      if patient in threshold_passed.index:
        patient_data = copy_numbers_on_same_chrom[['Location', patient]]
        patient_data = patient_data.reset_index().sort_values(by='Location') \
                                    .set_index('Location').drop('Symbol')
        continuous, total = find_continuous_region(patient_data[patient],
                                                  starting_at=gene_location,
                                                  alteration_type=gene['Type'])
      else:
        continuous, total = (None, None)
      results[patient] = pd.Series({'continuous_len': continuous,
                          'chr_len': total,
                          'fraction': continuous/total if continuous else None,
                          'copy number': gene_cnas[patient],
                          'time': clinical.loc[patient].time,
                          'censor': clinical.loc[patient].censor})
    results.transpose().to_csv(os.path.join(outdir, cancer_type + '_' + gene_name[1:] + '.cn_changes.csv'),
                            columns=['time', 'censor', 'copy number', 'fraction', 'continuous_len', 'chr_len'])


def multiprocess_copy_number_changes(args):
  cnv = args[0]
  clinical = args[1]
  outdir = args[2]
  cancer_type_genes = args[3]
  copy_number_changes(cnv, clinical, outdir, cancer_type_genes)

def main(argv=None):
  cnv_dir, clinical_dir, outdir, input_file = get_options()
  cnv_files = os.listdir(cnv_dir)
  cnv_files = util.remove_extraneous_files(cnv_files)
  cnv_files = [os.path.join(cnv_dir, i) for i in cnv_files]

  interesting_genes = pd.read_csv(input_file, comment='#')
  interesting_genes['Gene'] = '\'' + interesting_genes['Gene']

  zscore_inputs = []
  for cnv in cnv_files:
    cancer_type = util.get_cancer_type(cnv)
    cancer_type_genes = interesting_genes[interesting_genes['Cancer Type'] == cancer_type]
    if len(cancer_type_genes) == 0:
      continue

    clinical = glob.glob(os.path.join(clinical_dir, cancer_type + '*'))[0]
    multiprocess_copy_number_changes([cnv, clinical, outdir, cancer_type_genes])


if __name__ == "__main__":
  main()
