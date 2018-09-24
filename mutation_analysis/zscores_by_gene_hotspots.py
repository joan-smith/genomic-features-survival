#!/usr/bin/env python
# encoding: utf-8
'''
zscores_by_codon.py

Created by Joan Smith
on 2017-6-04.

Given a file with:
Gene, Chromosome, and a list of positions

For every codon with more than 10 mutations,
Use mutation data and clinical data to calculate zscores for every cancer type, for every row in the input file.
Patients who have a mutation in any of the list of codons associated with a gene are a "1."

Copyright (c) 2017. All rights reserved.
'''

import pdb
import pandas as pd
import numpy as np
import argparse
import sys
import os
import collections

import zscores_for_mutants
sys.path.append('../common/')
import utilities as util
import analysis

MUTATION_PERCENT = .02

def get_options():
  parser = argparse.ArgumentParser(
        description=('Use mutation and clinical data to calculate zscores for patients with mutations'
        'in particular codons'))
  parser.add_argument('-f', action='store', dest='infile', required=True, type=str)
  parser.add_argument('-i', action='store', dest='indir', required=True, type=str)
  parser.add_argument('-o', action='store', dest='outdir', required=False, type=str, default='.')
  namespace = parser.parse_args()

  return namespace.infile, namespace.indir, namespace.outdir

def calculate_cox_for_cancer_type(requested_data, mutation_data, outdir):
  cancer_type = util.get_cancer_type(mutation_data)
  clinical = os.path.join('.', 'clinical', cancer_type + '.clin.merged.txt')
  clinical_data = util.get_clinical_data(clinical)

  start_pos = None
  if cancer_type in ['COADREAD', 'OV']:
    folder = os.path.dirname(mutation_data)
    mutation_data = os.path.join(folder, 'HG36_HG37', cancer_type + '_hg36_hg37.txt')
    start_pos = u'hg37_start'


  df, clinical_with_sequenced_patients, num_patients = zscores_for_mutants.prep_data(mutation_data, clinical_data)
  if not start_pos:
    upper_columns =  [i.upper() for i in df.columns]
    start_pos_index = upper_columns.index('START_POSITION')
    start_pos = df.columns[start_pos_index]

  patients_with_gene = df.groupby(level=u'Hugo_Symbol')
  output_data = []
  for i, request in requested_data.iteritems():
    gene = i[1:]
    # print gene
    # print request
    if gene in patients_with_gene.groups.keys():
      patients_with_requested_gene = patients_with_gene.get_group(gene)
      mutated_at_positions = patients_with_requested_gene[start_pos].isin(request)
       # print mutated_at_positions
      patients_with_requested_positions = patients_with_requested_gene[mutated_at_positions]
      ids_with_requested_positions = patients_with_requested_positions.index.get_level_values('identifier')
      if len(ids_with_requested_positions) >= MUTATION_PERCENT*clinical_with_sequenced_patients.shape[0]:
        analysis_data = pd.DataFrame(
          {'mutated': np.ones(len(ids_with_requested_positions))}, index=ids_with_requested_positions)
        analysis_data = analysis_data.join(clinical_with_sequenced_patients, how='right')
        analysis_data['mutated'].fillna(0, inplace=True)
        cox_dict = analysis.do_cox(analysis_data['time'], analysis_data['censor'], analysis_data['mutated'])

        outdict = {cancer_type + ' p': cox_dict['p']}
        outdict[cancer_type + ' z'] = cox_dict['z']
        outdict[cancer_type + ' mutants'] = len(ids_with_requested_positions)
        outdict[cancer_type + ' n'] = cox_dict['n']
        outdict['gene'] = i
        outdict['positions'] = ':'.join(request)
        output_data.append(outdict)
  outdata = pd.DataFrame(output_data)
  print outdata
  if len(outdata):
    outdata = outdata.set_index(['gene', 'positions'])
  return outdata



def read_requested_data(infile):
  requested_data = collections.defaultdict(list)
  with open(infile) as f:
    f.readline()
    for line in f:
      sp = line.split(',')
      requested_data['gene_and_codon'].append(sp[0])
      requested_data['gene'].append(sp[1])
      requested_data['codon'].append(sp[2])
      requested_data['count'].append(sp[3])
      requested_data['chromosome'].append(sp[4])
      positions = sp[5:]
      positions = [i.strip() for i in positions if len(i.strip()) > 0]
      requested_data['positions'].append(positions)
  return pd.DataFrame(requested_data)

def main(argv=None):
  if argv is None:
    argv = sys.argv
    infile, indir, outdir = get_options()

    requested_data = read_requested_data(infile)
    requested_data = requested_data.groupby('gene')['positions'].apply(
        lambda l: [item for sublist in l for item in sublist])

    files = os.listdir(indir)
    files.remove('.DS_Store')
    files.remove('HG36_HG37')
    output_data = []
    for f in files:
      cancer_type = util.get_cancer_type(f)
      print cancer_type
      zscores = calculate_cox_for_cancer_type(requested_data, os.path.join(indir, f), outdir)
      output_data.append(zscores)
    df = pd.concat(output_data, axis=1)
    df.to_csv('scratch/zscores_by_gene_hotspot.csv')


if __name__ == "__main__":
  main()
