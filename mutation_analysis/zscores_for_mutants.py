#!/usr/bin/env python
# encoding: utf-8
'''
zscores_for_mutants.py

Created by Joan Smith
on 2016-2-20.

Given a clinical file and a mutations file, calculate zscores for genes with mutations in >10% of patients

Copyright (c) 2018. All rights reserved.
'''

import pandas as pd
import numpy as np
import argparse
import sys
import os
import glob
import multiprocessing

sys.path.append('../common/')
import utilities as util
import analysis
import mutation_base
import metagene as metagene_lib

MUTATION_PERCENT = .02

def calculate_cox(mutation, clinical, outdir, metagene_file=None, make_km=False):
  clinical_data = util.get_clinical_data(clinical)
  df = mutation_base.prep_mutation_data(mutation, clinical_data)
  clinical_and_data = df.join(clinical_data, how='inner')
  num_patients = len(clinical_and_data)

  #prep output file
  cancer_type = os.path.basename(mutation).split('_')[0].split('.')[0]
  if metagene_file:
    formatstring = '{0}, {1}, {2}, {3}, {4}, {5}\n'
    outfile = os.path.join(outdir, cancer_type + '_mutation-fraction-' + str(MUTATION_PERCENT) + '_metagene_zscores.csv')

    print "Processing metagene..."
    metagene = metagene_lib.get_metagene_data(metagene_file, cancer_type)
    print "Complete"
  else:
    outfile = os.path.join(outdir, cancer_type + '_mutation-fraction-' + str(MUTATION_PERCENT) + '.zscores.out.csv')
    formatstring = '{0}, {1}, {2}, {3}, {4}\n'

  with open(outfile, 'w') as out:
    if metagene_file:
      out.write('gene,zscore,pvalue,metagene-zscore,metagene-pvalue,num patients\n')
    else:
      out.write('gene,zscore,pvalue,num mutations,num patients\n')

    for gene in clinical_and_data:
      if gene in ['time', 'censor']:
        continue

      num_mutations = int(clinical_and_data[gene].sum())
      if num_mutations >= MUTATION_PERCENT * num_patients:
        time = clinical_and_data['time']
        censor = clinical_and_data['censor']
        data = clinical_and_data[gene]

        if metagene_file:
            cox_dict = analysis.do_metagene_cox(time,
                                                censor,
                                                data,
                                                metagene)
            out.write(formatstring.format(gene, cox_dict['z'], cox_dict['p'], cox_dict['metagene-z'],
                                          cox_dict['metagene-p'],
                                          cox_dict['n']))
        else:
          name = cancer_type+ '_' + gene
          if make_km:
            analysis.do_km(name, time, censor, data, outdir)
            clinical_and_data['time', 'censor', gene].to_csv(os.path.join(outdir, name + '_data.csv'),
                                columns=['time', 'censor', 'mutated'])

          cox_dict = analysis.do_cox(time, censor, data)
          out.write(formatstring.format(gene, cox_dict['z'], cox_dict['p'], num_mutations,cox_dict['n']))

def get_options(argv):
  parser = argparse.ArgumentParser(description='Get mutation, cnv, and clinical directories. optional output dir')
  parser.add_argument('-m', action='store', dest='mutation_file')
  parser.add_argument('-c', action='store', dest='clinical')
  parser.add_argument('-g', action='store', dest='metagene_file')
  parser.add_argument('-p', action='store', dest='parallel_workers')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  ns = parser.parse_args()

  return ns.mutation_dir, ns.clinical_dir, ns.outdir, ns.metagene_file


def multiprocess_zscores(args):
  calculate_cox(*args)

def all_cancer_types(mutation_dir, clinical_dir, outdir, metagene=None, parallel_workers=0):
  mutation_files = os.listdir(mutation_dir)
  mutation_files = util.remove_extraneous_files(mutation_files)
  mutation_files = [os.path.join(mutation_dir, f) for f in mutation_files]

  args = []
  for m in mutation_files:
    cancer_type = util.get_cancer_type(m)
    clinical = glob.glob(os.path.join(clinical_dir, '*' + cancer_type + '*'))[0]
    if parallel_workers == 0:
      calculate_cox(m, clinical, outdir, metagene_file=metagene)
    else:
      args.append([m, clinical, outdir, metagene])

  if parallel_workers > 0:
    p = multiprocessing.Pool(parallel_workers)
    p.map(multiprocess_zscores, args)


def main(argv=None):
  mutation_dir, clinical_dir, outdir, metagene_file = get_options(argv)
  all_cancer_types(mutation_dir, clinical_dir, outdir, metagene_file)


if __name__ == "__main__":
  main()
