#!/usr/bin/env python
# encoding: utf-8
'''

zscores_for_copy_number_when_mutated.py


Created by Joan Smith
on 2017-7-7.

Given a set of cnv files, mutation files and clinical files, all from TCGA
calculate multivariate z-scores with cnv for each gene, and whether p53 is
mutated in that patient.

Copyright (c) 2018. All rights reserved.

'''

import pandas as pd
import numpy as np
import argparse
import sys
import os
import pdb
import collections
import glob
import rpy2
from multiprocessing import Pool

sys.path.append('../common/')
import utilities as util
import analysis
import mutation_base

def get_options():
  parser = argparse.ArgumentParser(description='Get mutation, cnv, and clinical directories. optional output dir')
  parser.add_argument('-i', action='store', dest='cnv_directory')
  parser.add_argument('-m', action='store', dest='mutation_directory')
  parser.add_argument('-c', action='store', dest='clinical_directory')
  parser.add_argument('-f', action='store', dest='interesting_file')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  namespace = parser.parse_args()

  return (namespace.cnv_directory, namespace.mutation_directory,
          namespace.clinical_directory, namespace.output_directory,
          namespace.interesting_file)

def calculate_cox(data, gene):
  try:
    cox_dict = analysis.do_multivariate_cox(data.time,
                               data.censor,
                               data[gene],#'\'' + gene.split('_')[0]],
                               data[['TP53_mutation']],
                               float_vars=True)
    return cox_dict
  except rpy2.rinterface.RRuntimeError as e:
    print 'WARN: skipped', gene, 'due to R error'
    return {}


def make_zscores(copy_number, mutation, clinical, outdir, genes):
  cancer_type = util.get_cancer_type(copy_number)

  clinical_data = util.get_clinical_data(clinical)
  cnv = pd.read_csv(copy_number, index_col=0)

  mutation = mutation_base.prep_mutation_data(mutation, clinical_data)
  p53_mutation = mutation['\'TP53'].rename('TP53_mutation')

  cnv_by_patient = cnv.transpose()
  clinical_and_cnv = cnv_by_patient.join(clinical_data, how='inner')

  clinical_mutations_and_cnv = clinical_and_cnv.join(p53_mutation,
                                      how='inner')

  cox_dicts = {}
  for gene in genes['Gene']:
    clinical_gene = clinical_mutations_and_cnv[
        [gene, 'TP53_mutation', 'time', 'censor']]
    cox_dict = calculate_cox(clinical_gene, gene)
    cox_dict['mutation_count'] = clinical_gene['TP53_mutation'].sum()

    clinical_gene.to_csv(os.path.join(outdir, cancer_type + '_' + gene[1:] + '_p53_and_cna_data.csv'))
    cox_dicts[gene[1:]] = cox_dict
  return cox_dicts


def multiprocess_zscores(args):
  copy_number = args[0]
  mutation = args[1]
  clinical = args[2]
  outdir = args[3]
  genes = args[4]
  cox_dicts = make_zscores(copy_number, mutation, clinical, outdir, genes)
  cancer_type = util.get_cancer_type(copy_number)
  return {cancer_type: cox_dicts}


def main(argv=None):
  cnv_dir, mutation_dir, clinical_dir, outdir, input_file = get_options()
  cnv_files = os.listdir(cnv_dir)
  cnv_files = util.remove_extraneous_files(cnv_files)
  cnv_files = [os.path.join(cnv_dir, i) for i in cnv_files]

  interesting_genes = pd.read_csv(input_file, comment='#')
  print interesting_genes
  interesting_genes['Gene'] = '\'' + interesting_genes['Gene']

  zscore_inputs = []
  for cnv in cnv_files:
    cancer_type = util.get_cancer_type(cnv)
    cancer_type_genes = interesting_genes[interesting_genes['Cancer Type'] == cancer_type]
    if len(cancer_type_genes) == 0:
      continue
    print cancer_type
    print cancer_type_genes
    mutation = glob.glob(os.path.join(mutation_dir, cancer_type + '*'))[0]
    clinical = glob.glob(os.path.join(clinical_dir, '*' + cancer_type + '*'))[0]


    zscore_inputs.append([cnv, mutation, clinical, outdir, cancer_type_genes])
    #multiprocess_zscores([cnv, mutation, clinical, outdir, cancer_type_genes])

  p = Pool(4)
  results = p.map(multiprocess_zscores, zscore_inputs)
  with open(os.path.join(outdir, 'cox_results.csv'), 'w') as out:
    formatstr = '{},{},{},{},{},{},{},{}\n'
    out.write('Cancer Type,Gene,CNA Z Score, CNA P value, Mutation Z score, Mutation P Value, Mutation Count, n\n')
    for coxs in results:
      cancer_type = coxs.keys()[0]
      print cancer_type
      for gene, cox_dict in coxs[cancer_type].iteritems():
        print gene, cox_dict
        out.write(formatstr.format(cancer_type, gene,
                    cox_dict['var-z'], cox_dict['var-p'],
                    cox_dict['TP53_mutation-z'], cox_dict['TP53_mutation-p'],
                    cox_dict['mutation_count'], cox_dict['var-n']))


if __name__ == "__main__":
  main()
