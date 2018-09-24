#!/usr/bin/env python
# encoding: utf-8
'''
Created by Joan Smith
on 2018-09-07.

Given clinical files, mutation files, and cnv files, get zscores with confidence intervals and hazard
ratios for Figure 1 genes.

Copyright (c) 2018. All rights reserved.
'''

import pandas as pd
import numpy as np

import glob
import argparse
import sys
import os
import pdb
from multiprocessing import Pool

sys.path.append('../common/')
import mutation_base
import utilities as util
import analysis

MUTATION_PERCENT = 0.02

def get_options():
  parser = argparse.ArgumentParser(description='Get cnv, mutation, gene list,  and clinical directories. optional output dir')
  parser.add_argument('-i', action='store', dest='cnv_dir')
  parser.add_argument('-m', action='store', dest='mutation_dir')
  parser.add_argument('-c', action='store', dest='clinical_dir')
  parser.add_argument('-g', action='store', dest='gene_list')
  parser.add_argument('-o', action='store', dest='output_dir', default='.')
  ns = parser.parse_args()

  return (ns.cnv_dir, ns.mutation_dir, ns.clinical_dir,
          ns.gene_list, ns.output_dir)


def make_cnv_zscores(copy_number,clinical, gene_list):
  cancer_type = util.get_cancer_type(copy_number)

  cna = pd.read_csv(copy_number)
  cna_by_patient = cna.transpose()
  cna_by_patient.columns = cna_by_patient.loc['Symbol']
  cna_by_patient_gene_list_only = cna_by_patient[gene_list]

  clinical_data = util.get_clinical_data(clinical)
  clinical_and_cnv = cna_by_patient_gene_list_only.join(clinical_data, how='inner')

  results = pd.DataFrame()
  for gene in clinical_and_cnv:
    if gene in ['time', 'censor']:
      continue
    cox_dict = analysis.do_cox(clinical_and_cnv.time,
                               clinical_and_cnv.censor,
                               clinical_and_cnv[gene])
    cox_dict['cancer_type'] = cancer_type
    cox_dict['gene'] = gene
    results = results.append(cox_dict, ignore_index=True)
  return results


def make_mutation_zscores(mutation, clinical, gene_list):
  cancer_type = util.get_cancer_type(mutation)

  # get mutation patients
  clinical_data = util.get_clinical_data(clinical)
  mutation = mutation_base.prep_mutation_data(mutation, clinical_data)

  present_gene_list = list(set(gene_list.values) & set(mutation.columns.values))
  mutation_gene_list_only = mutation[present_gene_list]

  mutation_and_clinical = mutation_gene_list_only.join(clinical_data, how='inner')
  num_patients = len(mutation_and_clinical.index)

  results = pd.DataFrame()
  for gene in mutation_and_clinical:
    if gene in ['time', 'censor']:
      continue
    num_mutations = mutation_and_clinical[gene].sum()
    if num_mutations >= MUTATION_PERCENT * num_patients:
      cox_dict = analysis.do_cox(mutation_and_clinical.time,
                                 mutation_and_clinical.censor,
                                 mutation_and_clinical[gene])
      cox_dict['cancer_type'] = cancer_type
      cox_dict['gene'] = gene
      cox_dict['num_mutations'] = num_mutations
      results = results.append(cox_dict, ignore_index=True)
  print results
  return results

def main(argv=None):
  cnv_dir, mutation_dir, clinical_dir, gene_list, outdir = get_options()
  gene_list = pd.read_csv(gene_list, header=None)
  gene_list = '\'' + gene_list[0]

  cnv_files = os.listdir(cnv_dir)
  cnv_files = util.remove_extraneous_files(cnv_files)
  cnv_files = [os.path.join(cnv_dir, i) for i in cnv_files]

  mutation_results = pd.DataFrame()
  cnv_results = pd.DataFrame()
  for cnv in cnv_files:
    cancer_type = util.get_cancer_type(cnv)
    print cancer_type
    mutation = glob.glob(os.path.join(mutation_dir, cancer_type + '*'))[0]
    clinical = glob.glob(os.path.join(clinical_dir, '*' + cancer_type + '*'))[0]

    cnv_results = cnv_results.append(make_cnv_zscores(cnv, clinical, gene_list))
    mutation_results = mutation_results.append(make_mutation_zscores(mutation, clinical, gene_list))

  mutation_results.to_csv(os.path.join(outdir, 'mutation_zscores_w_hazards_fig1.csv'), index=False)
  cnv_results.to_csv(os.path.join(outdir, 'cnv_zscores_w_hazards_fig1.csv'), index=False)


if __name__ == "__main__":
  main()
