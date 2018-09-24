#!/usr/bin/env python
# encoding: utf-8
'''
Created by Joan Smith
on 2018-09-08.

Multivariate cox analysis to understand breast cancer CNA usefulness from CBioportal data.

Copyright (c) 2018. All rights reserved.
'''

import pandas as pd
import numpy as np
import argparse
import sys
import os
import pdb
import rpy2
import pprint

import cbioportal_util

sys.path.append('../common/')
import utilities as util
import analysis
import tumor_stage_util

cancer_type = 'BRCA'
age_r = 'AGE_AT_DIAGNOSIS'
er_r = 'ER_STATUS'
pr_r = 'PR_STATUS'
her2_r = 'HER2_STATUS'
stage_r = 'TUMOR_STAGE'


def get_options():
  parser = argparse.ArgumentParser(description='Tumor stage group counts')
  parser.add_argument('-c', action='store', dest='BRCA_clinical')
  parser.add_argument('-i', action='store', dest='clinical_variables')
  parser.add_argument('-d', action='store', dest='BRCA_cna')
  parser.add_argument('-o', action='store', dest='outdir', default='.')
  ns = parser.parse_args()

  return (ns.BRCA_clinical, ns.clinical_variables, ns.BRCA_cna,
          ns.outdir)

def make_clinical_data(clinical_file, clinical_variables, outdir):
  clinical = pd.read_csv(clinical_file, index_col=0)
  clinical = clinical[[age_r, er_r, pr_r, her2_r, stage_r, 'Time', 'Censor']]

  stage_groups = pd.read_csv(os.path.join(clinical_variables, 'BRCA_stage.csv'), dtype=str)
  er_groups = pd.read_csv(os.path.join(clinical_variables, 'BRCA_er.csv'), dtype=str)
  pr_groups = pd.read_csv(os.path.join(clinical_variables, 'BRCA_pr.csv'), dtype=str)
  her2_groups = pd.read_csv(os.path.join(clinical_variables, 'BRCA_her2.csv'), dtype=str)

  clinical = tumor_stage_util.group_discontinuous_vars(stage_r, 'stage', stage_groups, clinical)
  clinical = tumor_stage_util.group_discontinuous_vars(er_r, 'er', er_groups, clinical)
  clinical = tumor_stage_util.group_discontinuous_vars(pr_r, 'pr', pr_groups, clinical)
  clinical = tumor_stage_util.group_discontinuous_vars(her2_r, 'her2', her2_groups, clinical)

  clinical['combined_er_pr'] = np.where(clinical['er_0'] & clinical['pr_0'], 1, 0)
  clinical[age_r] = pd.to_numeric(clinical[age_r], errors='coerce')

  clinical.to_csv(os.path.join(outdir, cancer_type + '_clinical.csv'),
                               index_label='patient_id')
  return clinical

def do_cox_models(clinical, cn_file, outdir):
  cn = pd.read_csv(cn_file, sep='\t', index_col=0)
  cn_by_patient = cn.transpose()
  cn_by_patient = cn_by_patient.drop(['Entrez_Gene_Id'])
  cn = cn_by_patient[['MYC']]

  data = cn.join(clinical, how='inner')

  analyses = {
    'CNA only': [age_r, 'her2_0', 'combined_er_pr', 'stage_0', 'stage_1'],
  }
  results = pd.DataFrame()
  pp = pprint.PrettyPrinter(indent=2)
  for g in cn:
    for name, a in analyses.iteritems():
      cox_dict = analysis.do_multivariate_cox(data.Time,
                               data.Censor,
                               data[g],
                               data[a],
                               float_vars = True)
      cox_dict['gene'] = name + ' ' + g
      results = results.append(cox_dict, ignore_index=True)

  results = results.set_index('gene')
  results.T.to_csv(os.path.join(outdir, 'breast_analysis.csv'))


def main():
  brca_clinical, clinical_variables, cn_file, outdir = get_options()

  clinical = make_clinical_data(brca_clinical, clinical_variables, outdir)
  do_cox_models(clinical, cn_file, outdir)


if __name__ == "__main__":
  main()

