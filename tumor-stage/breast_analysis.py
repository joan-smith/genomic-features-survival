#!/usr/bin/env python
# encoding: utf-8
'''
Created by Joan Smith
on 2018-09-08.

Multivariate cox analysis to understand breast cancer CNA usefulness.

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

sys.path.append('../common/')
import utilities as util
import analysis
import tumor_stage_util
import mutation_base

cancer_type = 'BRCA'
age_r = 'patient.age_at_initial_pathologic_diagnosis'
er_r = 'patient.breast_carcinoma_estrogen_receptor_status'
pr_r = 'patient.breast_carcinoma_progesterone_receptor_status'
her2_r = 'patient.lab_proc_her2_neu_immunohistochemistry_receptor_status'
stage_r = 'patient.stage_event.pathologic_stage'


def get_options():
  parser = argparse.ArgumentParser(description='Tumor stage group counts')
  parser.add_argument('-c', action='store', dest='BRCA_clinical')
  parser.add_argument('-i', action='store', dest='clinical_variables')
  parser.add_argument('-m', action='store', dest='BRCA_mut')
  parser.add_argument('-d', action='store', dest='BRCA_cna')
  parser.add_argument('-o', action='store', dest='outdir', default='.')
  ns = parser.parse_args()

  return (ns.BRCA_clinical, ns.clinical_variables, ns.BRCA_cna, ns.BRCA_mut,
          ns.outdir)

def make_clinical_data(clinical_file, clinical_variables, outdir):
  clinical = util.get_clinical_data(clinical_file, extra_rows=[age_r, er_r, pr_r, her2_r, stage_r],
                                    extra_rows_numeric=False)


  stage_groups = pd.read_csv(os.path.join(clinical_variables, 'BRCA_stage.csv'), dtype=str)
  er_groups = pd.read_csv(os.path.join(clinical_variables, 'BRCA_er.csv'), dtype=str)
  pr_groups = pd.read_csv(os.path.join(clinical_variables, 'BRCA_pr.csv'), dtype=str)
  her2_groups = pd.read_csv(os.path.join(clinical_variables, 'BRCA_her2.csv'), dtype=str)

  clinical = tumor_stage_util.group_discontinuous_vars(stage_r, 'stage', stage_groups, clinical)
  clinical = tumor_stage_util.group_discontinuous_vars(er_r, 'er', er_groups, clinical)
  clinical = tumor_stage_util.group_discontinuous_vars(pr_r, 'pr', pr_groups, clinical)
  clinical = tumor_stage_util.group_discontinuous_vars(her2_r, 'her2', her2_groups, clinical)

  clinical['combined_er_pr'] = np.where(clinical['er_0'] & clinical['pr_0'], 1, 0)

  clinical.to_csv(os.path.join(outdir, cancer_type + '_clinical.csv'),
                               index_label='patient_id')
  clinical[age_r] = pd.to_numeric(clinical[age_r], errors='coerce')
  return clinical

def do_cox_models(clinical, cn_file, mut_file, outdir):
  cn = pd.read_csv(cn_file)
  cn_by_patient = cn.transpose()
  cn_by_patient = cn_by_patient.drop(['Chromosome', 'Location'])
  cn_by_patient.columns = cn_by_patient.loc['Symbol']
  cn = cn_by_patient[['\'MYC']]

  mut = mutation_base.prep_mutation_data(mut_file, clinical)
  p53_mut = mut[['\'TP53']]
  p53_mut.columns = ['TP53']

  data = cn.join(clinical, how='inner')
  data = data.join(p53_mut, how='inner')

  analyses = {
    'CNA only': [age_r, 'her2_0', 'combined_er_pr', 'stage_0', 'stage_1'],
    'CNA + P53': ['TP53', age_r, 'her2_0', 'combined_er_pr', 'stage_0', 'stage_1']
  }
  results = pd.DataFrame()
  pp = pprint.PrettyPrinter(indent=2)
  for g in cn:
    for name, a in analyses.iteritems():
      cox_dict = analysis.do_multivariate_cox(data.time,
                               data.censor,
                               data[g],
                               data[a],
                               float_vars = True)
      cox_dict['gene'] = name + ' ' + g
      results = results.append(cox_dict, ignore_index=True)

  cox_dict = analysis.do_multivariate_cox(data.time,
                           data.censor,
                           data['TP53'],
                           data[analyses['CNA only']],
                           float_vars = True)
  cox_dict['gene'] = 'TP53 mut'
  results = results.append(cox_dict, ignore_index=True)

  results = results.set_index('gene')
  results.T.to_csv(os.path.join(outdir, 'breast_analysis.csv'))


def main():
  brca_clinical, clinical_variables, cn_file, mut_file, outdir = get_options()

  clinical = make_clinical_data(brca_clinical, clinical_variables, outdir)
  do_cox_models(clinical, cn_file, mut_file, outdir)


if __name__ == "__main__":
  main()

