#!/usr/bin/env python
# encoding: utf-8
'''
Created by Joan Smith
on 2018-09-07.

Multivariate cox analysis to understand prostate cancer CNA usefulness.

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
import mutation_base
import tumor_stage_util

cancer_type = 'PRAD'
gleason_r = 'patient.clinical_cqcf.gleason_score_combined'
stage_r = 'patient.stage_event.tnm_categories.pathologic_categories.pathologic_t'
age_r = 'patient.age_at_initial_pathologic_diagnosis'
psa_r = 'patient.clinical_cqcf.psa_result_preop'
lymph_r = 'patient.number_of_lymphnodes_positive_by_he'

def get_options():
  parser = argparse.ArgumentParser(description='Tumor stage group counts')
  parser.add_argument('-c', action='store', dest='PRAD_clinical')
  parser.add_argument('-i', action='store', dest='clinical_variables')
  parser.add_argument('-m', action='store', dest='PRAD_mut')
  parser.add_argument('-d', action='store', dest='PRAD_cna')
  parser.add_argument('-o', action='store', dest='outdir', default='.')
  ns = parser.parse_args()

  return (ns.PRAD_clinical, ns.clinical_variables, ns.PRAD_cna, ns.PRAD_mut,
          ns.outdir)

def make_clinical_data(clinical_file, clinical_variables, outdir):
  clinical = util.get_clinical_data(clinical_file, extra_rows=[gleason_r, stage_r,
                                                              age_r, psa_r,
                                                              lymph_r],
                                    extra_rows_numeric=False)

  gleason_groups = pd.read_csv(os.path.join(clinical_variables, 'PRAD_gleason.csv'), dtype=str)
  stage_groups = pd.read_csv(os.path.join(clinical_variables, 'PRAD_stage.csv'), dtype=str)

  clinical = tumor_stage_util.group_discontinuous_vars(gleason_r, 'gleason', gleason_groups, clinical)

  clinical[psa_r] = pd.to_numeric(clinical[psa_r], errors='coerce')
  clinical[age_r] = pd.to_numeric(clinical[age_r], errors='coerce')
  clinical = clinical.dropna(subset=[psa_r])
  clinical['psa_0'] = np.where(clinical[psa_r] < 7, 0, 1)
  clinical.to_csv(os.path.join(outdir, cancer_type + '_clinical.csv'),
                               index_label='patient_id')

  return clinical

def do_cox_models(clinical, cn_file, mut_file, outdir):
  cn = pd.read_csv(cn_file)
  cn_by_patient = cn.transpose()
  cn_by_patient = cn_by_patient.drop(['Chromosome', 'Location'])
  cn_by_patient.columns = cn_by_patient.loc['Symbol']
  cn = cn_by_patient[['\'MDM4', '\'CCND1', '\'KDM5B', '\'FGF4']]

  mut = mutation_base.prep_mutation_data(mut_file, clinical)
  p53_mut = mut[['\'TP53']]
  p53_mut.columns = ['TP53']

  data = cn.join(clinical, how='inner')
  data = data.join(p53_mut, how='inner')

  analyses = {
    'CNA only': [age_r, 'psa_0', 'gleason_0'],
    'CNA + P53': ['TP53', age_r, 'psa_0', 'gleason_0'],
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
  results.T.to_csv(os.path.join(outdir, 'prostate_analysis.csv'))



def main():
  prad_clinical, clinical_variables, cn_file, mut_file, outdir = get_options()

  clinical = make_clinical_data(prad_clinical, clinical_variables, outdir)
  do_cox_models(clinical, cn_file, mut_file, outdir)


if __name__ == "__main__":
  main()

