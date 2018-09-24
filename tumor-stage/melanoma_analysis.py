#!/usr/bin/env python
# encoding: utf-8
'''
Created by Joan Smith
on 2018-09-07.

Multivariate cox analysis to understand melanoma cancer CNA usefulness.

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

cancer_type = 'SKCM'
age_r = 'patient.age_at_initial_pathologic_diagnosis'
breslow_r = 'patient.breslow_depth_value'
gender_r = 'patient.gender'
ulceration_r = 'patient.melanoma_ulceration_indicator'
stage_r = 'patient.stage_event.pathologic_stage'
clark_r = 'patient.melanoma_clark_level_value'


def get_options():
  parser = argparse.ArgumentParser(description='Tumor stage group counts')
  parser.add_argument('-c', action='store', dest='SKCM_clinical')
  parser.add_argument('-i', action='store', dest='clinical_variables')
  parser.add_argument('-m', action='store', dest='SKCM_mut')
  parser.add_argument('-d', action='store', dest='SKCM_cna')
  parser.add_argument('-o', action='store', dest='outdir', default='.')
  ns = parser.parse_args()

  return (ns.SKCM_clinical, ns.clinical_variables, ns.SKCM_cna, ns.SKCM_mut,
          ns.outdir)

def make_clinical_data(clinical_file, clinical_variables, outdir):
  clinical = util.get_clinical_data(clinical_file,
                                    extra_rows=[age_r, breslow_r, gender_r, ulceration_r,
                                    stage_r, clark_r],
                                    extra_rows_numeric=False)

  gender_groups = pd.read_csv(os.path.join(clinical_variables, 'SKCM_gender.csv'),
                                   dtype=str)
  stage_groups = pd.read_csv(os.path.join(clinical_variables, 'SKCM_stage.csv'), dtype=str)
  clark_groups = pd.read_csv(os.path.join(clinical_variables, 'SKCM_clark.csv'), dtype=str)
  ulceration_groups = pd.read_csv(os.path.join(clinical_variables, 'SKCM_ulceration.csv'),
                                  dtype=str)

  clinical = tumor_stage_util.group_discontinuous_vars(clark_r, 'clark', clark_groups, clinical)
  clinical = tumor_stage_util.group_discontinuous_vars(gender_r, 'gender', gender_groups, clinical)

  clinical.to_csv(os.path.join(outdir, cancer_type + '_clinical.csv'),
                               index_label='patient_id')
  clinical[age_r] = pd.to_numeric(clinical[age_r], errors='coerce')
  clinical[breslow_r] = pd.to_numeric(clinical[breslow_r], errors='coerce')
  clinical = clinical.dropna(subset=[breslow_r])
  clinical['breslow_0'] = np.where(clinical[breslow_r] <= 1, 0, 1)
  return clinical

def do_cox_models(clinical, cn_file, mut_file, outdir):
  cn = pd.read_csv(cn_file)
  cn_by_patient = cn.transpose()
  cn_by_patient = cn_by_patient.drop(['Chromosome', 'Location'])
  cn_by_patient.columns = cn_by_patient.loc['Symbol']
  cn = cn_by_patient[['\'NOTCH2', '\'MTAP']]

  mut = mutation_base.prep_mutation_data(mut_file, clinical)
  p53_mut = mut[['\'TP53']]
  p53_mut.columns = ['TP53']

  data = cn.join(clinical, how='inner')
  data = data.join(p53_mut, how='inner')

  analyses = {
    'CNA only': [age_r, 'breslow_0', 'gender_0', 'clark_0'],
    'CNA + P53': ['TP53', age_r, 'breslow_0', 'gender_0', 'clark_0'],
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
  results.T.to_csv(os.path.join(outdir, 'melanoma_analysis.csv'))



def main():
  clinical, clinical_variables, cn_file, mut_file, outdir = get_options()

  clinical = make_clinical_data(clinical, clinical_variables, outdir)
  do_cox_models(clinical, cn_file, mut_file, outdir)


if __name__ == "__main__":
  main()

