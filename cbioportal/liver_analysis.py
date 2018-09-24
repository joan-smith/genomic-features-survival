#!/usr/bin/env python
# encoding: utf-8
'''
Created by Joan Smith
on 2018-09-08.

Multivariate cox analysis to understand liver cancer CNA usefulness from CBioportal data.

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

cancer_type = 'LIHC'
age_r = 'AGE'
cirrhosis_r = 'CIRRHOSIS'
grade_r = 'GRADE'
hep_r = 'HISTOLOGICAL_SUBTYPE'


def get_options():
  parser = argparse.ArgumentParser(description='Tumor stage group counts')
  parser.add_argument('-c', action='store', dest='LIHC_clinical')
  parser.add_argument('-i', action='store', dest='clinical_variables')
  parser.add_argument('-d', action='store', dest='LIHC_cna')
  parser.add_argument('-o', action='store', dest='outdir', default='.')
  ns = parser.parse_args()

  return (ns.LIHC_clinical, ns.clinical_variables, ns.LIHC_cna,
          ns.outdir)

def make_clinical_data(clinical_file, clinical_variables, outdir):
  clinical = pd.read_csv(clinical_file, index_col=0, dtype=str)
  clinical = clinical[[age_r, cirrhosis_r, grade_r, hep_r, 'Time', 'Censor']]

  cirrhosis_groups = pd.read_csv(os.path.join(clinical_variables, 'LIHC_cirrhosis.csv'), dtype=str)
  clinical = tumor_stage_util.group_discontinuous_vars(cirrhosis_r, 'cirrhosis',
                                                       cirrhosis_groups, clinical)

  clinical[age_r] = pd.to_numeric(clinical[age_r], errors='coerce')
  clinical['Time'] = pd.to_numeric(clinical['Time'], errors='coerce')
  clinical['Censor'] = pd.to_numeric(clinical['Censor'], errors='coerce')
  clinical['HBV'] = np.where(clinical[hep_r] == 'HBV', 1, 0)
  clinical['HCV'] = np.where(clinical[hep_r] == 'HCV', 1, 0)
  clinical.to_csv(os.path.join(outdir, cancer_type + '_clinical.csv'),
                               index_label='patient_id')
  return clinical

def do_cox_models(clinical, cn_file, outdir):
  cn = pd.read_csv(cn_file, sep='\t', index_col=0)
  cn_by_patient = cn.transpose()
  cn_by_patient = cn_by_patient.drop(['Entrez_Gene_Id'])
  cn = cn_by_patient[['KLF6', 'FBLN1']]

  data = cn.join(clinical, how='inner')

  analyses = {
    'CNA only': [age_r, 'cirrhosis_0', 'HBV', 'HCV'],
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
  results.T.to_csv(os.path.join(outdir, 'liver_analysis.csv'))


def main():
  lihc_clinical, clinical_variables, cn_file, outdir = get_options()

  clinical = make_clinical_data(lihc_clinical, clinical_variables, outdir)
  do_cox_models(clinical, cn_file, outdir)


if __name__ == "__main__":
  main()

