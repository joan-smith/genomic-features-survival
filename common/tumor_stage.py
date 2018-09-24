"""
Helper methods for producing tumor stage output files
"""
import pandas as pd
import numpy as np
import os
import glob
import sys

import utilities as util

TUMOR_STAGE = {
    'BLCA': 'patient.stage_event.pathologic_stage',
    'BRCA': 'patient.stage_event.pathologic_stage',
    'COADREAD': 'patient.stage_event.pathologic_stage',
    'GBMLGG': 'admin.disease_code',
    'HNSC': 'patient.stage_event.pathologic_stage',
    'KIPAN': 'patient.stage_event.pathologic_stage',
    'LIHC': 'patient.stage_event.pathologic_stage',
    'LUAD': 'patient.stage_event.pathologic_stage',
    'LUSC': 'patient.stage_event.pathologic_stage',
    'OV': 'patient.stage_event.clinical_stage',
    'PAAD': 'patient.stage_event.pathologic_stage',
    'PRAD': 'patient.stage_event.tnm_categories.pathologic_categories.pathologic_t',
    'SARC': None,
    'SKCM': 'patient.stage_event.pathologic_stage',
    'STES': 'patient.stage_event.pathologic_stage',
    'UCEC': 'patient.stage_event.clinical_stage',
}

TUMOR_GRADE = 'patient.neoplasm_histologic_grade'
TUMOR_GRADE_TYPES = ['BLCA', 'GBMLGG', 'HNSC', 'KIPAN', 'LIHC', 'OV', 'PAAD', 'STES', 'UCEC']


# slightly messy method to take a data frame of groups and
# dicotimize the relevant clinical data
# e.g:

# row 1: stage a
# row 2: stage b, stage c
#
# produces clinical data with one new column that has called "name_0" that
# has 1s if the row contains stage b or c, and 0 otherwise.
def group_discontinuous_vars(row, name, groups, clinical):
  groups = groups.dropna(how='all')
  included_values = []
  for i, group in groups.iterrows():
    g = group.dropna().values
    if len(g) > 0:
      # print ', '.join(g) +  ': ', clinical[clinical[row].isin(g)][row].count()
      included_values.extend(g)

      clinical[name + '_' + str(i)] = np.where(
                          clinical[row].isin(included_values),
                          0, 1)

  clinical = clinical.drop(name + '_' + str(i), axis=1) # remove last group, since it's empty
  dropped_clinical = clinical[~clinical[row].isin(included_values)] # note excluded values
  if len(dropped_clinical) > 0:
    print 'WARN: Excluded some clinical data for',  name, dropped_clinical[row].value_counts().index.values

  clinical = clinical[clinical[row].isin(included_values)]
  return clinical

def prep_tumor_stage_data(tumor_stage_data_dir, cancer_type):
  tumor_stage_path = os.path.join(tumor_stage_data_dir, cancer_type + '_clinical.csv')
  if not os.path.isfile(tumor_stage_path):
    return None, None
  tumor_stage_data = pd.read_csv(tumor_stage_path, sep=',')
  tumor_stage_data = util.add_identifier_column(tumor_stage_data, 'patient_id')
  tumor_stage_data = tumor_stage_data.drop('patient_id', axis=1)
  tumor_stage_data = tumor_stage_data.set_index('identifier')
  tumor_stage_cols =  [i for i in tumor_stage_data.columns if 'group' in i]

  return tumor_stage_data, tumor_stage_cols


def zscores_for_tumor_stage_cols(cox_dict, tumor_stage_cols):
  ordered_output = []
  for group in tumor_stage_cols:
    ordered_output.extend([cox_dict[group + '-z'], cox_dict[group + '-p']])
  return ordered_output

def tumor_stage_output_header_and_format(num_initial_columns, tumor_stage_cols):
  base_formatstring = ', '.join(['{' + str(i) + '}' for i in range(num_initial_columns)])
  i += 1
  header = ''
  for var in tumor_stage_cols:
    header += ',' + var + '-zscore'
    header += ',' + var + '-pvalue'
    base_formatstring += ' ,{' + str(i) + '}, {' + str(i+1) + '}'
    i += 2
  base_formatstring += '\n'
  return header, base_formatstring


