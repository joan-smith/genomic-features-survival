#!/usr/bin/env python
# encoding: utf-8
'''
make pancan.py

Created by Joan Smith
on 2017-4-1.

Given a set of zscore files for each cancer type, make a pancan file, including stouffers

Copyright (c) 2018. All rights reserved.
'''

import pandas as pd
import numpy as np
import argparse
import sys
import os

sys.path.append('../common/')
import analysis


CANCER_TYPES = ['BLCA', 'BRCA', 'COADREAD', 'GBMLGG', 'HNSC', 'KIPAN', 'LIHC', 'LUAD',
                'LUSC', 'OV', 'PAAD', 'PRAD', 'SARC', 'SKCM', 'STES', 'UCEC']

CANCER_TYPES_TUMOR_STAGE = ['BLCA', 'BRCA', 'COADREAD', 'GBMLGG', 'HNSC', 'KIPAN', 'LIHC', 'LUAD',
                'LUSC', 'OV', 'PAAD', 'PRAD', 'SKCM', 'STES', 'UCEC']

CANCER_TYPES_GRADE = ['BLCA', 'GBMLGG', 'HNSC', 'KIPAN', 'LIHC', 'OV', 'PAAD', 'STES', 'UCEC']

CANCER_TYPES_ABRIDGED =  ['BLCA', 'BRCA',  'GBMLGG', 'HNSC', 'KIPAN', 'LIHC', 'LUAD',
                'LUSC', 'OV', 'PAAD', 'PRAD', 'SARC', 'STES']

ICGC_TYPES = ['COCA-CN', 'EOPC-DE', 'ESAD-UK', 'ESCA-CN', 'GACA-CN', 'LICA-FR', 'LINC-JP', 'LIRI-JP',
              'GBM-CN', 'MELA-AU', 'ORCA-IN', 'OV-AU', 'PACA-AU',
             'PACA-AU', 'PACA-CA', 'PBCA-DE', 'PRAD-UK', 'RECA-EU', 'SKCA-BR']

ICGC_CNA_TYPES = ['ESAD-UK', 'OV-AU', 'PACA-AU', 'SKCA-BR']

CBIOPORTAL_TYPES = ['BLCA_MSKCC', 'BRCA_METABRIC', 'LIHC', 'LUAD_BROAD', 'PRAD_MSKCC']

CBIOPORTAL_MUTATION_TYPES = ['BRCA_METABRIC', 'CCRC_UTOKYO', 'LUAD_BROAD', 'PRAD_CPCG', 'EGC_MUNICH']

VAF_CANCER_TYPES = ['BLCA', 'BRCA', 'HNSC', 'LIHC', 'LUAD', 'LUSC', 'PRAD', 'SARC', 'SKCM', 'UCEC']

TCGA_ADDITIONAL_CANCER_TYPES = ['ACC', 'CESC', 'CHOL', 'MESO', 'PCPG', 'TGCT', 'THCA', 'UCS', 'UVM']

MSKCC_CANCER_TYPES = ['Bladder-Urothelial-Carcinoma', 'Breast-Invasive-Ductal-Carcinoma', 'Colorectal-Cancer',
                      'Cutaneous-Melanoma',
                      'Esophagogastric-Cancer', 'Glioma', 'Lung-Adenocarcinoma', 'Lung-Squamous-Cell-Carcinoma',
                      'Pancreatic-Adenocarcinoma', 'Prostate-Adenocarcinoma', 'Renal-Cell-Carcinoma']


CANCER_TYPES_HISTOLOGICAL_SUBTYPES = ['BLCA_non-papillary',
      'BLCA_papillary', 'BRCA_ER|PR+_HER2-', 'BRCA_HER2+', 'BRCA_infiltrating_ductal_carcinoma',
      'BRCA_infiltrating_lobular_carcinoma', 'BRCA_Triple_Negative', 'COADREAD_coad',
      'COADREAD_read', 'GBMLGG_gbm', 'GBMLGG_lgg', 'HNSC_larynx', 'HNSC_tongue',
      'KIPAN_kich', 'KIPAN_kirc', 'KIPAN_kirp', 'SARC_dedifferentiated_liposarcoma',
      'SARC_leiomyosarcoma_lms', 'STES_esca', 'STES_stad']


def get_options():
  parser = argparse.ArgumentParser(description='Create pancan file for dataset')
  parser.add_argument('-i', action='store', dest='input_directory', default='.')
  parser.add_argument('-o', action='store', dest='outdir', default='NONE')
  parser.add_argument('-s', action='store', dest='suffix')

  parser.add_argument('-m', action='store_true', dest='metagene', help='Make pancan from metagene columns')
  parser.add_argument('-r', action='store_true', dest='multivariate', help='Make pancan from var-z columns')
  parser.add_argument('-g', action='store_true', dest='icgc', help='Use ICGC cancer types')
  parser.add_argument('-c', action='store_true', dest='cbioportal', help='Use CBioPortal cancer types')
  parser.add_argument('-v', action='store_true', dest='vaf', help='Use VAF cancer types')
  parser.add_argument('-t', action='store_true', dest='tcga_additional', help='Use TCGA additional cancer types')
  parser.add_argument('-k', action='store_true', dest='mskcc', help='Use MSKCC cancer types')

  ns = parser.parse_args()

  output_dir = '.'
  if ns.outdir == 'NONE':
    output_dir = ns.input_directory
  else:
    output_dir = ns.outdir

  return (output_dir, ns.input_directory, ns.suffix,
          ns.metagene, ns.multivariate, ns.icgc,
          ns.cbioportal, ns.vaf, ns.tcga_additional, ns.mskcc)

def make_path(indir, cancer, suffix):
  if os.path.isdir(os.path.join(indir, cancer)):
    return os.path.join(indir, cancer, cancer + suffix)
  else:
    return os.path.join(indir, cancer + suffix)

def make_pancan_df(outdir, indir, suffix, metagene=False, multivariate=False,
                   icgc=False, cbioportal=False, vaf=False,
                   tcga_additional=False, mskcc=False):
  pancan = {}
  zscore_header = 'zscore'
  outfile = 'pancan.csv'

  if metagene:
    print 'Calculating for metagene'
    zscore_header = 'metagene-zscore'
    outfile = 'metagene_pancan.csv'
  if multivariate:
    print 'multivariate'
    zscore_header = 'var-z'
    outfile = 'multivariate_pancan.csv'

  cancer_types = CANCER_TYPES
  if icgc:
    cancer_types = ICGC_TYPES
  elif cbioportal:
    cancer_types = CBIOPORTAL_TYPES
  elif vaf:
    cancer_types = VAF_CANCER_TYPES
  elif tcga_additional:
    cancer_types = TCGA_ADDITIONAL_CANCER_TYPES
  elif mskcc:
    cancer_types = MSKCC_CANCER_TYPES

  for cancer in cancer_types:
    make_path(indir, cancer, suffix)
    path = make_path(indir, cancer, suffix)
    df = pd.read_csv(path, index_col=0, na_values=[' NA']) # sometimes rpy2 gives back this monstrosity of a NaN value.
    if df.shape[0] > 0:
      if '\'' not in df.index[0]:
        # ADD ' to index
        df = df.reset_index()
        print 'ADD apostrophe'
        df['gene'] = '\'' + df['gene']
        df = df.set_index('gene')
      pancan[cancer] = df[zscore_header].astype(float)

  pancan_df = pd.DataFrame(pancan)
  pancan_df['stouffer unweighted'] = analysis.stouffer_unweighted(pancan_df)
  outpath = os.path.join(outdir, outfile)
  pancan_df.to_csv(outpath, index_label='gene')

def main(argv=None):
  if argv is None:
    argv = sys.argv
    outdir, indir, suffix, metagene, multivariate, icgc, cbioportal, vaf, tcga_additional, mskcc = get_options()
    make_pancan_df(outdir, indir, suffix, metagene, multivariate, icgc, cbioportal, vaf, tcga_additional, mskcc)

if __name__ == "__main__":
  main()
