#!/usr/bin/env python
# encoding: utf-8
'''
zscores.py

Calculate cna and mutation zscores for immuno therapy responses from prior papers

Created by Joan Smith
on 2017-11-04.

Copyright (c) 2018. All rights reserved.
'''

import pandas as pd
import numpy as np
import argparse
import sys
import os
import glob
import rpy2


sys.path.append('../common/')
import utilities as util
import analysis

MUTATION_PERCENT = 0.02
MUTATION_COLUMNS = {
  'MSKCC_LUAD_2015': 'SAMPLE_ID',
  'Roh_STM_2017': 'Sample',
  'SKCM_vanderbilt_MSKCC_2017': 'Tumor_Sample_Barcode',
  'Van_Allen_Science_2015': 'patient',
}

def get_options():
  parser = argparse.ArgumentParser(description='Calculate mutation and cna zscores')
  parser.add_argument('-i', action='store', dest='input_directory')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  namespace = parser.parse_args()

  return (namespace.input_directory, namespace.output_directory)

def mutations_for_gene(df, patient_col='PATIENT'):
  mutated_patients = df[patient_col].unique()
  return pd.DataFrame({'mutated': np.ones(len(mutated_patients))}, index=mutated_patients)

def prep_mutations(dir, mut, clinical):
  patient_col = MUTATION_COLUMNS[dir]
  number_barcodes_in_mutation_data = mut[patient_col].unique().size
  print 'Number of total sequenced barcodes:   ', number_barcodes_in_mutation_data

  mut['Hugo_Symbol'] = '\'' + mut['Hugo_Symbol'].astype(str)

  gene_mutation = mut.groupby(['Hugo_Symbol']).apply(mutations_for_gene, patient_col=MUTATION_COLUMNS[dir])
  gene_mutation.index.set_names(['Hugo_Symbol', 'patient'], inplace=True)
  gene_mutation = gene_mutation.reset_index()
  gene_patient_mutations = gene_mutation.pivot(index='Hugo_Symbol', columns='patient', values='mutated')
  gene_patient_mutations = gene_patient_mutations.fillna(0)

  return gene_patient_mutations.transpose()


def do_single_cancer_type_cna(name, clinical, cna, outdir):
  cna = cna.T
  if 'Chromosome' in cna.columns:
    cna = cna.drop(['Chromosome', 'Location'])
  print 'Patient count for CNAs:', cna.shape
  cnas_and_clinical = cna.join(clinical, how='inner')
  num_patients = cnas_and_clinical.shape[0]
  print 'num patients:', num_patients
  formatstring = '{0}, {1}, {2}, {3}\n'

  outfile = os.path.join(outdir, name.replace(' ', '-') + '.cnas.out.csv')
  print outfile
  with open(outfile, 'w') as out:
    out.write('gene,zscore,pvalue,num patients\n')
    for gene in cnas_and_clinical:
        if gene in ['Time', 'Censor']:
          continue
        number_non_zero = cnas_and_clinical[cnas_and_clinical[gene] != 0][gene].shape[0]
        try:
          cox_dict = analysis.do_cox(cnas_and_clinical.Time,
                                     cnas_and_clinical.Censor,
                                     cnas_and_clinical[gene], )
          out.write(formatstring.format(
                        gene, cox_dict['z'], cox_dict['p'], cox_dict['n']))
        except rpy2.rinterface.RRuntimeError as e:
          print 'Skipped ', gene, 'due to R error.'

def do_single_cancer_type_mutation(name, clinical, mutations, outdir):
  mutations_and_clinical = mutations.join(clinical, how='inner')
  formatstring = '{0}, {1}, {2}, {3}, {4}\n'
  num_patients = mutations_and_clinical.shape[0]
  print 'Number of patients in both:', num_patients
  outfile = os.path.join(outdir, name.replace(' ', '-') + '.mutations.out.csv')
  print outfile
  with open(outfile, 'w') as out:
    out.write('gene,zscore,pvalue,num mutations,num patients\n')
    for gene in mutations_and_clinical:
        if gene in ['Time', 'Censor']:
          continue
        if mutations_and_clinical[gene].sum() >= MUTATION_PERCENT*num_patients:
          try:
            cox_dict = analysis.do_cox(mutations_and_clinical.Time,
                                       mutations_and_clinical.Censor,
                                       mutations_and_clinical[gene])
            out.write(formatstring.format(
                          gene, cox_dict['z'], cox_dict['p'], mutations_and_clinical[gene].sum(), cox_dict['n']))
            mutations_and_clinical.to_csv(os.path.join(outdir, 'raw_mutations/', name + '_' +  gene + '_mutations.csv'), columns=['Time', 'Censor', gene])
          except rpy2.rinterface.RRuntimeError as e:
            print 'Skipped ', gene, 'due to R error.'

def main():
  indir, outdir = get_options()
  directories = os.listdir(indir)
  print indir, outdir
  directories = util.remove_extraneous_files(directories)
  directories.remove('output')
  for d in directories[2:]:
    print d
    cna_glob = os.path.join(indir, d, '*.cnv.*')
    print cna_glob
    cna_file = glob.glob(cna_glob)[0]
    cna = pd.read_csv(cna_file, index_col=0, sep=util.get_sep_from_filename(cna_file))

    clinical_glob = os.path.join(indir, d, '*clinical.*')
    clinical_file = glob.glob(clinical_glob)[0]
    clinical = pd.read_csv(clinical_file, sep=util.get_sep_from_filename(clinical_file), index_col=0)
    clinical = clinical[['Time', 'Censor']]

    mut_glob = os.path.join(indir, d, '*mutations*')
    mut_file = glob.glob(mut_glob)[0]
    mut = pd.read_csv(mut_file, sep=util.get_sep_from_filename(mut_file), low_memory=False)
    mutations = prep_mutations(d, mut, clinical)

    do_single_cancer_type_cna(d, clinical, cna, outdir)
    do_single_cancer_type_mutation(d, clinical, mutations, outdir)




if __name__ == "__main__":
  main()

