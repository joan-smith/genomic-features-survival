import pandas as pd
import numpy as np
import sys
import os

from matplotlib import pyplot
pyplot.switch_backend('agg')

PRIMARY_TUMOR_PATIENT_ID_REGEX = '^.{4}-.{2}-.{4}-01.*'
METASTASIS =  '^.{4}-.{2}-.{4}-06.*'

SHORTEN_PATIENT_REGEX = '^(.{4}-.{2}-.{4}).*'

PROCUREMENT_DAYS_BY_CANCER_TYPE = {
        'BLCA': 'patient.biospecimen_cqcf.days_to_sample_procurement',
        'BRCA': 'patient.biospecimen_cqcf.days_to_sample_procurement',
        'COADREAD': 'patient.biospecimen_cqcf.days_to_sample_procurement',
        'GBMLGG': None,
        'HNSC': 'patient.biospecimen_cqcf.days_to_sample_procurement',
        'KIPAN': 'patient.biospecimen_cqcf.days_to_sample_procurement',
        'LIHC': 'patient.biospecimen_cqcf.days_to_sample_procurement',
        'LUAD': 'patient.biospecimen_cqcf.days_to_sample_procurement',
        'LUSC': 'patient.biospecimen_cqcf.days_to_sample_procurement',
        'OV': 'patient.biospecimen_cqcf.days_to_sample_procurement',
        'PAAD': 'patient.biospecimen_cqcf.days_to_sample_procurement',
        'PRAD': 'patient.biospecimen_cqcf.days_to_sample_procurement',
        'SARC': 'patient.biospecimen_cqcf.days_to_sample_procurement',
        'SKCM': 'patient.biospecimen_cqcf.days_to_sample_procurement',
        'STES': 'patient.tumor_samples.tumor_sample.days_to_sample_procurement',
        'UCEC': 'patient.tumor_samples.tumor_sample.days_to_sample_procurement',
        }

PRAD_ENDPOINT = 'patient.days_to_first_biochemical_recurrence'

def get_sep_from_filename(filename):
  if '.csv' in filename:
    return ','
  if '.txt' in filename:
    return '\t'

  print 'Unknown separator for filename', filename
  sys.exit(1)


def maybe_clear_non_01s(df, patient_column, cancer_type):
  """ Returns a data frame with appropriately cleaned up patients.

  If the input file is an SKCM (melanoma) file, include both primary tumor and metatastasis data.
  Otherwise only include primary tumor data.

  This utility should be called before doing other processing on the data.

  Args:
    df: genomic features data
    patient_column: string that indexes patients in the df
    cancer_type: determines whether file is SKCM

  Returns:
    df: same dataframe, with inelligble rows removed
  """
  non_01_barcodes = df[~df[patient_column].str.contains(PRIMARY_TUMOR_PATIENT_ID_REGEX)][patient_column]
  print 'Number of barcodes ending in not -01: ', non_01_barcodes.unique().shape[0]

  if non_01_barcodes.count() > 0:
    if non_01_barcodes.unique().shape[0] < 10:
      print 'WARN WARN WARN: non-01 barcodes ', non_01_barcodes.unique()
    else:
      print 'WARN WARN WARN: many non 01-barcodes'

    if 'SKCM' in cancer_type:
      return process_SKCM_df(df, patient_column)
    else:
      print 'Dropping non-01 barcodes'
      df = df.drop(non_01_barcodes.index)
  return df

def add_identifier_column(df, patient_column):
  """Given a data frame and a patient column, add 'identifier' column with the shortened, standard patient id."""
  shortened_patients = df[patient_column].str.extract(SHORTEN_PATIENT_REGEX, expand=False)
  df['identifier'] = shortened_patients
  return df

def get_clinical_data(clinical_file):
  """Given a clincal file, return a dataframe with columns: time and censor, indexed by the short patient id"""

  cancer_type = get_cancer_type(clinical_file)


  identifiers = []
  days_to_death = []
  days_to_last_followup = []
  days_to_procurement = []

  procurement_label = PROCUREMENT_DAYS_BY_CANCER_TYPE[cancer_type]
  endpoint_label = 'patient.days_to_death'
  if cancer_type == 'PRAD':
    endpoint_label = PRAD_ENDPOINT
  print cancer_type, endpoint_label, procurement_label

  with open(clinical_file, 'rU') as f:
    for line in f.readlines():
      if 'patient.bcr_patient_barcode' in line:
        identifiers = [i.upper().strip() for i in line.split('\t')[1:]]
      elif endpoint_label in line:
        days_to_death = convert_to_num(line.split('\t')[1:])
      elif 'patient.days_to_last_followup' in line:
        days_to_last_followup = convert_to_num(line.split('\t')[1:])
      elif procurement_label and procurement_label in line:
        days_to_procurement = convert_to_num(line.split('\t')[1:])
        days_to_procurement = np.nan_to_num(days_to_procurement)

  if len(days_to_procurement) == 0:
    days_to_procurement = [0]*len(days_to_death)

  data = zip(identifiers, days_to_death, days_to_last_followup, days_to_procurement)
  orig_identifiers = identifiers

  # This is where clinical data gets cleaned up from being full of silliness to being usable.
  del_ind = []
  for i,d, in enumerate(data):
    if (np.isnan(d[1]) and np.isnan(d[2])):
      # days_to_death and days_to_last followup are both nan
      del_ind.append(i)

  bad_identifiers = [data[i][0] for i in del_ind] # grab the bad data so it can be spot checked.

  data = [d for i,d in enumerate(data) if i not in del_ind]
  return make_clinical_df(data)

def get_cancer_type(f):
  return os.path.basename(f).split('.')[0]

def remove_extraneous_files(files):
  if '.DS_Store' in files:
    files.remove('.DS_Store')
  return [f for f in files if not 'pancan' in f]


def make_histogram(histogram_data, outdir, permutation_count=None, show=False):
  minimum = min(histogram_data)
  maximum = max(histogram_data)

  hist_counts, bins, _ = pyplot.hist(histogram_data, bins=100)
  hist_counts = np.append(hist_counts, np.nan)
  bucketed_df = pd.DataFrame({'bins': bins, 'count': hist_counts})

  with open(os.path.join(outdir, 'pancan_histogram_data.csv'), 'w') as out:
    out.write('permutations,' + str(permutation_count) + '\n')
    out.write('min,' + str(minimum) + '\n')
    out.write('max,' + str(maximum) + '\n')
    bucketed_df.to_csv(out, index=False)

  pyplot.title(outdir)
  pyplot.savefig(os.path.join(outdir, 'pancan_histogram.png'))

  print outdir
  print minimum
  print maximum
  if show:
    pyplot.show()



"""Private Helpers"""
def make_clinical_df(processed_clinical_data):
  """Private. Transforms the processed clinical data with days to death/last followup to a DF with time/censor
  Args: processed_clinical_data: list of tuples
    tuple[0]: identifier
    tuple[1]: days_to_death
    tuple[2]: days_to_last_followup
    tuple[3]: days_to_procurement
  """
  times = []
  censors = []
  identifiers = []
  for d in processed_clinical_data:
    censor =  0 if np.isnan(d[1]) else 1
    time_since_surgery = d[1] - d[3]
    time = d[1] if censor else d[2]
    time = time - d[3] # subtract days since procurement, to get time_since_surgery
    times.append(time)
    censors.append(censor)
    identifiers.append(d[0])

  df = pd.DataFrame({'time': times, 'censor': censors}, columns=['time', 'censor'], index=identifiers)
  negative_times = df[df['time'] <= 0].index
  return df.drop(negative_times)


def convert_to_num(split_text):
  """Private. Convert a list of (mostly) numbers, convert each one by one to an int, or replace with nan."""
  nums = []
  for i in split_text:
    try:
      nums.append(int(i))
    except ValueError:
      nums.append(np.nan)
  return nums

def process_SKCM_df(df, patient_column):
  """ Returns the SKCM data frame with appropriate rows removed:
  1. remove a row if it is neither a primary tumor nor a metastasis
  2. remove a row if it is a metastasis, but there exists data for the primary tumor from the same patient
  """
  print 'Processing for SKCM'
  # Remove all non-01 and non-06 rows
  non_01_barcodes = df[~df[patient_column].str.contains(PRIMARY_TUMOR_PATIENT_ID_REGEX)][patient_column]
  non_01_non_06 = non_01_barcodes[~non_01_barcodes.str.contains(METASTASIS)]
  df = df.drop(non_01_non_06.index).copy()

  # This is complicated and a bit gross. We need to find all the metastasis samples that also have
  # primary tumor samples, and remove them.
  df = add_identifier_column(df, patient_column) # so we can compare easily
  # pull out the short ids of metastasis samples
  metastasis_ids = df[df[patient_column].str.contains(METASTASIS)]
  metastasis_uniques = set(metastasis_ids['identifier'].unique())
  # similarly, pull out the short ids of primary samples
  primary_tumor_ids = df[patient_column].str.contains(PRIMARY_TUMOR_PATIENT_ID_REGEX)
  primary_uniques = set(df[primary_tumor_ids]['identifier'].unique())
  # ew. Find the metastasis samples that aren't also in primary. Then take the opposite,
  # getting the metastasis sampes that *are* also in primary samples
  metastasis_with_primaries = metastasis_uniques - (metastasis_uniques -  primary_uniques)
  index_of_metastasis_with_primaries = metastasis_ids[metastasis_ids['identifier'].isin(metastasis_with_primaries)].index
  # yay. 5 alarming lines later, drop the duplicate metastasis.
  df = df.drop(index_of_metastasis_with_primaries)
  print "With only non-duplicated metastasis", len(df['identifier'].unique())
  df = df.drop('identifier', 1) # hacks. Drop the identifier index so we can add it back later :|
  return df

