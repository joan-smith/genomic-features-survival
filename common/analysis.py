import numpy as np
import pandas as pd
import os

from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
import rpy2.robjects as robjects
from rpy2.robjects import r


import rpy2.rinterface as ri

rpy2.robjects.numpy2ri.activate()

def do_km(name, time, censor, split, outdir):
  """Given three clean (pre-processed) lists, make a kmplot of the data, and save it to outdir"""
  data = {'time': robjects.IntVector(np.array(time)),
          'censor': robjects.IntVector(np.array(censor)),
          'split': robjects.IntVector(np.array(split))}
  df = robjects.DataFrame(data)

  surv = importr('survival')
  grdevices = importr('grDevices')
  km = surv.survfit(robjects.Formula('Surv(time, censor) ~ split'),
              data=df)
  grdevices.png(file=os.path.join(outdir, name+'_km_from_mutations.png'), width=512, height=512)
  surv.plot_survfit(km, xlab='Time', ylab='Cumulative Hazard')
  grdevices.dev_off()

def do_cox(time, censor, split, float_time=False):
  """Given three clean (pre-processed) lists, do a simple cox analysis."""
  surv = importr('survival')
  if float_time:
    r.assign('time',robjects.FloatVector(np.array(time)))
  else:
    r.assign('time',robjects.IntVector(np.array(time)))
  r.assign('censor', robjects.IntVector(np.array(censor)))
  r.assign('split', robjects.FloatVector(np.array(split)))

  coxuh_output = r('summary( coxph(formula = Surv(time, censor) ~ split, model=FALSE, x=FALSE, y=FALSE))')

  coef_ind = list(coxuh_output.names).index('coefficients')
  coeffs = coxuh_output[coef_ind]

  patient_count_ind = list(coxuh_output.names).index('n')
  patient_count = coxuh_output[patient_count_ind][0]
  zscore = get_zscore('split', coeffs)
  pvalue = get_pvalue('split', coeffs)

  cox_dict = {
      'n': patient_count,
      'z': zscore,
      'p': pvalue,

      }
  return cox_dict

def do_metagene_cox(time, censor, split, metagene):
  df = pd.DataFrame({
    'time': time,
    'censor': censor,
    'split': split})
  df = df.join(metagene, how='inner')

  surv = importr('survival')
  r.assign('time',robjects.IntVector(np.array(df['time'])))
  r.assign('censor', robjects.IntVector(np.array(df['censor'])))
  r.assign('split', robjects.FloatVector(np.array(df['split'])))
  r.assign('metagene', robjects.FloatVector(np.array(df['metagene'])))

  coxuh_output = r('summary( coxph(formula = Surv(time, censor) ~ split + metagene, model=FALSE, x=FALSE, y=FALSE))')

  coef_ind = list(coxuh_output.names).index('coefficients')
  coeffs = coxuh_output[coef_ind]

  patient_count_ind = list(coxuh_output.names).index('n')
  patient_count = coxuh_output[patient_count_ind][0]

  split_zscore = get_zscore('split', coeffs)
  split_pvalue = get_pvalue('split', coeffs)
  metagene_zscore = get_zscore('metagene', coeffs)
  metagene_pvalue = get_pvalue('metagene', coeffs)

  cox_dict = {
      'n': patient_count,
      'z': split_zscore,
      'p': split_pvalue,
      'metagene-z': metagene_zscore,
      'metagene-p': metagene_pvalue,
      }
  return cox_dict

# These two functions exist because in some cases rpy returns the silly type NARealType instead
# of a float with value nan. This str(NARealType) == NA (not nan), which makes for weird
# inconsistencies.
def get_zscore(variate, coeffs):
  zscore = coeffs.rx(variate, 'z')[0]
  if type(zscore) == ri.NARealType:
    zscore = np.nan
  return zscore

def get_pvalue(variate, coeffs):
  pvalue = coeffs.rx(variate, 'Pr(>|z|)')[0]
  if type(pvalue) == ri.NARealType:
    pvalue = np.nan
  return pvalue


def stouffer_unweighted(pancan_df):
  zscore_count = pancan_df.count(axis=1)
  sqrt_zscore_count = zscore_count.apply(np.sqrt)
  zscore_sums = pancan_df.sum(axis=1, skipna=True)
  unweighted_stouffer = zscore_sums / sqrt_zscore_count
  return unweighted_stouffer

