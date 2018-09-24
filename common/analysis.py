import numpy as np
import pandas as pd
import os
import pdb
from scipy.stats.distributions import chi2

from rpy2.robjects.packages import importr
import rpy2.robjects.numpy2ri
import rpy2.robjects.pandas2ri as pandas2ri
import rpy2.robjects as robjects
from rpy2.robjects import r

from matplotlib import pyplot

import rpy2.rinterface as ri


rpy2.robjects.numpy2ri.activate()

class KMPlot:
  def __init__(self):
    self.surv = importr('survival')

  def make_plottable_kms(self, time, percent):
    time =  list(time)
    time = [val for val in time for _ in (0, 1)]
    time.insert(0, 0)
    time.insert(-1, time[-1])

    percent = list(percent)
    percent = [val for val in percent for _ in (0, 1)]
    percent.insert(0, 1)
    percent.insert(1, 1)

    return time, percent


  def km_plot_data(self, name, time, censor, values):
    values_df = pd.DataFrame({'time': time, 'censor': censor, 'value': values}, dtype=float)
    mean_value = values_df.value.mean()
    values_df['high'] = values_df.value >= mean_value

    data = {'time': robjects.FloatVector(values_df['time']),
            'censor': robjects.IntVector(values_df['censor']),
            'high': robjects.IntVector(values_df['high'])}
    df = robjects.DataFrame(data)

    # p value
    km_diff = self.surv.survdiff(robjects.Formula('Surv(time, censor) ~ high'),
                data=df)
    chisq_ind = list(km_diff.names).index('chisq')
    pvalue = chi2.sf(km_diff[chisq_ind][0],1)


    km = self.surv.survfit(robjects.Formula('Surv(time, censor) ~ high'),
                data=df)
    summary = pandas2ri.ri2py(r.summary(km, extend=True))
    r.assign('km', km)
    r.assign('times', data['time'])
    r.assign('res', r('summary(km, times=times)'))
    cols = r('lapply(c(2:6, 8:11), function(x) res[x])')
    r.assign('cols', cols)
    km_results = r('do.call(data.frame, cols)')
    km_results = pd.DataFrame(km_results)

    low_km = km_results[km_results['strata']=='high=0']
    high_km = km_results[km_results['strata']=='high=1']

    high_time, high_percent = self.make_plottable_kms(high_km['time'], high_km['surv'])
    low_time, low_percent = self.make_plottable_kms(low_km['time'], low_km['surv'])

    high = [{'percent': i[0], 'time': i[1]} for i in zip(high_percent, high_time)]
    low = [{'percent': i[0], 'time': i[1]} for i in zip(low_percent, low_time)]

    return {'high': high, 'low': low, 'p': float('%.4g' % pvalue)}




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
  grdevices.png(file=os.path.join(outdir, name+'_km.png'), width=512, height=512)

  r.plot(km, xlab='Time', ylab='Cumulative Hazard',
                     col=robjects.StrVector(['Red', 'Blue']))
  r.legend(1000, 1, robjects.StrVector(['<= Mean', '> Mean']), lty=robjects.IntVector([1,1]), col=robjects.StrVector(['Red', 'Blue']))
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

  conf_int_ind = list(coxuh_output.names).index('conf.int')
  conf_int = coxuh_output[conf_int_ind]
  hazard_ratio = conf_int.rx('split', 'exp(coef)')[0]
  lower_conf = conf_int.rx('split', 'lower .95')[0]
  upper_conf = conf_int.rx('split', 'upper .95')[0]

  cox_dict = {
      'n': patient_count,
      'z': zscore,
      'p': pvalue,
      'hazard_ratio': hazard_ratio,
      'lower_conf': lower_conf,
      'upper_conf': upper_conf,
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

def do_multivariate_cox(time, censor, var, additional_vars, float_vars=False):
  df = pd.DataFrame({
    'time': time,
    'censor': censor,
    'var': var})
  df = df.join(additional_vars, how='inner')

  surv = importr('survival')
  r.assign('time',robjects.IntVector(np.array(df['time'])))
  r.assign('censor', robjects.IntVector(np.array(df['censor'])))
  r.assign('var', robjects.FloatVector(np.array(df['var'])))

  for i in additional_vars:
    if float_vars:
      r.assign(i, robjects.FloatVector(np.array(df[i])))
    else:
      r.assign(i, robjects.IntVector(np.array(df[i])))

  additional_vars_str = ' + '.join(additional_vars.columns)

  formula = 'Surv(time, censor) ~ var + ' + additional_vars_str
  coxuh_output = r('summary( coxph(formula = ' + formula + ', model=FALSE, x=FALSE, y=FALSE))')

  coef_ind = list(coxuh_output.names).index('coefficients')
  coeffs = coxuh_output[coef_ind]

  patient_count_ind = list(coxuh_output.names).index('n')
  patient_count = coxuh_output[patient_count_ind][0]

  var_zscore = get_zscore('var', coeffs)
  var_pvalue = get_pvalue('var', coeffs)
  hazards_dict = get_hazards('var', coxuh_output)

  cox_dict = {
      'var-n': patient_count,
      'var-z': var_zscore,
      'var-p': var_pvalue,
      }
  cox_dict.update(hazards_dict)

  for i in additional_vars.columns:
    cox_dict[i + '-z'] = get_zscore(i, coeffs)
    cox_dict[i + '-p'] = get_pvalue(i, coeffs)
    cox_dict.update(get_hazards(i, coxuh_output))

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

def get_hazards(variate, cox_output):
  conf_int_ind = list(cox_output.names).index('conf.int')
  conf_int = cox_output[conf_int_ind]
  hazard_ratio = conf_int.rx(variate, 'exp(coef)')[0]
  lower_conf = conf_int.rx(variate, 'lower .95')[0]
  upper_conf = conf_int.rx(variate, 'upper .95')[0]

  return {
    variate + '_hazard_ratio': hazard_ratio,
    variate + '_lower_conf': lower_conf,
    variate + '_upper_conf': upper_conf,
  }

def stouffer_unweighted(pancan_df):
  zscore_count = pancan_df.count(axis=1)
  sqrt_zscore_count = zscore_count.apply(np.sqrt)
  zscore_sums = pancan_df.sum(axis=1, skipna=True)
  unweighted_stouffer = zscore_sums / sqrt_zscore_count
  return unweighted_stouffer

