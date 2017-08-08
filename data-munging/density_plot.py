#!/usr/bin/env python
# encoding: utf-8
'''
density_plot.py

Created by Joan Smith
on 2017-7-15

Given a file, make a density plot

Copyright (c) 2017 . All rights reserved.
'''

import argparse
import sys
import os
import glob
import pandas as pd
import numpy as np
from scipy.stats import gaussian_kde
from scipy.stats import linregress
import matplotlib.pyplot as plt
import matplotlib as mpl


def get_options():
  parser = argparse.ArgumentParser(description='Produce intermediate file for analysis.')
  parser.add_argument('-i', action='store', dest='input')
  parser.add_argument('-o', action='store', dest='output_directory', default='.')
  namespace = parser.parse_args()
  return (namespace.input, namespace.output_directory)

def main():
  infile, outdir = get_options()
  font = {'fontname':'Arial',
          'size': 27}
  density_plot_data = pd.read_csv(infile, header=None, dtype=None, low_memory=True)
  title = density_plot_data.iloc[1][0].split(':', 1)[1].strip()
  xaxis = density_plot_data.iloc[2][0].split(':', 1)[1].strip()
  yaxis = density_plot_data.iloc[3][0].split(':', 1)[1].strip()
  lims = density_plot_data.iloc[5][0].split(':', 1)[1].strip().split(' ')
  ylims = None
  if 'Y' in str(density_plot_data.iloc[6][0]):
    ylims = density_plot_data.iloc[6][0].split(':')[1].strip().split(' ')
    ylims = [float(ylims[0]), float(ylims[2])]
    print ylims
  lims = [int(lims[0]), int(lims[2])]
  print lims
  print title
  print xaxis
  print yaxis

  density_plot_data = density_plot_data[density_plot_data[1].str.contains('RNASeq') == False]
  density_plot_data = density_plot_data.dropna(subset=[1,2], how='any')
  x = density_plot_data[1][1:].astype(float)
  y = density_plot_data[2][1:].astype(float)

  # calculate point density
  xy = np.vstack([x,y])
  z = gaussian_kde(xy)(xy)

  mpl.rcParams['xtick.labelsize'] = 22
  mpl.rcParams['ytick.labelsize'] = 22

  fig, ax = plt.subplots()
  ax.scatter(x, y, c=z, s=25, edgecolor='')
  if ylims:
    ax.set_ylim(ylims)
  else:
    ax.set_ylim(lims)
  ax.set_xlim(lims)
  plt.title(title, y=1.05, **font)
  plt.xlabel(xaxis, labelpad=20, **font)
  plt.ylabel(yaxis, **font)
  plt.savefig(infile.split('.')[0] + '.unsorted_output.png', bbox_inches='tight', pad_inches=0.5)
  plt.show()


if __name__ == "__main__":
  main()
