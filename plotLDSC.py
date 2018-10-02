# Packages
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
from pprint import pprint
import argparse

# Functions
def plotHeatmap(data, fname, xlabs, ylabs):
	fig = plt.figure(figsize=(3*len(xlabs), len(ylabs)))
	sns.set(style="white", font_scale = 2)
	sns.heatmap(data, xticklabels = xlabs, yticklabels = ylabs, cmap="inferno", cbar_kws = {"fraction":0.5, "shrink":0.5})
	plt.xticks(rotation=90)
	fig.savefig(fname = fname, bbox_inches = 'tight', pad_inches = 1)

def plotClustermap(data, fname, xlabs, ylabs):
	sns.set(style="white", font_scale = 2)
	fig = sns.clustermap(data, figsize = (3*len(xlabs), len(ylabs)), xticklabels = xlabs, yticklabels = ylabs, cmap="inferno", cbar_kws = {"fraction":0.5, "shrink":0.5})
	fig.savefig(fname = fname, bbox_inches = 'tight', pad_inches = 1)

def readLDSCfiles(fnames):
	alldat = [pd.read_csv(f, delim_whitespace=True) for f in fnames]
	enr = np.concatenate([d[['Enrichment']].fillna(0).values for d in alldat], 1)
	pvals = -np.log(np.concatenate([d[['Enrichment_p']].fillna(1).values for d in alldat], 1))
	op = np.argsort(-np.mean(pvals,1))
	pvals = pvals[op]
	enr = enr[op]
	ylabs = [l.split("_0")[0] for l in alldat[0]['Category'][op]]
	return alldat, enr, pvals, ylabs

# Code to run
parser = argparse.ArgumentParser(description='Plot results of LDSC',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--results', nargs=1, type=str, default="", help='comma separated list of result files from LDSC')
parser.add_argument('--labels', nargs=1, type=str, default="", help='comma separated list of labels for each LDSC analysis')
parser.add_argument('--outpref', nargs='?', type=str, default="", help='output file preference')

args = parser.parse_args()

fnames = args.results[0].replace(" ", "").split(",")

labels = args.labels[0].replace(" ", "").split(",")
alldat, enr, pvals, ylabs = readLDSCfiles(fnames)

if np.shape(enr)[1] == 1:
	plotHeatmap(enr, args.outpref+".png", args.labels, ylabs)
else:
	plotClustermap(enr, args.outpref+".png", labels, ylabs)