# Packages
import pandas as pd
import numpy as np
import matplotlib as mpl
mpl.use('Agg')
import matplotlib.pyplot as plt
import seaborn as sns
from pprint import pprint
import argparse

# Functions
def plotHeatmap(data, fname, xlabs, ylabs):
	cm = "YlOrRd"
	fig = plt.figure(figsize=(3*len(xlabs), len(ylabs)))
	sns.set(style="white", font_scale = 2)
	sns.heatmap(data, xticklabels = xlabs, yticklabels = ylabs, cmap=cm, cbar_kws = {"fraction":0.5, "shrink":0.5})
	plt.xticks(rotation=90)
	fig.savefig(fname = fname, bbox_inches = 'tight', pad_inches = 1)

def plotClustermap(enr, pvals, fname, xlabs, ylabs):
	cm = "YlOrRd"
	sns.set(style="white", font_scale = 4)
	fig = sns.clustermap(enr, figsize = (5*len(xlabs), len(ylabs)), xticklabels = xlabs, yticklabels = ylabs, cmap=cm, cbar_kws = {"fraction":0.5, "shrink":0.5})
	row_order = fig.dendrogram_row.reordered_ind
	col_order = fig.dendrogram_col.reordered_ind
	fig.savefig(fname = fname+".enrichment.png", bbox_inches = 'tight', pad_inches = 1)
	pvals = pvals[:, col_order][row_order]
	fig2 = plt.figure(figsize=(5*len(xlabs), len(ylabs)))
	sns.heatmap(pvals, xticklabels = [xlabs[c] for c in col_order], yticklabels = [ylabs[r] for r in row_order], cmap=cm, cbar_kws = {"fraction":0.5, "shrink":0.5})
	plt.xticks(rotation=90)
	fig2.savefig(fname = fname+".pvalues.png", bbox_inches = 'tight', pad_inches = 1)

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
parser.add_argument('--pthresh', nargs='?', type=float, default=0.05, help='pvalue threshold for plotting')

args = parser.parse_args()

fnames = args.results[0].replace(" ", "").split(",")

labels = args.labels[0].replace(" ", "").split(",")
alldat, enr, pvals, ylabs = readLDSCfiles(fnames)

nominal = np.min(pvals,1) > -np.log(args.pthresh)
pvals = pvals[nominal,:]
enr = enr[nominal,:]
ylabs = [y for i, y in enumerate(ylabs) if nominal[i]]


if np.shape(enr)[1] == 1:
	plotHeatmap(enr, args.outpref+"enrichment.png", args.labels, ylabs)
	plotHeatmap(pvals, args.outpref+"pvalues.png", args.labels, ylabs)
else:
	plotClustermap(enr, pvals, args.outpref, labels, ylabs)