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
	cm = "hot_r"
	fig = plt.figure(figsize=(3*len(xlabs), len(ylabs)))
	sns.set(style="white", font_scale = 2)
	sns.heatmap(data, xticklabels = xlabs, yticklabels = ylabs, cmap=cm, cbar_kws = {"fraction":0.5, "shrink":0.5})
	plt.xticks(rotation=90)
	fig.savefig(fname = fname, bbox_inches = 'tight', pad_inches = 1)

def plotClustermap(enr, pvals, fname, xlabs, ylabs):
	cm = "hot_r"
	sns.set(style="white", font_scale = 2)
	fig0 = sns.clustermap(enr, figsize = (5*len(xlabs), len(ylabs)), xticklabels = xlabs, yticklabels = ylabs, cmap=cm, cbar_kws = {"fraction":0.5, "shrink":0.8})
	row_order = fig0.dendrogram_row.reordered_ind
	col_order = fig0.dendrogram_col.reordered_ind
	pvals = pvals[:, col_order][row_order]
	fig = sns.clustermap(enr, figsize = (5*len(xlabs), len(ylabs)), xticklabels = xlabs, yticklabels = ylabs, cmap=cm, cbar_kws = {"fraction":0.5, "shrink":0.8}, annot = pvals)
	fig.savefig(fname = fname+".enrichment.png", bbox_inches = 'tight', pad_inches = 1)


def readLDSCfiles(fnames):
	alldat = [pd.read_csv(f, delim_whitespace=True) for f in fnames]
	enr = np.concatenate([d[['Enrichment']].fillna(0).values for d in alldat], 1)
	pvals = -np.log(np.concatenate([d[['Enrichment_p']].fillna(1).values for d in alldat], 1))
	enr_er = np.concatenate([d[['Enrichment_std_error']].fillna(0).values for d in alldat], 1)
	op = np.argsort(-np.mean(pvals,1))
	pvals = pvals[op]
	enr = enr[op]
	enr_er = enr_er[op]
	ylabs = [l.split("_0")[0] for l in alldat[0]['Category'][op]]
	return alldat, enr, pvals, enr_er, ylabs

# Code to run
parser = argparse.ArgumentParser(description='Plot results of LDSC',formatter_class=argparse.ArgumentDefaultsHelpFormatter)
parser.add_argument('--results', nargs=1, type=str, default="", help='comma separated list of result files from LDSC')
parser.add_argument('--labels', nargs=1, type=str, default="", help='comma separated list of labels for each LDSC analysis')
parser.add_argument('--outpref', nargs='?', type=str, default="", help='output file preference')
parser.add_argument('--pthresh', nargs='?', type=float, default=0.05, help='pvalue threshold for plotting')
parser.add_argument('--ethresh', nargs='?', type=float, default=np.inf, help='enrichment std error threshold for plotting')

args = parser.parse_args()

fnames = args.results[0].replace(" ", "").split(",")

labels = args.labels[0].replace(" ", "").split(",")
alldat, enr, pvals, enr_er, ylabs = readLDSCfiles(fnames)

nominal = np.max(pvals,1) > -np.log(args.pthresh)
pvals = pvals[nominal,:]
enr = enr[nominal,:]
enr_er = enr_er[nominal,:]
ylabs = [y for i, y in enumerate(ylabs) if nominal[i]]

abv_error = np.max(enr_er,1) < args.ethresh
pvals = pvals[abv_error,:]
enr = enr[abv_error,:]
enr_er = enr_er[abv_error,:]
ylabs = [y for i, y in enumerate(ylabs) if abv_error[i]]


if np.shape(enr)[1] == 1:
	plotHeatmap(enr, args.outpref+"enrichment.png", args.labels, ylabs)
	plotHeatmap(pvals, args.outpref+"pvalues.png", args.labels, ylabs)
else:
	plotClustermap(enr, pvals, args.outpref, labels, ylabs)