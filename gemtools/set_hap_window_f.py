import os
import sys
import __main__ as main
import pandas as pd
import numpy as np
import pysam
import vcf


## DEFINE FUNCTIONS TO CREATE WINDOWS AROUND BREAKPOINTS

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Parse SV input file to desired format                                       #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

def get_dist(r):
	if r['chrom1']!=r['chrom2']:
		dist="na"
	else:
		dist = r['start2'] - r['stop1']
	return dist

def make_window(s,e,w):
	cur_size = int(e)-int(s)
	adj_val = (int(w)-cur_size)/2
	adj_val = int(round(adj_val,0))
	new_start = max(0,s - adj_val)
	new_end = e + adj_val
	return [new_start,new_end]

def window_rows(r):
	wndw_row = [ r['name'],r['chrom1'],r['start1'],r['stop1'],r['chrom2'],r['start2'],r['stop2'],r['name1'],r['chrom1']] + make_window(r['start1'],r['stop1'],r['window_size']) + [r['name2'],r['chrom2']] + make_window(r['start2'],r['stop2'],r['window_size']) + [r['dist'],r['window_size']]
	return wndw_row

	
## READ IN SV FILE + PARSE TO DESIRED FORMAT

def set_hap_window(**kwargs):
	
	if 'bedpe' in kwargs:
		sv_input = kwargs['bedpe']
	if 'window' in kwargs:
		window_size = kwargs['window']
	else:
		window_size = 1000000
	if 'out' in kwargs:
		outpre = kwargs['out']

	df_sv = pd.read_table(sv_input, sep="\t", comment="#", header=None)
	df_sv.columns = ['chrom1','start1','stop1','chrom2','start2','stop2','name'] + list(df_sv.columns)[7:] 

	df_sv['name1'] = df_sv['name'].apply(lambda x: str(x) + "_1")
	df_sv['name2'] = df_sv['name'].apply(lambda x: str(x) + "_2")

	df_sv = df_sv[['name','name1','chrom1','start1','stop1','name2','chrom2','start2','stop2']]

	# Get distance between breakpoints -- return "na" if on different chromosomes
	df_sv['dist'] = df_sv.apply(lambda row: get_dist(row), axis=1)

	df_sv['window_size'] = window_size
	
	sv_wndw = df_sv.apply(lambda row: window_rows(row), axis=1)
	df_wndw = pd.DataFrame(list(sv_wndw))
	df_wndw.columns = ['name','chrom1','start1','stop1','chrom2','start2','stop2','name1','chrom1_w','start1_w','stop1_w','name2','chrom2_w','start2_w','stop2_w','dist','window_size']

	df_wndw.to_csv(str(outpre), sep="\t", index=False)
	return df_wndw
