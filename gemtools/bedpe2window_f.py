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
    wndw_row = [ r['name'],r['chrom1'],r['start1'],r['stop1'],r['chrom2'],r['start2'],r['stop2'],r['name1'],r['chrom1']] + make_window(r['start1'],r['stop1'],r['window_size']) + [r['name2'],r['chrom2']] + make_window(r['start2'],r['stop2'],r['window_size']) + [r['dist'],r['status'],r['window_size']]
    return wndw_row

# Below, try to automatically generate windows
def calc_window_size(d, m):
	if (str(d)=='na' or int(d)>=(2*int(m))):  	# If event breakpoints are distant, make large windows around breakpoint (200 kb)
		window_sz = int(m)
		status = "pass"	
	else:									# Otherwise				
		test_size = int(d) - int(m)			
		if test_size>=10000:				# Try to make windows as large as possible (lower limit 10kb), while maintaining 100 kb distance between breakpoints
			window_sz = test_size
			status="pass"
		else: # If can't do the above, then if distance between breakpoints is very small (i.e. less than 30kb), make windows 5kb
			window_sz = 10000
			status="fail"			# Otherwise, make windows 10kb
	return [window_sz,status]
	
## READ IN SV FILE + PARSE TO DESIRED FORMAT

def bedpe2window(**kwargs):
	
	if 'bedpe' in kwargs:
		sv_input = kwargs['bedpe']
	if 'window' in kwargs:
		window_size = kwargs['window']
	if 'out' in kwargs:
		outpre = kwargs['out']
	if 'mol_size' in kwargs:
		mol_size = kwargs['mol_size']

	df_sv = pd.read_table(sv_input, sep="\t", comment="#", header=None)
	df_sv.columns = ['chrom1','start1','stop1','chrom2','start2','stop2','name'] + list(df_sv.columns)[7:] 

	df_sv['name1'] = df_sv['name'].apply(lambda x: str(x) + "_1")
	df_sv['name2'] = df_sv['name'].apply(lambda x: str(x) + "_2")

	df_sv = df_sv[['name','name1','chrom1','start1','stop1','name2','chrom2','start2','stop2']]

	# Get distance between breakpoints -- return "na" if on different chromosomes
	df_sv['dist'] = df_sv.apply(lambda row: get_dist(row), axis=1)

	if str(window_size)!="None":
		df_sv['window_size'] = window_size
		df_sv['status'] = df_sv.apply(lambda row: "pass" if row['window_size']>row['dist'] else "fail", axis=1)
	elif str(mol_size)!="None":
		df_sv['window_size'] = df_sv['dist'].apply(lambda x: calc_window_size(x, mol_size)[0])
		df_sv['status'] = df_sv['dist'].apply(lambda x: calc_window_size(x, mol_size)[1])
	else:
		"Must specify either a window size or a molecule size"

	sv_wndw = df_sv.apply(lambda row: window_rows(row), axis=1)
	df_wndw = pd.DataFrame(list(sv_wndw))
	df_wndw.columns = ['name','chrom1','start1','stop1','chrom2','start2','stop2','name1','chrom1_w','start1_w','stop1_w','name2','chrom2_w','start2_w','stop2_w','dist','status','window_size']

	df_wndw.to_csv(str(outpre), sep="\t", index=False)
	return df_wndw
