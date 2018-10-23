import sys
import os
import __main__ as main
import argparse
import ast
import pandas as pd
import pysam
import numpy as np


MIN_MAPQ = 0


## DEFINE FUNCTION TO CREATE WINDOWS AROUND BREAKPOINTS

def make_window(s,e,w):
	cur_size = e-s
	adj_val = (w-cur_size)/2
	adj_val = int(round(adj_val,0))
	new_start = s - adj_val
	new_end = e + adj_val
	return [new_start,new_end]

## DEFINE FUNCTION TO OBTAIN BARCODES FROM BAM FILE FOR SPECIFIC REGIONS

def get_barcode_ids(bam_in, chrom, start, end, min_mapq):
	bcs = []
	for r in bam_in.fetch(chrom, start, end):
	  if r.mapq >= min_mapq:
		  if r.has_tag("BX"):
			  bc_id=r.get_tag("BX")
			  bcs.append(bc_id)
	return list(bcs)

#def count_bcs(sv_in, bam_file, full_w_size, small_w_size):
def count_bcs(full_w_size=500000, small_w_size=1000,bc_subset='shared',sv_n="None",outpre="out",**kwargs):

	if 'in_window' in kwargs:
		small_w_size = kwargs['in_window']
	if 'out_window' in kwargs:
		full_w_size = kwargs['out_window']	
	if 'sv' in kwargs:
		sv_input = kwargs['sv']
	if 'bam' in kwargs:
		bam_input = kwargs['bam']
	if 'sv_name' in kwargs:
		sv_n = kwargs['sv_name']
	if 'bcs' in kwargs:
		bc_subset = kwargs['bcs']
	if 'out' in kwargs:
		outpre = kwargs['out']

	full_w_size = int(full_w_size)  # -l #500000
	small_w_size = int(small_w_size)
	
	sv_df = pd.read_table(sv_input, sep="\t")
	sv_df[['name']] = sv_df[['name']].astype(str)
	
	svs_in_df = list(set(sv_df['name'].tolist()))
	
	#print locals()
	#print str(sv_n)
	
	sv_list = list(set(sv_n.split(",")))
	sv_list = [str(s) for s in sv_list]
	
	for s in sv_list:
		if not s in svs_in_df:
			print str(s) + " not present in table"
			sys.exit()
	
	print sv_list
	
	if str(sv_n)=="None":
		sys.exit()
	else:
		sv_df = sv_df.loc[sv_df['name'].isin(sv_list)]
	
	print sv_df
	
	# Make list of SV-specific barcodes
		
	if bc_subset=="shared":
		bc_list = sv_df[['bc_overlap_id']].sum()
		#bc_list = ast.literal_eval(row['bc_overlap_id'])
		if len(bc_list)<1:
			print "No shared barcodes -- exiting"
			sys.exit()
		
	elif bc_subset=="all":
		sv_df['all_bcs'] = sv_df.apply(lambda row: tuple(set(row['bc_1_id'] + row['bc_2_id'])), axis=1)
		bc_list = sv_df[['all_bcs']].sum()
		#bc_list = list(set(ast.literal_eval(row['bc_1_id']) + ast.literal_eval(row['bc_2_id'])))
		if len(bc_list)<1:
			print "No barcodes -- exiting"
			sys.exit()
	else:
		print "bcs -- must be either 'shared' or 'all'"
		sys.exit()
	
	print bc_list
	print len(bc_list)
	sys.exit()
	
	bam_open = pysam.Samfile(bam_input)

	for s in sv_list:
		sv_df_cur = sv_df.loc[sv_df['name']==str(s)]

		sv_df1 = sv_df_cur[['name','name1','chrom1','start1','stop1','bc_1_id','bc_2_id','bc_overlap_id']]
		sv_df2 = sv_df_cur[['name','name2','chrom2','start2','stop2','bc_1_id','bc_2_id','bc_overlap_id']]
		full_names = ['id','name','chrom','start','stop','bc_1_id','bc_2_id','bc_overlap_id']
		sv_df1.columns = full_names
		sv_df2.columns = full_names
		sv_df_full = pd.concat([sv_df1,sv_df2])


		w_start_list = [x[0] for x in sv_df_full.apply(lambda row: make_window(row['start'],row['stop'], full_w_size), axis=1)]
		w_stop_list = [x[1] for x in sv_df_full.apply(lambda row: make_window(row['start'],row['stop'], full_w_size), axis=1)]
		sv_df_full['w_start'] = w_start_list
		sv_df_full['w_stop'] = w_stop_list



		df_list = []
		for index,row in sv_df_full.iterrows():

			df_name = str(row['name']) + "_df"
			# Create data frame of 1kb windows
			start = int(row['w_start'])
			stop = int(row['w_stop'])
			window_start = np.arange(start, start+full_w_size, small_w_size+1)
			window_end = window_start+small_w_size
			df_name = pd.DataFrame([window_start, window_end]).transpose()
			df_name.columns = ['window_start','window_end']
			df_name['chrom'] = str(row['chrom'])
			df_name['id'] = str(row['id'])
			df_name['name'] = str(row['name'])
			df_name = df_name[['id','name','chrom','window_start','window_end']]

			# Make a list of barcodes in each of the 1kb windows
			region_bcs = []
			for i,r in df_name.iterrows():  
				region_bc_list = get_barcode_ids(bam_open, r['chrom'], r['window_start'], r['window_end'], MIN_MAPQ)
				region_bcs.append(region_bc_list)


			#bc_list = ast.literal_eval(row['bc_overlap_id'])

			# For each SV-specific barcode, count the number of times it occurs in each region
			for bc in bc_list:
				df_name[bc] = [x.count(bc) for x in region_bcs]

			df_list.append(df_name)
		
	# Write output to file
	df_full = pd.concat(df_list)
	df_full.to_csv(str(outpre), sep="\t", index=False)
	return df_full

