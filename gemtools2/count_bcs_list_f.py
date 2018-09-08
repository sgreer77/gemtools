import sys
import os
import __main__ as main
import argparse
import ast
import pandas as pd
import pysam
import numpy as np


MIN_MAPQ = 0

## DEFINE FUNCTION TO OBTAIN BARCODES FROM BAM FILE FOR SPECIFIC REGIONS

def get_barcode_ids(bam_in, chrom, start, end, min_mapq):
	bcs = []
	for r in bam_in.fetch(chrom, start, end):
	  if r.mapq >= min_mapq:
		  if r.has_tag("BX"):
			  bc_id=r.get_tag("BX")
			  bcs.append(bc_id)
	return list(bcs)


def count_bcs_list(outpre='out',region_subset='None',small_w_size=1000,bc_subset='None',**kwargs):

	if 'region' in kwargs:
		region_subset = kwargs['region']
	if 'in_window' in kwargs:
		small_w_size = kwargs['in_window']
	if 'bam' in kwargs:
		bam_input = kwargs['bam']
	if 'bcs' in kwargs:
		bc_subset = kwargs['bcs']
	if 'out' in kwargs:
		outpre = kwargs['out']

	small_w_size = int(small_w_size)

	if str(bc_subset)=="None":
		sys.exit()
	if str(region_subset)=="None":
		sys.exit()

	if os.path.isfile(str(bc_subset)):
		with open(bc_subset) as f:
			bc_list = f.read().splitlines()
			print bc_list
	else:
		sys.exit()
		
	bam_open = pysam.Samfile(bam_input)

	# Create data frame of 1kb windows
	
	region_list = region_subset.split(';')
	df_r_list = []
	for reg in region_list:
		start = int(str(reg).split(',')[1])
		stop = int(str(reg).split(',')[2])
		full_w_size = stop-start
		window_start = np.arange(start, start+full_w_size, small_w_size+1)
		window_end = window_start+small_w_size
		df_out = pd.DataFrame([window_start, window_end]).transpose()
		df_out.columns = ['window_start','window_end']
		df_out['chrom'] = str(reg).split(',')[0]
		df_out = df_out[['chrom','window_start','window_end']]
		df_r_list.append(df_out)
	
	df_r_out = pd.concat(df_r_list)

	# Make a list of barcodes in each of the 1kb windows
	region_bcs = []
	for i,r in df_r_out.iterrows():  
		region_bc_list = get_barcode_ids(bam_open, r['chrom'], r['window_start'], r['window_end'], MIN_MAPQ)
		region_bcs.append(region_bc_list)

	# For each SV-specific barcode, count the number of times it occurs in each region
	for bc in bc_list:
		df_r_out[bc] = [x.count(bc) for x in region_bcs]
	
# Write output to file

	df_r_out.to_csv(str(outpre), sep="\t", index=False)
	return df_r_out

