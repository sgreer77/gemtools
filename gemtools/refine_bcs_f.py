import sys
import os
import __main__ as main
import argparse
import ast
import pandas as pd
import pysam
import numpy as np


MIN_MAPQ = 0

def get_barcode_ids(bam_in, chrom, start, end, min_mapq):
	bcs = []
	for r in bam_in.fetch(chrom, start, end):
	  if r.mapq >= min_mapq:
		  if r.has_tag("BX"):
			  bc_id=r.get_tag("BX")
			  bcs.append(bc_id)
	return tuple(list(set(list(bcs))))

def refine_bcs(**kwargs):

	if 'bed_in' in kwargs:
		bed_in = kwargs['bed_in']
	#if 'shared_in' in kwargs:
	#	shared_in = kwargs['shared_in']
	if 'bam' in kwargs:
		bam_input = kwargs['bam']
	if 'out' in kwargs:
		outpre = kwargs['out']

	bam_open = pysam.Samfile(bam_input)

	bed_df = pd.read_table(bed_in, sep="\t", comment="#", header=None, names=['chrom','start','stop','name','status'])
	
	bed_df['bcs'] = bed_df.apply(lambda row: get_barcode_ids(bam_open,str(row['chrom']),int(row['start']),int(row['stop']),MIN_MAPQ),axis=1)
	print bed_df
	
	bed_grouped = bed_df.groupby('name')
	
	for name, group in bed_grouped:
		print name
		print group
	
	sys.exit()
	
	for name, group in bed_grouped:
	
		group_in = bed_grouped.loc[bed_grouped['status']=="in"]
		if len(group_in.index)>0:
			bc_list_in = group_in['bcs'].tolist()
			if len(bc_list_in)==1:
				common_bcs = bc_list_in[0]
			elif len(bc_list_in)>1:
				common_bcs = set(bc_list_in[0])
				for s in bc_list_in[1:]:
					common_bcs.intersection_update(s)

		group_out = bed_grouped.loc[bed_grouped['status']=="out"]
		if len(group_out.index)==0:
			f = open(outpre + ".bc_list_in.txt", 'w')
			for b in common_bcs:
				f.write(str(b)+"\n")
		
		else:
			bc_list_out=group_out['bcs'].tolist()
			bc_list_out = []	
			bc_list_out_flat = [item for sublist in bc_list_out for item in sublist]
			bc_list_out_flat_uq = list(set(bc_list_out_flat))
			
			bc_final = [x for x in common_bcs if x not in bc_list_out_flat_uq]
			
			f = open(outpre, 'w')
			for b in bc_final:
				f.write(str(b)+"\n")
	
		bam_open.close()
		f.close()
			
		return bc_final
