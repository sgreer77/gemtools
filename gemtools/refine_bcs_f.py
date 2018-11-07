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
	if 'shrd_file' in kwargs:
		shared_in = kwargs['shrd_file']
	if 'bam' in kwargs:
		bam_input = kwargs['bam']
	if 'out' in kwargs:
		outpre = kwargs['out']

	bam_open = pysam.Samfile(bam_input)

	bed_df = pd.read_table(bed_in, sep="\t", comment="#", header=None, names=['chrom','start','stop','name','status'])
	
	bed_df['bcs'] = bed_df.apply(lambda row: get_barcode_ids(bam_open,str(row['chrom']),int(row['start']),int(row['stop']),MIN_MAPQ),axis=1)
	#print bed_df
	
	bed_grouped = bed_df.groupby('name')

	out_data = []

	for name, group in bed_grouped:
	
		#print name
		group_in = group.loc[group['status']=="in"]

		if len(group_in.index)>0:
			bc_list_in = group_in['bcs'].tolist()

			if len(bc_list_in)==1:
				common_bcs = bc_list_in[0]
			elif len(bc_list_in)>1:
				common_bcs = set(bc_list_in[0])
				for s in bc_list_in[1:]:
					common_bcs.intersection_update(s)
		
		common_bcs = list(common_bcs)
		common_bcs_len = len(common_bcs)

		group_out = group.loc[group['status']=="out"]
		if len(group_out.index)==0:
			out_data.append([name,common_bcs_len,"na","na",tuple(common_bcs),"na","na"])
		
		else:
			bc_list_out=group_out['bcs'].tolist()
			bc_list_out_flat = [item for sublist in bc_list_out for item in sublist]
			bc_list_out_flat_uq = list(set(bc_list_out_flat))
			bc_list_out_flat_uq_len = len(bc_list_out_flat_uq)
			
			bc_final = [x for x in common_bcs if x not in bc_list_out_flat_uq]
			bc_final_len = len(bc_final)
			
			out_data.append([name,common_bcs_len,bc_list_out_flat_uq_len,bc_final_len,tuple(common_bcs),tuple(bc_list_out_flat_uq),tuple(bc_final)])
	
	bam_open.close()
	out_df = pd.DataFrame(out_data, columns = ['name','num_in_bcs','num_out_bcs','num_select_bcs','in_bcs','out_bcs','select_bcs'])
	
	if str(shared_in)=="None":
		out_df.to_csv(outpre, sep="\t", index=False)
	else:
		shared_df = pd.read_table(shared_in, sep="\t")
		merged = pd.merge(shared_df, out_df, on="name", how="left")
		merged.to_csv(outpre, sep="\t", index=False)
