import sys
import os
import __main__ as main
import argparse
import ast
import pandas as pd
import pysam
import numpy as np




def get_barcode_ids(bam_in, chrom, start, end, min_mapq):
	bcs = []
	for r in bam_in.fetch(chrom, start, end):
	  if r.mapq >= min_mapq:
		  if r.has_tag("BX"):
			  bc_id=r.get_tag("BX")
			  bcs.append(bc_id)
	return tuple(list(set(list(bcs))))

def get_shared_bcs(**kwargs):

	if 'bed_in' in kwargs:
		bed_in = kwargs['bed_in']
	if 'bam' in kwargs:
		bam_input = kwargs['bam']
	if 'out' in kwargs:
		outpre = kwargs['out']
	if 'map_qual' in kwargs:
		map_qual = kwargs['map_qual']

	MIN_MAPQ = map_qual

	bam_open = pysam.Samfile(bam_input)

	bed_df = pd.read_table(bed_in, sep="\t", comment="#", header=None, names=['chrom','start','stop','name','sub_name','status'])
	
	bed_df['bcs'] = bed_df.apply(lambda row: get_barcode_ids(bam_open,str(row['chrom']),int(row['start']),int(row['stop']),MIN_MAPQ),axis=1)
	#print bed_df
	
	bed_grouped = bed_df.groupby('sub_name')

	out_data = []

	for name, group in bed_grouped:
	
		sv_name = list(group['name'])[0]
		print sv_name
	
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
			out_data.append([sv_name,name,tuple(common_bcs)])
		
		else:
			bc_list_out=group_out['bcs'].tolist()
			bc_list_out_flat = [item for sublist in bc_list_out for item in sublist]
			bc_list_out_flat_uq = list(set(bc_list_out_flat))
			bc_list_out_flat_uq_len = len(bc_list_out_flat_uq)
			
			bc_final = [x for x in common_bcs if x not in bc_list_out_flat_uq]
			bc_final_len = len(bc_final)
			
			out_data.append([sv_name,name,tuple(bc_final)])
	
	bam_open.close()
	out_df = pd.DataFrame(out_data, columns = ['name','sub_name','select_bcs'])
	
	out_df = out_df[['name','select_bcs']]
	out_df = out_df.loc[out_df['select_bcs']!="na"]
	out_df2 = out_df.groupby('name')['select_bcs'].sum().reset_index()
	out_df2['num_select_bcs'] = out_df2['select_bcs'].apply(lambda x: len(x))
	out_df2 = out_df2[['name','num_select_bcs','select_bcs']]
	out_df2.columns = ['name','num_bcs','bcs']
	

	out_df2.to_csv(outpre, sep="\t", index=False)

