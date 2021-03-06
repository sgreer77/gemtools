import sys
import os
import __main__ as main
import ast
import pandas as pd
import pysam
import numpy as np
import vcf

def barcodeSplit(bc_str):
		bc_list = str(bc_str).split(";")
		bc_list_pre = [value.partition('_')[0] for value in bc_list if '_' in value]
		if "n/a" in bc_list_pre:
			nona_list = [value for value in bc_list_pre if value!="n/a"]
			return nona_list
		else:
			return bc_list_pre

def get_phased_bcs(phased_basic_file='None',outpre='out',phase_block='None',**kwargs):

### SETTING INPUT ###

	if 'infile_basic' in kwargs:
		phased_basic_file = kwargs['infile_basic']
	if 'out' in kwargs:
		outpre = kwargs['out']
	if 'ps' in kwargs:
		phase_block = kwargs['ps']
	
	if str(phased_basic_file)!='None':
		df=pd.read_table(phased_basic_file, sep="\t", dtype={'#chrom':str,'gt':str,'block_id':str,'bc1':str,'bc2':str})

	else:
		print "Must specify input"
		sys.exit()	

### END OF INPUT ###

	df.sort_values(by='block_id', inplace=True)
	grouped=df.groupby('block_id')

	phase_data=[]

	for name, group in grouped:
		
		if str(name)!="n/a":
			name = str(int(float(name)))
			
			if str(name)==str(phase_block):
				chr=group['#chrom'].unique()[0]
				beg_pos=group['pos'].min()
				end_pos=group['pos'].max()
				dist=end_pos-beg_pos+1
				all_SNVs=len(group)
	
				#group_phsd=group.loc[group['gt'].isin(['0|1','1|0'])]
				group_phsd = group.loc[(group['phase_status']=="phased") & (group['hom_status']=="het") & (group['var_type']=="snv") & (group['filter']=="[]")]
				phased_het=len(group_phsd)
	
				group_phsd['bc1_ls']=group_phsd['bc1'].apply(lambda x: barcodeSplit(x))
				bc1_all=group_phsd['bc1_ls'].tolist()
				bc1_all_zip=sum(bc1_all,[])
				bc1_total=len(bc1_all_zip)	
				bc1_unique=len(set(bc1_all_zip))

				group_phsd['bc2_ls']=group_phsd['bc2'].apply(lambda x: barcodeSplit(x))
				bc2_all=group_phsd['bc2_ls'].tolist()
				bc2_all_zip=sum(bc2_all,[])
				bc2_total=len(bc2_all_zip)
				bc2_unique=len(set(bc2_all_zip))	
	
				bcs_all=bc1_all_zip+bc2_all_zip
				total=len(bcs_all)
				unique=len(set(bcs_all))

				phase_data.append([chr, beg_pos, end_pos, dist, name, all_SNVs, phased_het, total, unique, bc1_total, bc1_unique, bc2_total, bc2_unique, tuple(list(set(bc1_all_zip))), tuple(list(set(bc2_all_zip)))])

	df=pd.DataFrame(phase_data)
	df.columns=["chr", "beg_pos", "end_pos", "dist", "PS", "all_SNVs", "phased_het", "total", "unique", "hap1_total", "hap1_unique", "hap2_total", "hap2_unique","hap1_uq_bcs","hap2_uq_bcs"]

	df.to_csv(str(outpre), sep="\t", index=False)
	return df


