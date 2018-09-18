import sys
import os
import __main__ as main
import argparse
import ast
import pandas as pd
import pysam
import numpy as np
import vcf

def barcodeSplit(bc_list):
        new_list = [value.partition('_')[0] for value in bc_list if '_' in value]
        return new_list

def naDrop(bc_list):
        if "n/a" in bc_list:
                new_list = [value for value in bc_list if value!="n/a"]
                return new_list
        else:
                return bc_list

def get_phase_blocks(phased_basic='None',phased_basic_file='None',outpre='out',**kwargs):

### SETTING INPUT ###
	
	if 'infile_basic' in kwargs:
		phased_basic_file = kwargs['infile_basic']
	if 'in_basic' in kwargs:
		phased_basic = kwargs['in_basic']
	if 'out' in kwargs:
		outpre = kwargs['out']

	if (str(phased_basic_file)=='None' and str(phased_basic)=='None'):
		print "Must specify input"
		sys.exit()
	elif (str(phased_basic_file)!='None' and str(phased_basic)!='None'):
		print "Either specify df or file"
		sys.exit()
	
	if str(phased_basic_file)!='None':
		df=pd.read_table(phased_basic_file, sep="\t", header=None, names=['chr','pos','hap1_bc_ct', 'hap1_bcs','hap2_bc_ct','hap2_bcs','gt','PS'], dtype={'chr':str,'gt':str,'PS':str})
	elif str(phased_basic)!='None':
		df = phased_basic
		df.columns=['chr','pos','hap1_bc_ct', 'hap1_bcs','hap2_bc_ct','hap2_bcs','gt','PS']
		df[['chr','gt','PS']] = df[['chr','gt','PS']].astype(str)
	else:
		print "Must specify input"
		sys.exit()	
	
### END OF INPUT ###
	
	df.sort_values(by='PS', inplace=True)
	grouped=df.groupby('PS')

	phase_data=[]

	for name, group in grouped:
		chr=group['chr'].unique()[0]
		beg_pos=group['pos'].min()
		end_pos=group['pos'].max()
		dist=end_pos-beg_pos+1
		all_SNVs=len(group)
	
		group_phsd=group.loc[group['gt'].isin(['0|1','1|0'])]

		phased_het=len(group_phsd)
	
		group_phsd['hap1_bcs_ls']=group_phsd['hap1_bcs'].apply(lambda x: str(x).split(";"))	
		group_phsd['hap1_bcs_rm']=group_phsd['hap1_bcs_ls'].apply(lambda x: barcodeSplit(x))
		group_phsd['hap1_bcs_na']=group_phsd['hap1_bcs_rm'].apply(lambda x: naDrop(x))
		hap1_bcs_all=group_phsd['hap1_bcs_na'].tolist()
		hap1_bcs_all_zip=sum(hap1_bcs_all,[])
		hap1_total=len(hap1_bcs_all_zip)	
		hap1_unique=len(set(hap1_bcs_all_zip))

		group_phsd['hap2_bcs_ls']=group_phsd['hap2_bcs'].apply(lambda x: str(x).split(";"))
		group_phsd['hap2_bcs_rm']=group_phsd['hap2_bcs_ls'].apply(lambda x: barcodeSplit(x))
		group_phsd['hap2_bcs_na']=group_phsd['hap2_bcs_rm'].apply(lambda x: naDrop(x))
		hap2_bcs_all=group_phsd['hap2_bcs_na'].tolist()
		hap2_bcs_all_zip=sum(hap2_bcs_all,[])
		hap2_total=len(hap2_bcs_all_zip)
		hap2_unique=len(set(hap2_bcs_all_zip))	

		bcs_all=hap1_bcs_all_zip+hap2_bcs_all_zip
		total=len(bcs_all)
		unique=len(set(bcs_all))	

		phase_data.append([chr, beg_pos, end_pos, dist, name, all_SNVs, phased_het, total, unique, hap1_total, hap1_unique, hap2_total, hap2_unique])

	df=pd.DataFrame(phase_data)
	df.columns=["chr", "beg_pos", "end_pos", "dist", "PS", "all_SNVs", "phased_het", "total", "unique", "hap1_total", "hap1_unique", "hap2_total", "hap2_unique"]
	#df_sub=df[df['phased_het']>1]

	#df_sub.to_csv(outfile, sep="\t", index=False)
	df.to_csv(str(outpre), sep="\t", index=False)
	return df


