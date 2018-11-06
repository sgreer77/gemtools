#!/usr/bin/env python

import pandas as pd
from pybedtools import BedTool
import pysam

ps_size=50

file_ps_1 = "878.phase_blocks" 
file_ps_2 = "891.phase_blocks"

file_tbx_1 = "878.phased_basic.gz"
file_tbx_2 = "891.phased_basic.gz"

df_ps_1 = pd.read_table(file_ps_1, sep="\t")
df_ps_2 = pd.read_table(file_ps_2, sep="\t")

tbx_1 = pysam.TabixFile(file_tbx_1)
tbx_2 = pysam.TabixFile(file_tbx_2)

## Get headers from the bed files
for h in tbx_1.header:
	tbx_1_header=h.split('\t')
	tbx_1_header=[x.replace('#','') for x in tbx_1_header]

for h in tbx_2.header:
	tbx_2_header=h.split('\t')
	tbx_2_header=[x.replace('#','') for x in tbx_2_header]

##### FUNCTIONS ######

def ps_to_bed(d,s):
	df = d.loc[d['phased_het']>=int(ps_size)]
	df = df[["chr","beg_pos","end_pos","PS"]]
	df.columns = ["#chr","beg_pos","end_pos",str(s) + "_PS"]
	return df

def make_bed_obj(b):
	str_obj = b.to_string(index=False)
	bed_obj=BedTool(str_obj, from_string=True)
	return bed_obj

def get_tbx_info(r,h,i):
	tbx_row_list = r.split('\t')
	tbx_row_dict = dict(zip(h, tbx_row_list))
	
	tbx_chrom = str(tbx_row_dict.get('chrom'))
	tbx_pos = int(tbx_row_dict.get('pos'))
	tbx_id = str(tbx_row_dict.get('block_id'))
	tbx_base1 = str(tbx_row_dict.get('base_1'))
	tbx_base2 = str(tbx_row_dict.get('base_2'))
	tbx_var = str(tbx_row_dict.get('var_type'))
	tbx_phase = str(tbx_row_dict.get('phase_status'))
	tbx_hom = str(tbx_row_dict.get('hom_status'))
	tbx_filt = str(tbx_row_dict.get('filter'))
	
	if (tbx_var=='snv' and tbx_filt=='[]' and tbx_id==str(i) and ((tbx_phase=='phased' and tbx_hom=='het') or (tbx_hom=='hom'))):
		return [tbx_chrom, tbx_pos, tbx_id, tbx_base1, tbx_base2]

def count_alleles(a1,a2,a3,a4):
	allele_list = [a1,a2,a3,a4]
	allele_list = [str(i) for i in allele_list]
	allele_list = list(set(allele_list))
	allele_count = len(allele_list)
	return allele_count
	

def compare_haps(b1,b2):
	num_match = sum(len(set(i))==1 for i in zip(b1,b2))
	return num_match
	

#######################

# BEDTOOLS INTERSECTIONS TO FIND SHARED PHASED REGIONS

# Create bedtools objects in preparation for intersection
df_bed_1 = ps_to_bed(df_ps_1,"1")
df_bed_2 = ps_to_bed(df_ps_2,"2")
bed_1 = make_bed_obj(df_bed_1)
bed_2 = make_bed_obj(df_bed_2)

# First intersection: get the PS values for parent + child
bed_ps = bed_1.intersect(bed_2, wa=True, wb=True)
df_ps = bed_ps.to_dataframe()
df_ps.columns = df_bed_1.columns.tolist() + df_bed_2.columns.tolist()
ps_cols = [c for c in df_ps.columns.tolist() if c.endswith("_PS")]
df_ps = df_ps[ps_cols]

# Second intersection: Get the intersection coordinates
bed_overlap = bed_1.intersect(bed_2)
df_overlap = bed_overlap.to_dataframe()
df_overlap = df_overlap.iloc[:,0:3]
df_overlap.columns = ['#overlap_chr','overlap_start','overlap_end']

# Concatenate the 2 intersection results
df_regions = pd.concat([df_overlap, df_ps], axis=1)
df_regions.sort_values(by=['#overlap_chr',"overlap_start"], inplace=True)

df_regions['dist'] = df_regions['overlap_end']-df_regions['overlap_start']
df_regions.to_csv("overlaps.txt", sep="\t", index=False)

##### COMPARE ALLELES IN EACH SHARED REGION

block_info = []

for index, row in df_regions.iterrows():

	# SUMMARIZE VARIANT INFO
	c = row['#overlap_chr']#.replace("chr","")
	s,e = row['overlap_start'],row['overlap_end']
	ps_1 = str(row['1_PS'])
	ps_2 = str(row['2_PS'])

	shared_block_size = e-s

	## loop through first phased basic
	tbx_list_1 = []
	for tbx_row in tbx_1.fetch(str(c), int(s), int(e)):
		tbx_info = get_tbx_info(tbx_row, tbx_1_header, ps_1)
		if str(tbx_info)!="None":
			tbx_list_1.append(tbx_info)
	df_tbx_1 = pd.DataFrame(tbx_list_1, columns=['chrom','pos','ps_1','base1_i1','base2_i1',])

	## loop through second phased basic
	tbx_list_2 = []
	for tbx_row in tbx_2.fetch(str(c), int(s), int(e)):
		tbx_info = get_tbx_info(tbx_row, tbx_2_header, ps_2)
		if str(tbx_info)!="None":
			tbx_list_2.append(tbx_info)
	df_tbx_2 = pd.DataFrame(tbx_list_2, columns=['chrom','pos','ps_2','base1_i2','base2_i2',])
	
	## merge the two
	df_merge = pd.merge(df_tbx_1, df_tbx_2, on=['chrom','pos'], how='inner')
	
	## get rid of positions with no variant at all
	df_merge['allele_count'] = df_merge.apply(lambda r: count_alleles(r['base1_i1'],r['base2_i1'],r['base1_i2'],r['base2_i2']), axis=1)
	num_pos = len(df_merge.index)
	df_merge = df_merge.loc[df_merge['allele_count']>1]
	df_merge.to_csv(str(c) + "_" + str(s) + "_" + str(e) + ".txt", sep="\t", index=False)
	
	num_var = len(df_merge.index)
	print num_var
	
	if num_var>ps_size:

		## compare each haplotype
		base1_i1_list = df_merge['base1_i1'].tolist()
		base1_i2_list = df_merge['base1_i2'].tolist()
		base2_i1_list = df_merge['base2_i1'].tolist()
		base2_i2_list = df_merge['base2_i2'].tolist()
		
		base1_i1_v_base1_i2 = compare_haps(base1_i1_list,base1_i2_list)
		base1_i1_v_base2_i2 = compare_haps(base1_i1_list,base2_i2_list)
		base2_i1_v_base1_i2 = compare_haps(base2_i1_list,base1_i2_list)
		base2_i1_v_base2_i2 = compare_haps(base2_i1_list,base2_i2_list)
		
		perc_base1_i1_v_base1_i2 = base1_i1_v_base1_i2/float(num_var)
		perc_base1_i1_v_base2_i2 = base1_i1_v_base2_i2/float(num_var)
		perc_base2_i1_v_base1_i2 = base2_i1_v_base1_i2/float(num_var)
		perc_base2_i1_v_base2_i2 = base2_i1_v_base2_i2/float(num_var)
		
		num_gt = sum(i > 0.90 for i in [perc_base1_i1_v_base1_i2,perc_base1_i1_v_base2_i2,perc_base2_i1_v_base1_i2,perc_base2_i1_v_base2_i2])
		
		block_info.append([c,s,e,shared_block_size,ps_1,ps_2,num_pos,num_var,perc_base1_i1_v_base1_i2,perc_base1_i1_v_base2_i2,perc_base2_i1_v_base1_i2,perc_base2_i1_v_base2_i2,num_gt])
	
	
block_df = pd.DataFrame(block_info, columns = ['chrom','start','end','block_size','ps_1','ps_2','num_pos','num_var_pos','hap1_i1_v_hap1_i2','hap1_i1_v_hap2_i2','hap2_i1_v_hap1_i2','hap2_i1_v_hap2_i2','num_gt'])

block_df.to_csv('block_891.txt', sep="\t", index=False)

