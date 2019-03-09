import os
import sys
import argparse
import __main__ as main
import pandas as pd
import numpy as np
import pysam
import vcf
from ast import literal_eval

pd.options.mode.chained_assignment = None

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Parse SV input file to desired format                                       #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

### DEFINE FUNCTIONS TO PARSE VCF FILES

def vcf_info_norm(n,bn,c,s,e,olist):
	for record in vcf_reader_norm.fetch(c, s, e):
		if record.is_snp:
			format_list = (record.FORMAT).split(":")
			geno_list = set(record.genotype(norm_smpl)['GT'].split('|'))
			if 'PS' in format_list and len(geno_list)>1:
				fields = "n,bn,record.genotype(norm_smpl)['PS'], record.CHROM, record.POS, record.REF, tuple(record.ALT), record.genotype(norm_smpl)['GT']"
				olist.append(eval(fields))

def vcf_info_tum(n,bn,c,s,e,olist):
	for record in vcf_reader_tum.fetch(c, s, e):
		if record.is_snp:
			format_list = (record.FORMAT).split(":")
			geno_list = set(record.genotype(tum_smpl)['GT'].split('|'))
			if 'PS' in format_list and 'BX' in format_list and len(geno_list)>1:
				bc_1 = record.genotype(tum_smpl)['BX'][0]
				bc_2 = record.genotype(tum_smpl)['BX'][1]
				fields = "n, bn,record.genotype(tum_smpl)['PS'], record.CHROM, record.POS, record.REF, tuple(record.ALT), record.genotype(tum_smpl)['GT'], bc_1, bc_2"
				olist.append(eval(fields))     

def make_window_haps(s,e,w):
    cur_size = int(e)-int(s)
    adj_val = (int(w)-cur_size)/2
    adj_val = int(round(adj_val,0))
    new_start = max(0,s - adj_val)
    new_end = e + adj_val
    return [new_start,new_end]

def window_rows_haps(r):
    wndw_row = [ r['name'],r['chrom1'],r['start1'],r['stop1'],r['chrom2'],r['start2'],r['stop2'],r['name1'],r['chrom1']] + make_window_haps(r['start1'],r['stop1'],r['window_size']) + [r['name2'],r['chrom2']] + make_window_haps(r['start2'],r['stop2'],r['window_size']) + [r['window_size']]
    return wndw_row

	
def assign_sv_haps(**kwargs):

	if 'sv' in kwargs:
		sv_input = kwargs['sv']
	if 'vcf_control' in kwargs:
		vcf_norm_input = kwargs['vcf_control']
	if 'vcf_test' in kwargs:
		vcf_tum_input = kwargs['vcf_test']
	if 'out' in kwargs:
		outpre = kwargs['out']
	if 'window' in kwargs:
		window_size = kwargs['window']
	if 'bcs' in kwargs:
		bc_subset = kwargs['bcs']
		
	#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
	# GET INFO FROM VCF FILES                                                                #
	#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

	## OPEN VCF FILES

	if str(vcf_norm_input)=="None":
		vcf_norm_input = vcf_tum_input

	global vcf_reader_norm 
	vcf_reader_norm = vcf.Reader(filename=vcf_norm_input)
	global norm_smpl 
	norm_smpl = vcf_reader_norm.samples[0]
	vcf_data_norm = []

	global vcf_reader_tum
	vcf_reader_tum = vcf.Reader(filename=vcf_tum_input)
	global tum_smpl 
	tum_smpl = vcf_reader_tum.samples[0]
	vcf_data_tum = []


	## OBTAIN SNVs + INFO FROM NORMAL VCF FILE -- will use normal file to define phase blocks + phase of variants 

	df_sv = pd.read_table(sv_input, sep="\t")
	#df_sv = df_sv.rename(columns = {'#chrom1':'chrom1'})

	if str(window_size)=='None':
	## Generate list of columns to loop through
		sv_wndw = df_sv[['name','name1','chrom1_w','start1_w','stop1_w','name2','chrom2_w','start2_w','stop2_w']].values.tolist()
	else:
		df_sv_sub = df_sv[['name','name1','chrom1','start1','stop1','name2','chrom2','start2','stop2']]
		df_sv_sub['window_size'] = window_size
		wndw_out = df_sv_sub.apply(lambda row: window_rows_haps(row), axis=1)
		df_wndw = pd.DataFrame(list(wndw_out))
		df_wndw.columns = ['name','chrom1','start1','stop1','chrom2','start2','stop2','name1','chrom1_w','start1_w','stop1_w','name2','chrom2_w','start2_w','stop2_w','window_size']
		sv_wndw = df_wndw[['name','name1','chrom1_w','start1_w','stop1_w','name2','chrom2_w','start2_w','stop2_w']].values.tolist()
		#df_wndw.to_csv("hap_wndws.txt", sep="\t", index=False)

		if bc_subset=="shared":
			df_sv_sub2 = df_sv[['name','bc_overlap_num','bc_overlap_id']]
			#df_sv_sub2 = df_sv.loc[df_sv['bc_overlap_id']!="na"]
		elif bc_subset=="select":
			df_sv_sub2 = df_sv[['name','num_select_bcs','select_bcs']]
			#df_sv_sub2 = df_sv.loc[df_sv['select_bcs']!="na"]
		else:
			print "bcs -- must be either 'shared' or 'select'"
			sys.exit(1)
		
		df_sv_sub2[['name']]=df_sv_sub2[['name']].astype(str)
		df_wndw[['name']]=df_wndw[['name']].astype(str)

		df_sv = pd.merge(df_wndw, df_sv_sub2, on="name")

		#df_sv.to_csv("testing_df_sv.txt", sep="\t", index=False)

	for (name,name_1,chrom_1,start_1,end_1,name_2,chrom_2,start_2,end_2) in sv_wndw:
		name,name_1,chrom_1,name_2,chrom_2 = str(name),str(name_1),str(chrom_1),str(name_2),str(chrom_2)
		start_1,end_1,start_2,end_2 = int(start_1),int(end_1),int(start_2),int(end_2)

		vcf_info_norm(name,name_1,chrom_1,start_1,end_1,vcf_data_norm)
		vcf_info_norm(name,name_2,chrom_2,start_2,end_2,vcf_data_norm)

	## OBTAIN SNVs + INFO FROM TUMOR VCF FILE -- will use tumor file to obtain barcodes of phased variants

		vcf_info_tum(name,name_1,chrom_1,start_1,end_1,vcf_data_tum)
		vcf_info_tum(name,name_2,chrom_2,start_2,end_2,vcf_data_tum)


	## MERGE VCF INFO DATA FRAMES THEN ADJUST BARCODE COLUMN ORDER TO REFLECT HAPLOTYPES

	df_norm = pd.DataFrame(vcf_data_norm)
	df_norm.columns = ['name','bp_name','phase_id_norm','chr','pos','ref_norm','alt_norm','gt_norm']
	print df_norm.head()
	df_norm.drop_duplicates(inplace=True)

	df_tum = pd.DataFrame(vcf_data_tum)
	df_tum.columns = ['name','bp_name','phase_id_tum','chr','pos','ref_tum','alt_tum','gt_tum', 'bc_1', 'bc_2']
	print df_tum.head()
	df_tum.drop_duplicates(inplace=True)

	vcf_merge = pd.merge(df_norm, df_tum, on=['name','bp_name','chr','pos'], how="inner")
	print vcf_merge.head()
	vcf_merge.drop_duplicates(inplace=True)

	vcf_merge['bc_1_phased'] = vcf_merge.apply(lambda row: tuple(row['bc_1'].split(";")) if row['gt_norm']=='0|1' else tuple(row['bc_2'].split(";")), axis=1)
	vcf_merge['bc_2_phased'] = vcf_merge.apply(lambda row: tuple(row['bc_1'].split(";")) if row['gt_norm']=='1|0' else tuple(row['bc_2'].split(";")), axis=1)

	vcf_merge.to_csv("test1.txt", sep="\t", index=False)

	#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
	# Overlap SV-specific barcodes with barcodes of phased SNVs                   #
	#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#


	#df.groupby('A').agg({'B': ['min', 'max'], 'C': 'sum'})
	vcf_merge_grp = vcf_merge.groupby(['name','bp_name','phase_id_norm']).agg({'bc_1_phased' : 'sum', 'bc_2_phased' : 'sum'}).reset_index()

	vcf_merge_grp.to_csv("test2.txt", sep="\t", index=False)

	## CREATE A LIST OF HAPLOTYPE-SPECIFIC BARCODES FOR EACH PHASED VARIANT 

	bc_1_list = vcf_merge['bc_1_phased'].tolist()
	bc_1_list_mod = []
	for bc_ls in bc_1_list:
	    new_bc_ls = bc_ls.split(';')
	    new_bc_ls_mod = [b.split('_')[0] for b in new_bc_ls]
	    bc_1_list_mod.append(new_bc_ls_mod)

	bc_2_list = vcf_merge['bc_2_phased'].tolist()
	bc_2_list_mod = []
	for bc_ls in bc_2_list:
	    new_bc_ls = bc_ls.split(';')
	    new_bc_ls_mod = [b.split('_')[0] for b in new_bc_ls]
	    bc_2_list_mod.append(new_bc_ls_mod)

	## INTERSECT "SV-SPECIFIC" BARCODES WITH LIST OF BCS FOR EACH PHASED VARIANT (GENERATED ABOVE) -- COUNT NUMBER OF OVERLAPPING BARCODES

	for index, row in df_sv.iterrows():    

	    #print row['bc_overlap_id']
	    #print type(row['bc_overlap_id'])
	    
	    if bc_subset=="shared":
	    	if str(row['bc_overlap_id'])!="None":
	    		sv_bcs = literal_eval(row['bc_overlap_id'])
	    elif bc_subset=="select":
	    	if int(row['num_select_bcs'])!=0:
	    		sv_bcs = literal_eval(row['select_bcs'])
	    	else:
	    		print "No shared barcodes for: " + str(row["name"])
	    else:
	    	print "bcs -- must be either 'shared' or 'select'"
	    	sys.exit(1)
	    
	    num_bcs = len(sv_bcs) ## new_add
	    
	    
	    #sv_bcs = row['bc_overlap_id']
	    sv_id = row['name']

	    overlap_bcs_1 = []
	    overlap_counts_1 = []
	    for b in bc_1_list_mod:
			bc_overlap_1 = list(set(sv_bcs) & set(b))
			num_overlaps = len(bc_overlap_1)
			overlap_bcs_1.append(tuple(bc_overlap_1))
			overlap_counts_1.append(num_overlaps)

	    overlap_bcs_2 = []
	    overlap_counts_2 = []
	    for b in bc_2_list_mod:
			bc_overlap_2 = list(set(sv_bcs) & set(b))
			num_overlaps = len(bc_overlap_2)
			overlap_bcs_2.append(tuple(bc_overlap_2))
			overlap_counts_2.append(num_overlaps)

	    vcf_merge[str(sv_id) + "_num_bcs_checked"] = num_bcs
		
	    vcf_merge[str(sv_id) + "_hap1_overlap_count"] = overlap_counts_1
	    vcf_merge[str(sv_id) + "_hap2_overlap_count"] = overlap_counts_2
	   
	    vcf_merge[str(sv_id) + "_hap1_overlap_bcs"] = overlap_bcs_1
	    vcf_merge[str(sv_id) + "_hap2_overlap_bcs"] = overlap_bcs_2
		

	## ADD COUNTS OF UNIQUE BARCODES FOR EACH BREAKPOINT

	#vcf_merge.to_csv("vcf_merge_testing.txt", sep="\t", index=False) #debug
	
	bp_list = list(set(vcf_merge['bp_name']))
	sv_list = list(set(vcf_merge['name']))

	df_bp_list = []

	for sv in sv_list:
	    df_bp_name = sv + "_bc_counts_bp"
	    hap1_bc_col = sv + "_hap1_overlap_bcs"
	    hap2_bc_col = sv + "_hap2_overlap_bcs"
	    sv_num_bcs = sv + "_num_bcs_checked"
	    df_bp_name = vcf_merge.groupby(['name','bp_name','phase_id_norm',sv_num_bcs]).agg({hap1_bc_col:'sum', hap2_bc_col: 'sum'}).reset_index()
	    df_bp_name[sv + "_hap1_overlap_count_bp"] = df_bp_name[hap1_bc_col].apply(lambda x: len(set(x)))
	    df_bp_name[sv + "_hap2_overlap_count_bp"] = df_bp_name[hap2_bc_col].apply(lambda x: len(set(x)))
	    df_bp_name.rename(columns = {hap1_bc_col: sv + "_hap1_overlap_bcs_bp", hap2_bc_col: sv + "_hap2_overlap_bcs_bp", sv_num_bcs: sv + "_num_bcs_checked"}, inplace=True)
	    df_bp_list.append(df_bp_name)
	    #df_bp_name.to_csv(sv + "_bp_list.txt", sep="\t", index=False)

	df_bp_counts = reduce(lambda x, y: pd.merge(x, y, on = ['name','bp_name','phase_id_norm']), df_bp_list)
	#df_bp_counts.to_csv("testing_counts.txt", sep="\t", index=False) #debug
	
	## CREATE FINAL SUMMARY OUTPUT

	sv_list = list(set(df_bp_counts['name']))

	general_cols = ['name','bp_name','phase_id_norm']
	all_cols = df_bp_counts.columns

	df_sub_list = []

	for sv in sv_list:
	    df_sub_name = sv + "_bc_counts_sub"
	    df_sub = df_bp_counts.loc[df_bp_counts['name']==sv]
	    sv_cols = [c for c in all_cols if c.startswith(str(sv) + "_")]
	    sub_cols = general_cols + sv_cols
	    df_sub = df_sub[sub_cols]
	    #df_sub_cols = general_cols + ['_'.join(n.split("_")[2:]) for n in sv_cols]
	    df_sub_cols = general_cols + [str(n).partition(str(sv) + "_")[2] for n in sv_cols]
	    df_sub.columns = df_sub_cols

	    df_sub_list.append(df_sub)

	summ_df = pd.concat(df_sub_list)
	#summ_df = summ_df[['name', 'bp_name', 'phase_id_norm', 'hap1_overlap_bcs_bp', 'hap2_overlap_bcs_bp', 'hap1_overlap_count_bp', 'hap2_overlap_count_bp', 'hap1_overlap_bcs_sv', 'hap2_overlap_bcs_sv', 'hap1_overlap_count_sv', 'hap2_overlap_count_sv','both_overlap_count_bp']]
	#summ_df.to_csv("testing_2.txt", sep="\t", index=False) #debug
	summ_df = summ_df[['name', 'bp_name', 'phase_id_norm', 'num_bcs_checked', 'hap1_overlap_bcs_bp', 'hap2_overlap_bcs_bp', 'hap1_overlap_count_bp', 'hap2_overlap_count_bp']]

	summ_df.to_csv(outpre, sep="\t", index=False)
	return summ_df
