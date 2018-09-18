import os, sys, argparse, __main__ as main
import pandas as pd, numpy as np
import pysam, vcf
from ast import literal_eval

pd.options.mode.chained_assignment = None

#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
# Parse SV input file to desired format                                       #
#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

### DEFINE FUNCTIONS TO PARSE VCF FILES

def vcf_info_norm(n,bn,c,s,e,olist):
	for record in vcf_reader_norm.fetch(c, s, e):
		print record
		if record.is_snp:
			format_list = (record.FORMAT).split(":")
			geno_list = set(record.genotype(norm_smpl)['GT'].split('|'))
			if 'PS' in format_list and len(geno_list)>1:
				fields = "n,bn,record.genotype(norm_smpl)['PS'], record.CHROM, record.POS, record.REF, record.ALT, record.genotype(norm_smpl)['GT']"
				olist.append(eval(fields))

def vcf_info_tum(n,bn,c,s,e,olist):
	for record in vcf_reader_tum.fetch(c, s, e):
		if record.is_snp:
			format_list = (record.FORMAT).split(":")
			geno_list = set(record.genotype(tum_smpl)['GT'].split('|'))
			if 'PS' in format_list and 'BX' in format_list and len(geno_list)>1:
				bc_1 = record.genotype(tum_smpl)['BX'][0]
				bc_2 = record.genotype(tum_smpl)['BX'][1]
				fields = "n, bn,record.genotype(tum_smpl)['PS'], record.CHROM, record.POS, record.REF, record.ALT, record.genotype(tum_smpl)['GT'], bc_1, bc_2"
				olist.append(eval(fields))     
		
def assign_sv_haps(outpre='out',**kwargs):

	if 'sv' in kwargs:
		sv_input = kwargs['sv']
	if 'vcf_control' in kwargs:
		vcf_norm_input = kwargs['vcf_control']
	if 'vcf_test' in kwargs:
		vcf_tum_input = kwargs['vcf_test']
	if 'out' in kwargs:
		outpre = kwargs['out']
		
	#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
	# GET INFO FROM VCF FILES                                                                #
	#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

	## OPEN VCF FILES

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

	## Generate list of columns to loop through
	sv_wndw = df_sv[['name','name1','chrom1_w','start1_w','stop1_w','name2','chrom2_w','start2_w','stop2_w']].values.tolist()

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

	df_tum = pd.DataFrame(vcf_data_tum)
	df_tum.columns = ['name','bp_name','phase_id_tum','chr','pos','ref_tum','alt_tum','gt_tum', 'bc_1', 'bc_2']

	vcf_merge = pd.merge(df_norm, df_tum, on=['name','bp_name','chr','pos'], how="inner")

	vcf_merge['bc_1_phased'] = vcf_merge.apply(lambda row: row['bc_1'] if row['gt_norm']=='0|1' else row['bc_2'], axis=1)
	vcf_merge['bc_2_phased'] = vcf_merge.apply(lambda row: row['bc_1'] if row['gt_norm']=='1|0' else row['bc_2'], axis=1)


	#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
	# Overlap SV-specific barcodes with barcodes of phased SNVs                   #
	#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

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
	    print type(row['bc_overlap_id'])
	    sv_bcs = literal_eval(row['bc_overlap_id'])
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

	    vcf_merge[str(sv_id) + "_hap1_overlap_count"] = overlap_counts_1
	    vcf_merge[str(sv_id) + "_hap2_overlap_count"] = overlap_counts_2
	   
	    vcf_merge[str(sv_id) + "_hap1_overlap_bcs"] = overlap_bcs_1
	    vcf_merge[str(sv_id) + "_hap2_overlap_bcs"] = overlap_bcs_2


	## ADD COUNTS OF UNIQUE BARCODES FOR EACH BREAKPOINT

	vcf_merge.to_csv("testing.txt", sep="\t", index=False) #debug
	
	bp_list = list(set(vcf_merge['bp_name']))
	sv_list = list(set(vcf_merge['name']))

	df_bp_list = []

	for sv in sv_list:
	    print sv
	    df_bp_name = sv + "_bc_counts_bp"
	    hap1_bc_col = sv + "_hap1_overlap_bcs"
	    hap2_bc_col = sv + "_hap2_overlap_bcs"
	    df_bp_name = vcf_merge.groupby(['name','bp_name','phase_id_norm']).agg({hap1_bc_col:'sum', hap2_bc_col: 'sum'}).reset_index()
	    df_bp_name[sv + "_hap1_overlap_count_bp"] = df_bp_name[hap1_bc_col].apply(lambda x: len(set(x)))
	    df_bp_name[sv + "_hap2_overlap_count_bp"] = df_bp_name[hap2_bc_col].apply(lambda x: len(set(x)))
	    df_bp_name.rename(columns = {hap1_bc_col: sv + "_hap1_overlap_bcs_bp", hap2_bc_col: sv + "_hap2_overlap_bcs_bp"}, inplace=True)
	    df_bp_list.append(df_bp_name)


	df_bp_counts = reduce(lambda x, y: pd.merge(x, y, on = ['name','bp_name','phase_id_norm']), df_bp_list)
	df_bp_counts.to_csv("testing_2.txt", sep="\t", index=False) #debug
	
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
	summ_df = summ_df[['name', 'bp_name', 'phase_id_norm', 'hap1_overlap_bcs_bp', 'hap2_overlap_bcs_bp', 'hap1_overlap_count_bp', 'hap2_overlap_count_bp']]

	summ_df.to_csv(outpre, sep="\t", index=False)
	return summ_df
