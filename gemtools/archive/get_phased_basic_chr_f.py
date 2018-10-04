import sys
import os
import __main__ as main
import ast
import pandas as pd
import pysam
import numpy as np
import vcf

def get_phased_basic_chr(inputvcf='None',outpre='out',chrom='None',**kwargs):

	if 'vcf' in kwargs:
		inputvcf = kwargs['vcf']
	if 'out' in kwargs:
		outpre = kwargs['out']
	if 'chr' in kwargs:
		chrom = kwargs['chr']
	
	if str(chrom)=='None':
		sys.exit

	vcf_reader = vcf.Reader(open(inputvcf,'r'))

	sample_list = vcf_reader.samples
	cur_sample = sample_list[0]

	vcf_data=[]
	for record in vcf_reader.fetch(str(chrom)):
		if record.is_snp and record.genotype(cur_sample)['GT'] in ["0|1","1|1","0/1","1/1","1|0","1/0"]:
			format_field=(record.FORMAT).split(":")
			if 'PS' in format_field:
				block_id=record.genotype(cur_sample)['PS']
			else:
				block_id="n/a"
			if 'BX' in format_field:
				bc1=record.genotype(cur_sample)['BX'][0]
				bc2=record.genotype(cur_sample)['BX'][1]
				if bc1=='':
					bc1="n/a"	
					bc1_ct=0
				else:
					bc1_ct=len(bc1.split(";"))
				if bc2=='':
					bc2="n/a"
					bc2_ct=0
				else:
					bc2_ct=len(bc2.split(";"))
				if record.genotype(cur_sample)["GT"]=="1|0":
					bc1,bc2=bc2,bc1
					bc1_ct,bc2_ct=bc2_ct,bc1_ct
			vcf_data.append((record.CHROM, record.POS, bc1_ct, bc1, bc2_ct, bc2, record.genotype(cur_sample)["GT"], block_id))

	df=pd.DataFrame(vcf_data)

	df.to_csv(str(outpre), sep="\t", index=False, header=False)

	return df
