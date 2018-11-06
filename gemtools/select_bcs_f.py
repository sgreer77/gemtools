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
	return list(set(list(bcs)))

def select_bcs(outpre='out',region_in_subset='None',region_out_subset='None',**kwargs):

	if 'region_in' in kwargs:
		region_in_subset = kwargs['region_in']
	if 'region_out' in kwargs:
		region_out_subset = kwargs['region_out']
	if 'bam' in kwargs:
		bam_input = kwargs['bam']
	if 'out' in kwargs:
		outpre = kwargs['out']

	bam_open = pysam.Samfile(bam_input)


	## Get barcodes common to all 'in' regions:
	region_list_in = region_in_subset.split(';')
	bc_list_in = []

	for reg_in in region_list_in:
		chr = str(reg_in).split(',')[0]
		start = int(str(reg_in).split(',')[1])
		stop = int(str(reg_in).split(',')[2])
		region_bc_list_in = get_barcode_ids(bam_open, chr, start, stop, MIN_MAPQ)
		bc_list_in.append(region_bc_list_in)

	if len(bc_list_in)==1:
		common_bcs = bc_list_in[0]
	elif len(bc_list_in)>1:
		common_bcs = set(bc_list_in[0])
		for s in bc_list_in[1:]:
			common_bcs.intersection_update(s)

	if str(region_out_subset)=='None':
		return common_bcs
		f = open(outpre + ".bc_list_in.txt", 'w')
		for b in common_bcs:
			f.write(str(b)+"\n")
	
		bam_open.close()
		f.close()
		sys.exit()
		
	else:
		
		## Get barcodes in out regions:
		region_list_out = region_out_subset.split(';')
		bc_list_out = []

		for reg_out in region_list_out:
			chr = str(reg_out).split(',')[0]
			start = int(str(reg_out).split(',')[1])+200
			stop = int(str(reg_out).split(',')[2])-200
			region_bc_list_out = get_barcode_ids(bam_open, chr, start, stop, MIN_MAPQ)
			bc_list_out.append(region_bc_list_out)		
		
		bc_list_out_flat = [item for sublist in bc_list_out for item in sublist]
		bc_list_out_flat_uq = list(set(bc_list_out_flat))
			
		bc_final = [x for x in common_bcs if x not in bc_list_out_flat_uq]
			
		f = open(outpre, 'w')
		for b in bc_final:
			f.write(str(b)+"\n")
	
		bam_open.close()
		f.close()
			
		return bc_final
