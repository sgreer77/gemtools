import sys
import os
import __main__ as main
import ast
import pandas as pd
import pysam
import numpy as np


MIN_MAPQ = 0

## DEFINE FUNCTION TO OBTAIN BARCODES FROM BAM FILE FOR SPECIFIC REGIONS

def get_barcode_ids(bam_in, chrom, start, end, min_mapq):
	bcs = []
	for r in bam_in.fetch(chrom, start, end):
	  if r.mapq >= min_mapq:
		  if r.has_tag("BX"):
			  bc_id=r.get_tag("BX")
			  bcs.append(bc_id)
	return list(bcs)


def get_bcs_in_region(outpre='out',region_subset='None',**kwargs):

	if 'region' in kwargs:
		region_subset = kwargs['region']
	if 'bam' in kwargs:
		bam_input = kwargs['bam']
	if 'out' in kwargs:
		outpre = kwargs['out']

	if str(region_subset)=="None":
		sys.exit()
		
	bam_open = pysam.Samfile(bam_input)

	# Create data frame of 1kb windows
	
	region_list = region_subset.split(';')
	bc_list = []
	
	for reg in region_list:
		chr = str(reg).split(',')[0]
		start = int(str(reg).split(',')[1])
		stop = int(str(reg).split(',')[2])
		region_bc_list = get_barcode_ids(bam_open, chr, start, stop, MIN_MAPQ, PERF_CIGAR)
		bc_list.append(region_bc_list)
	
	flat_bc_list = [item for sublist in bc_list for item in sublist]
	flat_bc_list_uq = list(set(flat_bc_list))

	f = open(outpre, 'w')
	for b in flat_bc_list_uq:
		f.write(str(b)+"\n")
	
	bam_open.close()
	f.close()
	return flat_bc_list_uq
