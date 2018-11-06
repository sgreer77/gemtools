import os, sys, argparse, __main__ as main
import pandas as pd, numpy as np
import pysam, vcf

## DEFINE FUNCTION TO OBTAIN BARCODES FROM BAM FILE FOR SPECIFIC REGIONS

def get_barcodes(bam_in, chrom, start, end, min_mapq):
        bcs = set()
        for r in bam_in.fetch(chrom, start, end):
                if r.has_tag("BX") and r.mapq >= min_mapq:
                    bc_id=r.get_tag("BX")
                    bcs.add(bc_id)
        return list(bcs)


#def get_shared_bcs(sv_input, bam_input, outpre):
def get_shared_bcs(outpre='out',**kwargs):
	
	if 'sv' in kwargs:
		sv_input = kwargs['sv']
	if 'bam' in kwargs:
		bam_input = kwargs['bam']
	if 'out' in kwargs:
		outpre = kwargs['out']

	df_sv = pd.read_table(sv_input, sep="\t")
	#df_sv = df_sv.rename(columns = {'#chrom1':'chrom1'})

	## Generate list of columns to loop through
	sv_wndw = df_sv[['name','name1','chrom1_w','start1_w','stop1_w','name2','chrom2_w','start2_w','stop2_w']].values.tolist()


	#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
	# Get SV-specific barcodes from bam file                                      #
	#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#


	## OBTAIN SV-SPECIFIC BARCODES + THEN INTERSECT BETWEEN EVENTS TO LOOK FOR SHARED BARCODES

	MIN_MAPQ = 0

	bam_file = pysam.AlignmentFile(bam_input, "rb")

	bam_data = []

	for (name,name_1,chrom_1,start_1,end_1,name_2,chrom_2,start_2,end_2) in sv_wndw:

	    name_1, chrom_1, name_2, chrom_2 = str(name_1), str(chrom_1), str(name_2), str(chrom_2)
	    start_1, end_1, start_2, end_2 = int(start_1), int(end_1), int(start_2), int(end_2)

	    # Obtain + count SV-specific barcodes
	    bc_1 = get_barcodes(bam_file,chrom_1,start_1,end_1,MIN_MAPQ)
	    bc_1_unq = set(bc_1)
	    bc_1_num = len(bc_1_unq)

	    bc_2 = get_barcodes(bam_file,chrom_2,start_2,end_2,MIN_MAPQ)
	    bc_2_unq = set(bc_2)
	    bc_2_num = len(bc_2_unq)

	    # Intersect SV-specific barcodes
	    bc_overlap = bc_1_unq & bc_2_unq
	    num_overlaps = len(bc_overlap)

	    bam_bc_list = [name, bc_1_num, bc_2_num, num_overlaps, tuple(bc_1_unq), tuple(bc_2_unq), tuple(bc_overlap)]

	    bam_data.append(bam_bc_list)

	# Convert list to data frame + write specified columns to output
	df_bam = pd.DataFrame(bam_data)
	df_bam.columns = ['name','bc_1_num','bc_2_num','bc_overlap_num','bc_1_id','bc_2_id','bc_overlap_id']

	## CREATE TABLE CONTAINING EVENTS ALONG WITH THEIR SV-SPECIFIC BARCODES

	df_sv[['name']], df_bam[['name']] = df_sv[['name']].astype(str), df_bam[['name']].astype(str)
	df_bc = pd.merge(df_sv, df_bam, on='name')


	#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#
	# Pairwise intersection of SV-specific barcodes                               #
	#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#+#

	## PERFORM PAIRWISE INTERSECTION OF SV CANDIDATES (i.e. INTERSECTION OF INTERSECTION)
	sv_inter = []

	sv_spec_bcs = df_bc['bc_overlap_id'].values.tolist()
	sv_names = list(df_bc['name'])

	for bcs in sv_spec_bcs:
	    sv_inter.append([set(bcs) & set(x) for x in sv_spec_bcs]) # '&' performs the intersection


	## COUNT BARCODES IN LISTS FROM PREVIOUS STEP + GENERATE ADJACENCY MATRIX

	sv_inter_matrix = []

	for inter in sv_inter:
	    sv_inter_matrix.append([len(x) for x in inter])
	sv_inter_matrix = np.matrix(sv_inter_matrix)

	## CONVERT TO DF + WRITE TO OUTPUT

	df_inter = pd.DataFrame(sv_inter_matrix, columns=sv_names)
	#print df_inter
	#print sv_names
	df_inter['name'] = sv_names

	df_isect = pd.merge(df_bc, df_inter, on='name')
	#print df_merge

	df_isect.to_csv(outpre, sep="\t", index=False)
	
	return df_isect
