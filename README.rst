gemtools2
---------

Installation
============

**Dependencies:**
"""""""""""""""""
**Python packages:** pandas, numpy, pysam, pyvcf, rpy2

**R package:** ggplot2

**To install gemtools2:**
"""""""""""""""""""""""""
pip install git+https://github.com/sgreer77/gemtools2.git


Running gemtools2
=================

gemtools2 is a collection of tools that use the **output of Long Ranger** (10X Genomics) to perform additional analyses

Long Ranger output files are indicated in the instructions below with an **LR** prefix

**General usage: gemtools2 -T [tool] [-options]**

Tools for getting basic information about the phase blocks:
""""""""""""""""""""""""""""""""""""""""""""""""""""

**get_phased_basic**: Obtain phasing information for all SNVs in the vcf file

	gemtools2 -T get_phased_basic -v [LR.vcf.gz] -o [output.phased_basic]
	
	Ex: gemtools2 -T get_phased_basic -v phased_variants.vcf.gz -o output.phased_basic
	
	Input:
		-v gzipped vcf file output from LR
	Output:
		-o output file: each row is an SNV; columns are: chr,pos,number of hap1 barcodes, hap1 barcodes, number of hap2 barcodes, hap2 barcodes, genotype, phase block id

**get_phase_blocks**: Summarize phase blocks -- coordinates, size, number of phased heterozygous SNVs per phase block etc.

	gemtools2 -T get_phase_blocks -i [output.phased_basic] -o [output.phase_blocks]
	
	Ex: gemtools2 -T get_phase_blocks -i output.phased_basic -o output.phase_blocks
	
	Input:
		-i output from 'get_phased_basic' tool
	Output:
		-o output file: each row is a phase block, columns are phase block coordinates, phase block size, number of phased heterozygotes etc.

Generally useful tools:
""""""""""""""""""""""""""

**get_phased_bcs**

	gemtools2 -T get_phased_bcs -i [output.phased_basic] -p [phase_block_id] -o [out.phased_bcs]

**select_bcs**

	gemtools2 -T select_bcs -b [LR.bam] -f [region_in] -g [region_out] -o [out.bcs]

**get_bcs_in_region**

	gemtools2 -T get_bcs_in_region -b [LR.bam] -f [region_in] -o [out.bcs]

**count_bcs_list**

	gemtools2 -T count_bcs_list -b [LR.bam] -f [region_in] -x [in_window] -b [bc_list.txt] -o [out.bc_count]

**plot_hmw**

	gemtools2 -T plot_hmw -i [out.bc_count] -o [out.pdf]

SV analysis tools:
"""""""""""""""""""""

**bedpe2window**

	gemtools2 -T bedpe2window -i [LR_input.bedpe] -w [window_size] -o [out.bed]

**get_shared_bcs**

	gemtools2 -T get_shared_bcs -i [out.bed] -b [LR_bam_file] -o [out.shared]

**assign_sv_haps**

	gemtools2 -T assign_sv_haps -i [out.shared] -c [LR_control.vcf.gz] -t [LR_test.vcf.gz] -o [out.haps]

**count_bcs**

	gemtools2 -T count_bcs -i [out.shared] -b [LR.bam] -x [in_window] -y [out_window] -s [sv_name] -q [all|shared] -o [out.bc_count] 

**plot_hmw** (see above also)

	gemtools2 -T plot_hmw -i [out.bc_count] -o [out.pdf]
