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

Long Ranger output files are indicated in the instructions below with an **LR**

General usage: **gemtools2 -T [tool] [-options]**
"""""""""""""""""""""""""""""""""""""""""""""""""

Tools for getting basic information about the phase blocks:
""""""""""""""""""""""""""""""""""""""""""""""""""""

**get_phased_basic**

	gemtools2 -T get_phased_basic -v [LR_vcf_file] -o [output.phased_basic]
	Ex: gemtools2 -T get_phased_basic -v phased_variants.vcf.gz -o output.phased_basic
	Input:
		-v gzipped vcf file output from LR
	Output:
		-o output.phased_basic

**get_phase_blocks**

	gemtools2 -T get_phase_blocks -i [output.phased_basic] -o [output.phase_blocks]
	Ex: gemtools2 -T get_phase_blocks -i output.phased_basic -o output.phase_blocks
	Input:
		-i output from 'get_phased_basic' tool
	Output:
		-o output file

Generally useful tools:
""""""""""""""""""""""""""

**get_phased_bcs**

	gemtools2 -T get_phased_bcs -i [output.phased_basic] -p [phase_block_id] -o [out.phased_bcs]

**select_bcs**

	gemtools2 -T select_bcs -b [LR_bam_file] -f [region_in] -g [region_out] -o [out.bcs]

**get_bcs_in_region**

	gemtools2 -T get_bcs_in_region -b [LR_bam_file] -f [region_in] -o [out.bcs]

**count_bcs_list**

	gemtools2 -T count_bcs_list -b [LR_bam_file] -f [region_in] -x [in_window] -b [bc_list.txt] -o [out.bc_count]

**plot_hmw**

	gemtools2 -T plot_hmw -i [out.bc_count] -o [out.pdf]

SV analysis tools:
"""""""""""""""""""""

**bedpe2window**

	gemtools2 -T bedpe2window -i [LR_input.bedpe] -w [window_size] -o [out.bed]

**get_shared_bcs**

	gemtools2 -T get_shared_bcs -i [out.bed] -b [LR_bam_file] -o [out.shared]

**assign_sv_haps**

	gemtools2 -T assign_sv_haps -i [out.shared] -c [LR_control.vcf] -t [LR_test.vcf] -o [out.haps]

**count_bcs**

	gemtools2 -T count_bcs -i [out.shared] -b [LR_bam_file] -x [in_window] -y [out_window] -s [sv_name] -q [all|shared] -o [out.bc_count] 

**plot_hmw** (see above also)

	gemtools2 -T plot_hmw -i [out.bc_count] -o [out.pdf]
