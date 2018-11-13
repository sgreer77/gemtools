#gemtools
---------

##Installation
============

**Dependencies:**
"""""""""""""""""
**Python packages:** pandas, numpy, pysam, pyvcf, pybedtools, rpy2

**R package:** ggplot2

**To install gemtools:**
"""""""""""""""""""""""""
**If you already have the dependencies:**
::
	git clone https://github.com/sgreer77/gemtools.git
	cd gemtools
	pip install .

**If you do not already have the dependencies, a yaml file is provided in the package to generate a conda virtual env with all of the dependencies:**
                                
For this, you must have anaconda installed. To check if you have anaconda installed, run:

::
conda

If you do not receive the usage information then you can install anaconda/miniconda, following the instructions at: 
<https://conda.io/docs/user-guide/install/download.html>`_

::
	git clone https://github.com/sgreer77/gemtools.git
	cd gemtools
	conda env create -f gemtools_env.yml
	source activate gemtools_env
	pip install .

**Test the gemtools installation:**
::
	cd test
	./test_script.sh


##Running gemtools
=================

gemtools is a collection of tools that use the **output of Long Ranger** (10X Genomics) to perform additional analyses      (Long Ranger output files are indicated in the instructions below with an **LR** prefix)

**General usage: gemtools -T [tool] [-options]**


Tools for getting basic information about the phase blocks:
""""""""""""""""""""""""""""""""""""""""""""""""""""

**get_phased_basic**: Obtain phasing information for all SNVs in the vcf file

	gemtools -T get_phased_basic -v [LR.vcf.gz] -o [output.phased_basic] -n [chr_num]
	
	Ex: gemtools -T get_phased_basic -v phased_variants.vcf.gz -o output.phased_basic
	
	Ex: gemtools -T get_phased_basic -v phased_variants.vcf.gz -o output.phased_basic -n 22
	
	Input:
		-v gzipped vcf file output from LR
		
		-n chromosome (optional)
	Output:
		-o output file: each row is an SNV; columns are phasing information for each SNV

**get_phase_blocks**: Summarize phase blocks -- coordinates, size, number of phased heterozygous SNVs per phase block etc.

	gemtools -T get_phase_blocks -i [out.phased_basic] -o [out.phase_blocks]
	
	Ex: gemtools -T get_phase_blocks -i out.phased_basic -o out.phase_blocks
	
	Input:
		-i output from 'get_phased_basic' tool
	Output:
		-o output file: each row is a phase block, columns summarize information for each phase block (size etc.)


Generally useful tools:
""""""""""""""""""""""""""

**get_phased_bcs:** For a particular phase block, return the haplotype 1 and haplotype 2 barcodes

	gemtools -T get_phased_bcs -i [out.phased_basic] -p [phase_block_id] -o [out.phased_bcs]
	
	Ex: gemtools -T get_phased_bcs -i out.phased_basic -p 1356780 -o out.phased_bcs

	Input:
		-i output from 'get_phased_basic' tool
		
		-p id number for phase block of interest (phase block ids are originally assigned in the LR.vcf.gz file)
	Output:
		-o output file: a table with the haplotype 1 and haplotype 2 barcodes indicated
	
**get_bcs_in_region:** Get all the barcodes that exist in a given region(s) of the genome

	gemtools -T get_bcs_in_region -b [LR.bam] -f [region_in] -o [out.bcs]
	
	Ex: gemtools -T get_bcs_in_region -b phased_possorted.bam -f chr1,1000,2000 -o out.bcs.txt
	
	Ex: gemtools -T get_bcs_in_region -b phased_possorted.bam -f chr1,1000,2000;chr1,3000,4000 -o out.bcs.txt

	Input:
		-b bam file generated by Long Ranger
		
		-f region(s) where barcodes must be located (does not require that barcode be in all regions)
		
	Output:
		-o output file: list of barcodes

**count_bcs_list:** Determine presence and quantity of given barcodes across a given region

	gemtools -T count_bcs_list -b [LR.bam] -f [region_in] -x [in_window] -l [bc_list] -o [out.bc_count]
	
	Ex: gemtools -T count_bcs_list -b phased_possorted.bam -f chr1,1000,2000 -x 100 -l bc_list.txt -o out.bc_count.txt

	Input:
		-b bam file generated by Long Ranger
		
		-f region(s) to assess barcodes
		
		-x size of windows to check for barcodes
		
		-b file containing list of barcodes (one barcode per line)
		
	Output:
		-o output file: rows are genomic window coordinates, columns are each barcode in bc_list file, entries are number of each barcode in each window

**plot_hmw:** Generate a plot of the mapping locations of reads with each barcode

	gemtools -T plot_hmw -i [out.bc_count] -o [out.pdf]

	Input:
		-i output file generated by 'count_bcs_list' tool
		
	Output:
		-o output file: plot of barcode mapping locations in a given region

**refine_bcs:** Obtain barcodes based on where they do and do NOT map. 

	gemtools -T refine_bcs -i [bed file in/out regions specified] -b [LR.bam] -e [out.shared] -o [out.refined]
	
	Input:
		-i bed file of regions; order of columns must be: ['chrom','start','stop','name','status']; header line must be commented
			
			ex: #chrom	start	stop	name	status
				chr1	1000	2000	call_1	in
				chr1	2000	3000	call_1	out
				chr2	4000	5000	call_2	in
	
		-b bam file generated by Long Ranger
		
		-e name of file generated by 'get_shared_bcs'
		
	Output:
		-o output file: barcode info summary for each event (specified by 'name')

SV analysis tools:
"""""""""""""""""""""

**bedpe2window:** Generate windows around SV breakpoints for SV analysis

	gemtools -T bedpe2window [OPTIONS] -i [LR_input.bedpe] -o [out.bedpe] -m [mean_mol_size]
	
	Ex: gemtools -T bedpe2window -i large_sv_calls.bedpe -o large_sv_calls.wndw.bedpe -m 50000

	NOTE: User can specify only one of -m and -w; -m is recommended

	Input:
		-i bedpe file of SV breakpoints; this is typically the Long Ranger output: large_sv_calls.bedpe OR large_sv_candidates.bedpe
		
		-m size of HMW molecules input (this can be obtained from the Long Ranger 'summary.csv' file as 'molecule_length_mean')
		
		-w size of window to generate around the breakpoints
		
	Output:
		-o output file: bedpe file with windows around breakpoints

**get_shared_bcs:** Determine barcodes shared between SV breakpoints

	gemtools -T get_shared_bcs -i [out.bedpe] -b [LR_bam_file] -o [out.shared]
	
	Ex: gemtools -T get_shared_bcs -i large_sv_calls.wndw.bedpe -b phased_possorted.bam -o out.shared.txt
	
	Input:
		-i output file from 'bedpe2window' tool
		
		-b bam file generated by Long Ranger
		
	Output:
		-o output file: List and count of SV-specific barcodes for each SV event

**assign_sv_haps:** Assign SV barcodes to existing haplotypes (SNVs)

	gemtools -T assign_sv_haps -i [out.shared] -c [LR_control.vcf.gz] -t [LR_test.vcf.gz] -o [out.haps] -q [shared|select]
	
	Ex: gemtools -T assign_sv_haps -i out.shared.txt -v phased_variants.vcf.gz -c phased_variants.vcf.gz -o out.haps.txt -q shared
	
	Input:
		-i output file from 'get_shared_bcs' or 'refine_bcs' tool

		-v vcf file generated by Long Ranger for test sample (ex: tumor sample)
		
		-c vcf file generated by Long Ranger for control sample (ex: normal sample) -- this is optional, if the user wants to use a different vcf to define phase blocks

		-q define whether to check the shared or select barcodes of an SV
				
	Output:
		-o output file: List of breakpoints with phase id and number of barcodes supporting assignment to each haplotype

**count_bcs:** Determine presence and quantity of given barcodes across a given region surrounding the SV breakpoints

	gemtools -T count_bcs -i [out.shared] -b [LR.bam] -x [in_window] -y [out_window] -s [sv_name] -q [all|shared|select] -o [out.bc_count]
	
	Ex: gemtools -T count_bcs -i out.shared.txt -b phased_possorted.bam -x 1000 -y 50000 -s call_110 -q shared -o out.bc_count.txt 
	
	Input:
		-i output file from 'get_shared_bcs' or 'refine_bcs' tool
		
		-b bam file generated by Long Ranger
		
		-x size of small windows to check for barcodes
		
		-y size of large windows around breakpoints to check for barcodes
		
		-s name(s) of the SV(s) to check; if multiple, use a comma-separated list
		
		-q define whether to check all barcodes of an SV, only the shared barcodes, or the select barcodes
		
	Output:
		-o output file: rows are genomic window coordinates, columns are each barcode in bc_list file, entries are number of each barcode in each window

**plot_hmw:** Generate a plot of the mapping locations of reads with each barcode (SAME AS ABOVE)

	gemtools -T plot_hmw -i [out.bc_count] -o [out.pdf]

	Input:
		-i output file generated by 'count_bcs_list' tool
		
	Output:
		-o output file: plot of barcode mapping locations in a given region


Tools for extracting subset barcoded reads from fastq files:
""""""""""""""""""""""""""""""""""""""""""""""""""""

**extract_reads_separate**: Obtain reads with particular barcodes from Long Ranger fastq files (where fastq output is R1,R2,I1)

	gemtools -T extract_reads_separate -l [bc_list] -z [fastq_output_dir] --read1 [LR_R1.fastq.gz] --read2 [LR_R2.fastq.gz] --index1 [LR_I1.fastq.gz]
	
	Ex: gemtools -T extract_reads_separate -l bc_list.txt -z fastq_subset --read1 SAMPLE_S1_L001_R1_001.fastq.gz --read2 SAMPLE_S1_L001_R2_001.fastq.gz --index1 SAMPLE_S1_L001_I1_001.fastq.gz
	
	Input:
		-l file containing list of barcodes (one barcode per line)
		
		--read1 Long Ranger read 1 fastq
		
		--read2 Long Ranger read 2 fastq
		
		--index1 Long Ranger index 1 fastq
	Output:
		-z Output directory for output fastq files; subsetted R1, R2 and I1 files will be generated here

**extract_reads_interleaved**: Obtain reads with particular barcodes from Long Ranger fastq files (where fastq output is RA,I1,I2)

	gemtools -T extract_reads_interleaved -l [bc_list] -z [fastq_output_dir] -d [LR_fastq_dir] -j [sample_barcodes] -k [sample_lanes]
	
	Ex: gemtools -T extract_reads_interleaved -l bc_list.txt -z fastq_subset -d fastq -j 'ACGACGCT,CGCCATTC,GTAGTCAG,TATTGAGA' -k '1,5'
	
	Input:
		-l file containing list of barcodes (one barcode per line)
		
		-d Long Ranger fastq directory, containing RA and I1 fastq files
		
		-j Long Ranger sample barcodes
		
		-k seq lanes to consider
	Output:
		-z Output directory for output fastq files; subsetted RA and I1 files will be generated here
