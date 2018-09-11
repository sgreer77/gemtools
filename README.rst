gemtools2
---------

Tools for working with linked-read sequencing data -- performing barcode analyses and SV structural analyses. 

Installation
============

Dependencies - the following python and R packages are required:

Python:
pandas
numpy
pysam
pyvcf
rpy2

R:
ggplot2

To install gemtools2:
pip install git+https://github.com/sgreer77/gemtools2.git


Running gemtools2
=================

1. Run longranger
"""""""""""""""""

gemtools2 uses the output of longranger to perform additional analyses


2. Running tools in gemtools2
"""""""""""""""""""""""""""""
General usage:

**gemtools2 -T [tool] [-options]**

a) Getting basic information about the phase blocks:

**get_phased_basic**

gemtools2 -T get_phased_basic -v [vcf_file] -o [output.phased_basic]

**get_phase_blocks**

gemtools2 -T get_phase_blocks -i [output.phased_basic] -o [output.phase_blocks]


b) Generally useful tools:

**get_phased_bcs**

gemtools2 -T get_phased_bcs -i [output.phased_basic] -p [phase_block_id] -o [out.phased_bcs]

**select_bcs**

gemtools2 -T select_bcs -b [bam_file] -f [region_in] -g [region_out] -o [out.bcs]

**get_bcs_in_region**

gemtools2 -T get_bcs_in_region -b [bam_file] -f [region_in] -o [out.bcs]

**count_bcs_list**

gemtools2 -T count_bcs_list -b [bam_file] -f [region_in] -x [in_window] -b [bc_list.txt] -o [out.bc_count]

**plot_hmw**

gemtools2 -T plot_hmw -i [out.bc_count] -o [out.pdf]

c) SV analysis tools:

**bedpe2window**

gemtools2 -T bedpe2window -i [input.bedpe] -w [window_size] -o [out.bed]

**get_shared_bcs**

gemtools2 -T get_shared_bcs -i [out.bed] -b [bam_file] -o [out.shared]

**assign_sv_haps**

gemtools2 -T assign_sv_haps -i [out.shared] -c [control.vcf] -t [test.vcf] -o [out.haps]

**count_bcs**

gemtools2 -T count_bcs -i [out.shared] -b [bam_file] -x [in_window] -y [out_window] -s [sv_name] -q [all|shared] -o [out.bc_count] 

**plot_hmw** (see above also)

gemtools2 -T plot_hmw
