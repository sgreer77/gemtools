gemtools2
---------

Tools for working with linked-read sequencing data -- performing barcode analyses and SV structural analyses. 

Installation
============

Dependencies - the following python and R packages are required:

Python:_
pandas_
numpy_
pysam_
pyvcf_
rpy2_

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

gemtools2 -T 

**select_bcs**

gemtools2 -T 

**get_bcs_in_region**

gemtools2 -T 

**count_bcs_list**

gemtools2 -T 

**plot_hmw**

gemtools2 -T 

c) SV analysis tools:

**bedpe2window**

gemtools2 -T 

**get_shared_bcs**

gemtools2 -T 

**assign_sv_haps**

gemtools2 -T 

**count_bcs**

gemtools2 -T 

**plot_hmw** (see above also)

gemtools2 -T 
