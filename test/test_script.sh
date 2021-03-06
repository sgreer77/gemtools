#!/bin/bash

SV_FILE="test_files/HCC1954_chr9_svs.bedpe"
BAM_FILE="test_files/HCC1954_subset.bam"
VCF_FILE="test_files/HCC1954_subset.vcf.gz"
REFINED_FILE="test_files/refined_regions.txt" # user must generate this if they want to refine barcode selection

echo -e "Testing package functions -- this should take about a minute.\n"

# Get shared bcs
echo "Testing set_bc_window..."
gemtools -T set_bc_window -i $SV_FILE -o svs.wndw.bed -w 100000 -m auto
echo "Testing get_shared_bcs..."
gemtools -T get_shared_bcs -i svs.wndw.bed -b $BAM_FILE -o svs.shared.txt

# Plotting
#echo "Testing count_bcs..."
#gemtools -T count_bcs -i svs.shared.txt -b $BAM_FILE -x 1000 -y 500000 -s 'call_2080,call_440,call_189' -q shared -o svs.bc_count.txt
#echo "Testing plot_hmw..."
#gemtools -T plot_hmw -i svs.bc_count.txt -o svs.pdf

# Assign haps with shared barcodes
#echo "Testing assign_sv_haps..."
#gemtools -T assign_sv_haps -i svs.shared.txt -v $VCF_FILE -o svs.haps.txt -w 1000000 -q shared

echo "Testing count_bcs..."
#gemtools -T count_bcs -i svs.shared_refined.txt -b $BAM_FILE -x 1000 -y 500000 -s 'call_2080,call_440,call_189' -q select -o svs.bc_count_refined.txt
gemtools -T count_bcs -i $SV_FILE -e svs.shared.txt -b $BAM_FILE -x 10000 -y 300000 -s call_2080 -o svs.bc_count.txt
echo "Testing plot_hmw..."
gemtools -T plot_hmw -i svs.bc_count.txt -o call_2080.png

# Assign haps with select_barcodes -- better than before, with shared barcodes!
echo "Testing set_hap_window..."
gemtools -T set_hap_window -i $SV_FILE -o svs.hap_wndw.txt -w 100000
echo "Testing assign_sv_haps..."
gemtools -T assign_sv_haps -i svs.hap_wndw.txt -e svs.shared.txt -v $VCF_FILE -o svs.haps.txt
#gemtools -T assign_sv_haps -i svs.shared_refined.txt -v $VCF_FILE -o svs.haps_refined.txt -w 1000000 -q select

# Get phase block information
echo "Testing get_phased_basic..."
gemtools -T get_phased_basic -v $VCF_FILE -o phased_basic.txt
echo "Testing get_phase_blocks..."
gemtools -T get_phase_blocks -i phased_basic.txt -o phase_blocks.txt

# General tools
echo "Testing get_phased_bcs..."
gemtools -T get_phased_bcs -i phased_basic.txt -p 123859090 -o phased_bcs.txt
echo "Testing get_bcs_in_region (single region)..."
gemtools -T get_bcs_in_region -b $BAM_FILE -f chr9,128200000,128300000 -o bcs_in_region.txt
echo "Testing get_bcs_in_region (multiple regions)..."
gemtools -T get_bcs_in_region -b $BAM_FILE -f 'chr9,128200000,128300000;chr9,128700000,128800000' -o bcs_in_regions.txt

echo "Testing count_bcs_list..."
gemtools -T count_bcs_list -b $BAM_FILE -f chr9,128200000,128300000 -x 1000 -l bcs_in_region.txt -o bc_count_from_list.txt
#echo "Testing plot_hmw (again)..."
#gemtools -T plot_hmw -i bc_count_from_list.txt -o list.pdf

# Subset fastq's
echo "Testing extract_reads_interleaved..."
gemtools -T extract_reads_interleaved --bc_list test_files/call_189_bcs.txt --fqdir test_files/fastq_subset --sample_bcs ACGACATT,CACGTCGG,GTATGTCA,TGTCAGAC --lanes 1,2 --outdir fastq_call_189

echo "Testing complete!"
