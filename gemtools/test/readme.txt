
SV_FILE="HCC1954_chr9_svs.bedpe"
BAM_FILE="HCC1954_subset.bam"
VCF_FILE="HCC1954_subset.vcf.gz"
REFINED_FILE="in_n_out.txt" # user must generate this if they want to refine barcode selection

# Get shared bcs
gemtools -T bedpe2window -i $SV_FILE -o svs.wndw.bedpe -m 85000
gemtools -T get_shared_bcs -i svs.wndw.bedpe -b $BAM_FILE -o svs.shared.txt

# Plotting
gemtools -T count_bcs -i svs.shared.txt -b $BAM_FILE -x 1000 -y 500000 -s 'call_2080,call_440,call_189' -q shared -o svs.bc_count.txt
gemtools -T plot_hmw -i svs.bc_count.txt -o svs.pdf

# Assign haps with shared barcodes
gemtools -T assign_sv_haps -i svs.shared.txt -v $VCF_FILE -o svs.haps.txt -w 1000000 -q shared

# Refine barcodes then rerun plotting
gemtools -T refine_bcs -i $REFINED_FILE -b $BAM_FILE -o svs.shared_refined.txt -e svs.shared.txt
gemtools -T count_bcs -i svs.shared_refined.txt -b $BAM_FILE -x 1000 -y 500000 -s 'call_2080,call_440,call_189' -q select -o svs.bc_count_refined.txt
gemtools -T plot_hmw -i svs.bc_count_refined.txt -o svs_refined.pdf

# Assign haps with select_barcodes -- better than before, with shared barcodes!
gemtools -T assign_sv_haps -i svs.shared_refined.txt -v $VCF_FILE -o svs.haps_refined.txt -w 1000000 -q select

# Get phase block information
gemtools -T get_phased_basic -v $VCF_FILE -o phased_basic.txt
gemtools -T get_phase_blocks -i phased_basic.txt -o phase_blocks.txt

# General tools
gemtools -T get_phased_bcs -i phased_basic.txt -p 123859090 -o phased_bcs.txt
gemtools -T get_bcs_in_region -b $BAM_FILE -f chr9,128200000,128300000 -o bcs_in_region.txt
gemtools -T get_bcs_in_region -b $BAM_FILE -f 'chr9,128200000,128300000;chr9,128700000,128800000' -o bcs_in_regions.txt

gemtools -T count_bcs_list -b $BAM_FILE -f chr9,128200000,128300000 -x 1000 -l bcs_in_region.txt -o bc_count_from_list.txt
gemtools -T plot_hmw -i bc_count_from_list.txt -o list.pdf

# Subset fastq's
gemtools -T extract_reads_interleaved -l call_189_bcs.txt -d fastq_subset -j ACGACATT,CACGTCGG,GTATGTCA,TGTCAGAC -k 1,2,3,4,5,6,7,8 -z fastq_call_189
