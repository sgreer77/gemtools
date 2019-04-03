#!/usr/bin/env python


"""
gemtools version %version
Copyright (C) 2018 Stephanie Greer <sgreer2@stanford.edu>
gemtools is a suite of tools to work with linked read sequencing data.
Usage:
    gemtools -T TOOL [options]
"""

import os
import sys
import argparse

from gemtools import __version__
from gemtools.set_bc_window_f import set_bc_window
from gemtools.get_shared_bcs_f import get_shared_bcs
from gemtools.set_hap_window_f import set_hap_window
from gemtools.assign_sv_haps_f import assign_sv_haps
from gemtools.count_bcs_f import count_bcs
from gemtools.plot_hmw_f import plot_hmw
from gemtools.get_phased_basic_f import get_phased_basic
from gemtools.get_phase_blocks_f import get_phase_blocks
from gemtools.get_bcs_in_region_f import get_bcs_in_region
from gemtools.get_phased_bcs_f import get_phased_bcs
from gemtools.count_bcs_list_f import count_bcs_list
from gemtools.extract_reads_interleaved_f import extract_reads_interleaved
from gemtools.extract_reads_separate_f import extract_reads_separate

from gemtools.get_hmw_summary_f import get_hmw_summary
from gemtools.align_contigs_f import align_contigs
from gemtools.assess_contigs_f import assess_contigs

def gt_usage_msg(name=None):                                                            
    return '''\tgemtools -T <sub-tool> [options]
        '''

gt_help_msg = """gemtools: flexible tools for linked read sequencing (10X Genomics) analysis.
usage:	gemtools -T <sub-tool> [options]\n
The gemtools sub-tools include:\n 
[ Phase blocks ] 
    get_phased_basic   Obtain phasing information for all SNVs in the vcf file. 
    get_phase_blocks   Summarize phase blocks -- coordinates, size, SNVs per phase block etc.\n
[ SV analysis tools ]
    set_bc_window	Generate windows around SV breakpoints for SV analysis.
    get_shared_bcs	Determine barcodes shared between SV breakpoints.
    set_hap_window	Generate windows around SV breakpoints for haplotype analysis.
    assign_sv_haps	Assign SV barcodes to existing haplotypes (SNVs).
    count_bcs		Determine presence and quantity of given barcodes across a given region surrounding the SV breakpoints.
    plot_hmw		Generate a plot of the mapping locations of reads with each barcode. \n
[ Subset reads by barcode ]
    extract_reads_separate	Obtain reads with particular barcodes from Long Ranger fastq files (R1,R2,I1).
    extract_reads_interleaved	Obtain reads with particular barcodes from Long Ranger fastq files (RA,I1,I2). \n
[ General tools ]
    get_phased_bcs	For a particular phase block, return the haplotype 1 and haplotype 2 barcodes.
    get_bcs_in_region	Get all the barcodes that exist in a given region of the genome.
    count_bcs_list	Determine presence and quantity of given barcodes across a given region.
    plot_hmw		Generate a plot of the mapping locations of reads with each barcode (SAME AS ABOVE).
        """
# from SV analysis tools
#    align_contigs	Align de novo assembled contigs to a genome reference (using mappy/minimap2).
#    assess_contigs	Determine which contigs may contain rearranged structures. \n


def get_option_parser():
	parser = argparse.ArgumentParser(description='Gemtools.', add_help=False, usage=gt_usage_msg())
	
	parser.add_argument("-T", "--tool", default=None, dest="tool",
		help="Name of tool to use ")
	parser.add_argument("-o", "--output", metavar="FILE",
		dest="outfile",
		help="Name of output file")
	parser.add_argument("-h", "--help",
		dest="help", help="Print this help menu", action="store_true")

	parser.add_argument("-i", "--input", metavar="FILE",
		dest="infile",
		help="Name of input file ")
	parser.add_argument("-b","--bam", metavar="FILE",
		dest="bam",
		help="bam file")
	parser.add_argument("-v","--vcf", metavar="FILE",
		dest="vcf",
		help="vcf file")
	parser.add_argument("-c","--vcf_control", metavar="FILE",
		dest="vcf_control",
		help="Control vcf file")
	parser.add_argument("-w", "--window", type=int,
		dest="window_size", metavar="WINDOW",
		help="Size of window to create around bkpts in bp      ")
	parser.add_argument("-x","--in_window", type=int, default=1000,
		dest="in_window",metavar="WINDOW",
		help="Size of small windows in bp                       "
			"default: 1000")
	parser.add_argument("-y","--out_window", type=int, default=50000,
		dest="out_window",metavar="WINDOW",
		help="Size of large window in bp                       "
			"default: 50000")
	parser.add_argument("-f","--region_in",
		dest="region_in", metavar='REGION',
		help="In region(s) in format: "
		"chr1,1000000,2000000 or "
		"chr1,1000000,2000000;chr2:3000000,4000000")
	parser.add_argument("-g","--region_out",
		dest="region_out", metavar='REGION',
		help="Out region(s) in format: "
		"chr1,1000000,2000000 or "
		"chr1,1000000,2000000;chr2:3000000,4000000")
	parser.add_argument("-l","--bc_list", metavar="FILE",
		dest="bcs",
		help="File with list of barcodes")			
	parser.add_argument("-q","--bc_select",metavar='(all|shared|select)',
		dest="bc_select",choices=('all','shared','select'), default="shared",
		help="BCs to consider: all bcs or shared bcs               "
			"default: shared")
	parser.add_argument("-n","--chrom", metavar="CHR",
		dest="chrom", default="None",
		help="Chromosome number; ex: 'chr22','22'")
	parser.add_argument("-p","--phase_block",
		dest="phase_block", metavar="PHASE_ID",
		help="Phase block id (from vcf)")
	parser.add_argument("-s","--sv_name",
		dest="sv_name", metavar="SV",
		help="Name of SV; ex: 'call_144', '144'")
	parser.add_argument("-d","--fqdir",
		dest="fqdir", metavar="FQ_DIR",
		help="Long Rnager fastq directory")
	parser.add_argument("-j","--sample_bcs",
		dest="s_bcs", metavar="SAMPLE_BCS",
		help="Sample barcodes")
	parser.add_argument("-k","--lanes",
		dest="lanes", metavar="LANES",
		help="Numbers of sequencing lanes; ex: '1,5'")
	parser.add_argument("-z","--outdir",
		dest="outdir", metavar="OUT_DIR",
		help="Name of output directory")
	parser.add_argument("--read1",
		dest="read1", metavar="READ1",
		help="read1 fastq file")
	parser.add_argument("--read2",
		dest="read2", metavar="READ2",
		help="read2 fastq file")
	parser.add_argument("--index1",
		dest="index1", metavar="INDEX1",
		help="index1 fastq file")
	parser.add_argument("-m","--region_mode",metavar='(auto|window)',
		choices=('auto','window'),
		dest="region_mode",
		help="Mode to generate regions: auto or window")
	parser.add_argument("-e", "--shared_bc_file", metavar="FILE",
		dest="shrd_file",
		help="File of shared bcs")
	parser.add_argument("-r", "--ref", metavar="FILE",
		dest="ref_file",
		help="File of genome reference")
        parser.add_argument("--sort",
                dest="sort", help="Sort the barcodes by start coordinate", action="store_true")

	return parser

def pipeline_from_parsed_args(args):
	if args.tool=="set_hap_window":
		pipeline = set_hap_window(bedpe=args.infile, window=args.window_size, out=args.outfile)
	if args.tool=="get_shared_bcs":
		pipeline = get_shared_bcs(bed_in=args.infile, bam=args.bam, out=args.outfile)
	if args.tool=="assign_sv_haps":
		pipeline = assign_sv_haps(sv=args.infile, window=args.window_size, vcf_control=args.vcf_control, vcf_test=args.vcf, out=args.outfile, shrd_file = args.shrd_file)
	if args.tool=="count_bcs":
		pipeline = count_bcs(bam=args.bam, sv=args.infile, in_window=args.in_window, out_window=args.out_window, sv_name=args.sv_name, out=args.outfile, shrd_file = args.shrd_file)
	if args.tool=="get_phased_basic":
		pipeline = get_phased_basic(vcf=args.vcf, out=args.outfile, chrom=args.chrom)
	if args.tool=="get_phase_blocks":
		pipeline = get_phase_blocks(infile_basic=args.infile, out=args.outfile)
	if args.tool=="get_phased_bcs":
		pipeline = get_phased_bcs(infile_basic=args.infile, ps=args.phase_block, out=args.outfile)
	if args.tool=="count_bcs_list":
		pipeline = count_bcs_list(region=args.region_in, in_window=args.in_window, bam=args.bam, bcs=args.bcs, out=args.outfile)
	if args.tool=="get_bcs_in_region":
		pipeline = get_bcs_in_region(region=args.region_in,bam=args.bam, out=args.outfile)
	if args.tool=="plot_hmw":
		pipeline = plot_hmw(in_windows=args.infile, out=args.outfile, sort_by_coord=args.sort)
	if args.tool=="extract_reads_interleaved":
		pipeline = extract_reads_interleaved(fqdir=args.fqdir, s_bcs=args.s_bcs, lanes=args.lanes, bcs=args.bcs, fq_outdir=args.outdir)
	if args.tool=="extract_reads_separate":
		pipeline = extract_reads_separate(bcs=args.bcs, fq_outdir=args.outdir, read1=args.read1, read2=args.read2, index1=args.index1)
	if args.tool=="align_contigs":
		pipeline = align_contigs(infile_fasta=args.infile, genome=args.ref_file, out=args.outfile)
	if args.tool=="assess_contigs":
		pipeline = assess_contigs(infile_aln=args.infile, out=args.outfile)
	if args.tool=="set_bc_window":
		pipeline = set_bc_window(bedpe=args.infile, window=args.window_size, out=args.outfile, mode=args.region_mode)
	return pipeline

def main(cmdlineargs=None):
	parser = get_option_parser()

	args = parser.parse_args()
	
	if not args.tool:
		print gt_help_msg
		sys.exit(1)
		if args.help:
			print gt_help_msg
			sys.exit(1)

	if args.tool not in ['get_phased_basic','get_phase_blocks','set_bc_window','get_shared_bcs','set_hap_window','assign_sv_haps','count_bcs','plot_hmw','extract_reads_separate','extract_reads_intererleaved','get_phased_bcs','get_bcs_in_region','count_bcs_list','plot_hmw']:
		print "Please provide a valid gemtools sub-tool.\n"
		print gt_help_msg
		sys.exit(1)

##########################################################################################		
	if args.tool=="set_hap_window":
		#print "gemtools -T set_bc_window -i [LR_input.bedpe] -w [window_size] -o [out.bedpe]"
		if args.help:
			print """
Tool:	gemtools -T set_hap_window
Summary: Generate windows around SV breakpoints for SV analysis\n
Usage:   gemtools -T set_hap_window [OPTIONS] -i <LR_input.bedpe> -o <out.bedpe> -w <window_size>
Input:
	-i  bedpe file of SV breakpoints (ex: sv_call.bedpe from Long Ranger)
Output:
	-o  output file: bedpe file with windows around breakpoints
Options:
	-w  size of window to generate around the breakpoints; if supplied, windows of this specified size will be generated for each SV breakpoint
			"""
			sys.exit(1)
		if not (args.infile or args.outfile):
			parser.error('Missing required input')
		if not args.window_size:
			parser.error('Missing required input')
		
		if not os.path.isfile(args.infile):
			parser.error(str(args.infile) + " does not exist")
	
##########################################################################################
	if args.tool=="get_shared_bcs":
		#print "gemtools -T get_shared_bcs -i [out.bedpe] -b [LR_bam_file] -o [out.shared]"
		if args.help:
			print """
Tool:	gemtools -T get_shared_bcs
Summary: Determine barcodes shared between SV breakpoints\n
Usage:   gemtools -T get_shared_bcs -i <out.bedpe> -b <LR_bam_file> -o <out.shared.txt>
Input:
	-i  output file from 'set_bc_window' tool
	-b  bam file generated by Long Ranger
Output:
	-o  output file: List and count of SV-specific barcodes for each SV event
			"""
			sys.exit(1)
		if not (args.infile or args.outfile or args.bam):
			parser.error('Missing required input')

		if not os.path.isfile(args.infile):
			parser.error(str(args.infile) + " does not exist")
		if not os.path.isfile(args.bam):
			parser.error(str(args.bam) + " does not exist")
		if not str(args.bam).endswith(".bam"):
			parser.error(str(args.bam) + " does not appear to be a bam file")

##########################################################################################
	if args.tool=="assign_sv_haps":
		#print "gemtools -T assign_sv_haps -i [out.shared] -c [LR_control.vcf.gz] -t [LR_test.vcf.gz] -o [out.haps]"
		if args.help:
			print """
Tool:	gemtools -T assign_sv_haps
Summary: Assign SV barcodes to existing haplotypes (SNVs)\n
Usage:   gemtools -T assign_sv_haps [OPTIONS] -i <out.shared.txt> -v <LR_test.vcf.gz> -o <out.haps.txt> -q <shared/select>
Input:
	-i  output file from 'get_shared_bcs' tool
	-v  vcf file generated by Long Ranger for test sample (ex: tumor sample)
	-e  name of file generated by 'get_shared_bcs'
Output:
	-o  output file: List of breakpoints with phase id and number of barcodes supporting assignment to each haplotype
Options:
	-c  vcf file generated by Long Ranger for control sample (ex: normal sample); only specify if you want to use a different vcf to define phase blocks
			"""
			sys.exit(1)
		if not (args.infile or args.outfile or args.vcf or args.shrd_file):
			parser.error('Missing required input')

		if str(args.vcf_control)!="None":
			if not str(args.vcf_control).endswith(".vcf.gz"):
				parser.error(str(args.vcf_control) + " does not appear to be a gzipped vcf file")
		if not str(args.vcf).endswith(".vcf.gz"):
			parser.error(str(args.vcf) + " does not appear to be a gzipped vcf file")

##########################################################################################
	if args.tool=="count_bcs":
		#print "gemtools -T count_bcs -i [out.shared] -b [LR.bam] -x [in_window] -y [out_window] -s [sv_name] -q [all|shared] -o [out.bc_count]"
		if args.help:
			print """
Tool:	gemtools -T count_bcs
Summary: Determine presence and quantity of given barcodes across a given region surrounding the SV breakpoints\n
Usage:   gemtools -T count_bcs [OPTIONS] -i <out.shared.txt> -b <LR.bam> -s <sv_name> -q <all/shared/select> -o <out.bc_count.txt>
Input:
	-i  output file from 'get_shared_bcs' tool (can also be merged output from 'refine_bcs' tool)
	-b  bam file generated by Long Ranger
	-s  name(s) of the SV(s) to check; if multiple: comma-separated list
Output:
	-o  output file: barcode counts in windows
Options:
	-q  define whether to check all barcodes for an SV, only the shared barcodes, or the select barcodes (default: shared)
	-x  size of small windows to check for barcodes (default: 1000 bp)
	-y  size of large windows around breakpoints to check for barcodes (default: 500000 bp)
			"""
			sys.exit(1)
		if not (args.infile or args.outfile or args.sv_name or args.bam or args.shrd_file):
			parser.error('Missing required input')

		if not os.path.isfile(args.infile):
			parser.error(str(args.infile) + " does not exist")
		if not(str(args.in_window).isdigit() and int(args.in_window)>0):
			parser.error(str(args.in_window) + " must be an integer >0")
		if not (str(args.out_window).isdigit() and int(args.out_window)>0):
			parser.error(str(args.out_window) + " must be an integer >0")
		if not str(args.bam).endswith(".bam"):
			parser.error(str(args.bam) + " does not appear to be a bam file")

##########################################################################################
	if args.tool=="get_phased_basic":
		#print "gemtools -T get_phased_basic -v [LR.vcf.gz] -o [output.phased_basic] -n [chr_num]"
		if args.help:
			print """
Tool:	gemtools -T get_phased_basic
Summary: Obtain phasing information for all SNVs in the vcf file\n
Usage:   gemtools -T get_phased_basic [OPTIONS] -v <LR.vcf.gz> -o <output.phased_basic.txt>
Input:
	-v  gzipped vcf file output from Long Ranger
Output:
	-o  output file: each row is an SNV; columns are phasing information for each SNV
Options:
	-n  chromosome number (ex: 22 or chr22)
			"""
			sys.exit(1)
		if not (args.outfile or args.vcf):
			parser.error('Missing required input')

		if not str(args.vcf).endswith(".vcf.gz"):
			parser.error(str(args.vcf) + " does not appear to be a gzipped vcf file")

##########################################################################################	
	if args.tool=="get_phase_blocks":
		#print "gemtools -T get_phase_blocks -i [out.phased_basic] -o [out.phase_blocks]"
		if args.help:
			print """
Tool:	gemtools -T get_phase_blocks
Summary: Summarize phase blocks\n
Usage:   gemtools -T get_phase_blocks -i <output.phased_basic.txt> -o <output.phased_blocks.txt>
Input:
	-i  output from 'get_phased_basic' tool
Output:
	-o  output file: each row is a phase block, columns summarize information for each phase block (size etc.)
			"""
			sys.exit(1)
		if not (args.infile or args.outfile):
			parser.error('Missing required input')

		if not os.path.isfile(args.infile):
			parser.error(str(args.infile) + " does not exist")

##########################################################################################
	if args.tool=="get_bcs_in_region":
		#print "gemtools -T get_bcs_in_region -b [LR.bam] -f [region_in] -o [out.bcs]"
		if args.help:
			print """
Tool:	gemtools -T get_bcs_in_region
Summary: Get all the barcodes that exist in a given region of the genome\n
Usage:   gemtools -T get_bcs_in_region -b <LR.bam> -f <region> -o <out.bcs.txt>
Input:
	-b  bam file generated by Long Ranger
	-f  region(s) where barcodes must be located; format 'chr1,1000,2000' or 'chr1,1000,2000;chr1,3000,4000'
Output:
	-o  output file: list of barcodes
			"""
			sys.exit(1)
		if not (args.outfile or args.bam or args.region_in):
			parser.error('Missing required input')
		
		if not str(args.bam).endswith(".bam"):
			parser.error(str(args.bam) + " does not appear to be a bam file")

##########################################################################################	
	if args.tool=="get_phased_bcs":
		#print "gemtools -T get_phased_bcs -i [out.phased_basic] -p [phase_block_id] -o [out.phased_bcs]"
		if args.help:
			print """
Tool:	gemtools -T get_phased_bcs
Summary: For a particular phase block, return the haplotype 1 and haplotype 2 barcodes\n
Usage:   gemtools -T get_phased_bcs -i <output.phased_basic.txt> -p <phase_block_id> -o <output.phased_bcs.txt>
Input:
	-i  output from 'get_phased_basic' tool
	-p  id number for phase block of interest
Output:
	-o  output file: a table with the haplotype 1 and haplotype 2 barcodes indicated
			"""
			sys.exit(1)
		if not (args.infile or args.outfile or args.phase_block):
			parser.error('Missing required input')

		if not os.path.isfile(args.infile):
			parser.error(str(args.infile) + " does not exist")

##########################################################################################
#	if args.tool=="select_bcs":
#		#print "gemtools -T select_bcs -b [LR.bam] -f [region_in] -g [region_out] -o [out.bcs]"
#		if args.help:
#			print """
#Tool:	gemtools -T select_bcs
#Summary: Get barcodes shared by ALL region_in's and present in NONE of the region_out's\n
#Usage:   gemtools -T select_bcs -b <LR.bam> -f <region_in> -g <region_out> -o <output.bcs.txt>
#Input:
#	-i  bam file generated by Long Ranger
#	-f  region(s) where barcode must be located
#	-g  region(s) where barcode must NOT be located
#Output:
#	-o  output file: list of barcodes
#			"""
#			sys.exit(1)
#		if not (args.bam or args.outfile or args.region_in or args.region_out):
#			parser.error('Missing required input')
#		
#		if not str(args.bam).endswith(".bam"):
#			parser.error(str(args.bam) + " does not appear to be a bam file")
#"""
##########################################################################################
	if args.tool=="count_bcs_list":
		#print "gemtools -T count_bcs_list -b [LR.bam] -f [region_in] -x [in_window] -l [bc_list] -o [out.bc_count]"
		if args.help:
			print """
Tool:	gemtools -T count_bcs_list
Summary: Determine presence and quantity of given barcodes across a given region\n
Usage:   gemtools -T count_bcs_list -b <LR.bam> -f <region_in> -x <in_window> -l <bc_list> -o <output.bc_count.txt>
Input:
	-b  bam file generated by Long Ranger
	-f  region(s) to assess barcodes
	-x  size of windows to check for barcodes
	-l  file containing list of barcodes (one barcode per line)
Output:
	-o  output file: barcode counts in windows 
			"""
			sys.exit(1)		
		if not (args.bam or args.region_in or args.outfile or args.bcs):
			parser.error('Missing required input')

		if not str(args.bam).endswith(".bam"):
			parser.error(str(args.bam) + " does not appear to be a bam file")
		if not (str(args.in_window).isdigit() and int(args.in_window)>0):
			parser.error(str(args.in_window) + " must be an integer >0")
		if not os.path.isfile(args.bcs):
			parser.error(str(args.bcs) + " does not exist")

##########################################################################################
	if args.tool=="plot_hmw":
		#print "gemtools -T plot_hmw -i [out.bc_count] -o [out.pdf]"
		if args.help:
			print """
Tool:	gemtools -T plot_hmw
Summary: Generate a plot of the mapping locations of reads with each barcode\n
Usage:   gemtools -T plot_hmw -i <output.bc_count.txt> -o <output.pdf>
Input:
	-i  output file generated by 'count_bcs' or 'count_bcs_list' tool
	--sort sort the barcodes by mapping coordinate (optional)
Output:
	-o  output file: plot of barcode mapping locations in a given region
			"""
			sys.exit(1)		
		if not (args.infile or args.outfile):
			parser.error('Missing required input')
		
		if not os.path.isfile(args.infile):
			parser.error(str(args.infile) + " does not exist")
                if not str(args.outfile).endswith(".png"):
                        parser.error(str(args.outfile) + " : Output file must be a png file, ex: out.png")

##########################################################################################
	if args.tool=="extract_reads_interleaved":
		#print "gemtools -T extract_reads_interleaved -l [bc_list] -z [fastq_output_dir] -d [LR_fastq_dir] -j [sample_barcodes] -k [sample_lanes]"
		if args.help:
			print """
Tool:	gemtools -T extract_reads_interleaved
Summary: Obtain reads with particular barcodes from Long Ranger fastq files (RA,I1,I2)\n
Usage:   gemtools -T extract_reads_interleaved -l <bc_list> -z <fastq_output_dir> -d <LR_fastq_dir> -j <sample_barcodes> -k <sample_lanes>
Input:
	-l  file containing list of barcodes (one barcode per line)
	-d  Long Ranger fastq directory, containing RA and I1 fastq files
	-j  Long Ranger sample barcodes
	-k  seq lanes to consider
Output:
	-z  Output directory for output fastq files; subsetted RA and I1 files will be generated here
			"""
			sys.exit(1)		
		if not (args.fqdir or args.s_bcs or args.lanes or args.bcs or args.outdir):
			parser.error('Missing required input')
		
		if not os.path.isfile(args.bcs):
			parser.error(str(args.bcs) + " does not exist")
		if os.path.isdir(args.outdir):
			parser.error(str(args.outdir) + " already exists")

##########################################################################################	
	if args.tool=="extract_reads_separate":
		#print "gemtools -T extract_reads_separate -l [bc_list] -z [fastq_output_dir] --read1 [LR_R1.fastq.gz] --read2 [LR_R2.fastq.gz] --index1 [LR_I1.fastq.gz]"
		if args.help:
			print """
Tool:	gemtools -T extract_reads_separate
Summary: Obtain reads with particular barcodes from Long Ranger fastq files (R1,R2,I1)\n
Usage:   gemtools -T extract_reads_separate -l <bc_list> -z <fastq_output_dir> --read1 <LR_R1.fastq.gz> --read2 <LR_R2.fastq.gz> --index1 <LR_I1.fastq.gz>
Input:
	-l  file containing list of barcodes (one barcode per line) 
	--read1  Long Ranger read 1 fastq
	--read2  Long Ranger read 2 fastq
	--index1  Long Ranger index 1 fastq  
Output:
	-z  Output directory for output fastq files; subsetted R1, R2 and I1 files will be generated here
			"""
			sys.exit(1)		
		if not (args.bcs or args.read1 or args.read2 or args.index1 or args.outdir):
			parser.error('Missing required input')
			
		if not os.path.isfile(args.bcs):
			parser.error(str(args.bcs) + " does not exist")

##########################################################################################	
	if args.tool=="align_contigs":
		if args.help:
			print """
Tool:	gemtools -T align_contigs
Summary: Aligns contigs to the genome\n
Usage:   gemtools -T align_contigs [OPTIONS] -i <de_novo.fasta.gz> -o <out.txt> -r <genome_reference_fasta>
Input:
	-i  fasta file of de novo contigs (ex: output of supernova)
	-r  genome reference fasta file
Output:
	-o  output file: alignment info
			"""
			sys.exit(1)
		if not (args.infile or args.outfile or args.ref_file):
			parser.error('Missing required input')
		
		if not os.path.isfile(args.infile):
			parser.error(str(args.infile) + " does not exist")
		if not os.path.isfile(args.ref_file):
			parser.error(str(args.ref_file) + " does not exist")
	
##########################################################################################
	if args.tool=="assess_contigs":
		if args.help:
			print """
Tool:	gemtools -T assess_contigs
Summary: Assess whether contigs have rearranged structures\n
Usage:   gemtools -T assess_contigs [OPTIONS] -i <aligned_contigs.txt> -o <out.txt>
Input:
	-i  aligned contigs file (output of align_contigs)
Output:
	-o  output file: alignment info marked with assessment
			"""
			sys.exit(1)
		if not (args.infile or args.outfile):
			parser.error('Missing required input')
		
		if not os.path.isfile(args.infile):
			parser.error(str(args.infile) + " does not exist")
	
##########################################################################################
	if args.tool=="set_bc_window":
		#print "gemtools -T bedpe2window -i [LR_input.bedpe] -w [window_size] -o [out.bedpe]"
		if args.help:
			print """
Tool:	gemtools -T set_bc_regions
Summary: Generate windows for SV analysis\n
Usage:   gemtools -T set_bc_regions [OPTIONS] -i <LR_input.bedpe> -o <out.bedpe> -w <window_size> -m <region_mode: auto or window>
Input:
	-i  bedpe file of SV breakpoints (ex: sv_call.bedpe from Long Ranger)
	-m	mode to run: auto or window (if auto: bedpe file must include an 'info' column with 'TYPE=' defined for each event)
Output:
	-o  output file: bedpe file with windows around breakpoints
Options:
	-w  size of window to generate around the breakpoints; if supplied, windows of this specified size will be generated for each SV breakpoint
			"""
			sys.exit(1)
		if not (args.infile or args.outfile):
			parser.error('Missing required input')
		if not args.window_size:
			parser.error('Missing required input')
		if not args.region_mode:
			parser.error('Missing required input')
		
		if not os.path.isfile(args.infile):
			parser.error(str(args.infile) + " does not exist")
	
##########################################################################################
	pipeline = pipeline_from_parsed_args(args)
	runner = pipeline

if __name__ == '__main__':
	main()
