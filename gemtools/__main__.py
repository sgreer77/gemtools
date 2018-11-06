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
from gemtools.bedpe2window_f import bedpe2window
from gemtools.get_shared_bcs_f import get_shared_bcs
from gemtools.assign_sv_haps_f import assign_sv_haps
from gemtools.count_bcs_f import count_bcs
from gemtools.plot_hmw_f import plot_hmw
from gemtools.get_phased_basic_f import get_phased_basic
from gemtools.get_phase_blocks_f import get_phase_blocks
from gemtools.refine_bcs_f import refine_bcs

from gemtools.get_bcs_in_region_f import get_bcs_in_region
from gemtools.get_phased_bcs_f import get_phased_bcs
#from gemtools.select_bcs_f import select_bcs
from gemtools.count_bcs_list_f import count_bcs_list

from gemtools.extract_reads_interleaved_f import extract_reads_interleaved
from gemtools.extract_reads_separate_f import extract_reads_separate

from gemtools.get_hmw_summary_f import get_hmw_summary


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
    bedpe2window	Generate windows around SV breakpoints for SV analysis.
    get_shared_bcs	Determine barcodes shared between SV breakpoints.
    refine_bcs		Select barcodes that are shared between regions and not in other regions.
    assign_sv_haps	Assign SV barcodes to existing haplotypes (SNVs).
    count_bcs		Determine presence and quantity of given barcodes across a given region surrounding the SV breakpoints.
    plot_hmw		Generate a plot of the mapping locations of reads with each barcode (SAME AS ABOVE). \n
[ Subset reads by barcode ]
    extract_reads_separate	Obtain reads with particular barcodes from Long Ranger fastq files (R1,R2,I1).
    extract_reads_interleaved	Obtain reads with particular barcodes from Long Ranger fastq files (RA,I1,I2). \n
[ General tools ]
    get_phased_bcs	For a particular phase block, return the haplotype 1 and haplotype 2 barcodes.
    get_bcs_in_region	Get all the barcodes that exist in a given region of the genome.
    count_bcs_list	Determine presence and quantity of given barcodes across a given region.
    plot_hmw		Generate a plot of the mapping locations of reads with each barcode.
    refine_bcs		Select barcodes that are shared between regions and not in other regions.
        """

#    select_bcs		Get barcodes shared by ALL region_in's and present in NONE of the region_out's.

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
	#parser.add_argument("-t","--vcf_test", metavar="FILE",
	#	dest="vcf_test",
	#	help="Test vcf file")
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
	parser.add_argument("-m","--mol_size",
		dest="mol_size", metavar="MOL_SIZE",
		help="Mean size of HMW molecules")
	parser.add_argument("-e", "--shared_bc_file", metavar="FILE",
		dest="shrd_file",
		help="File of shared bcs")

	return parser

def pipeline_from_parsed_args(args):
	if args.tool=="bedpe2window":
		pipeline = bedpe2window(bedpe=args.infile, window=args.window_size, out=args.outfile, mol_size=args.mol_size)
	if args.tool=="get_shared_bcs":
		pipeline = get_shared_bcs(sv=args.infile, bam=args.bam, out=args.outfile)
	if args.tool=="assign_sv_haps":
		pipeline = assign_sv_haps(sv=args.infile, window=args.window_size, vcf_control=args.vcf_control, vcf_test=args.vcf, out=args.outfile, bcs=args.bc_select)
	if args.tool=="count_bcs":
		pipeline = count_bcs(bam=args.bam, sv=args.infile, in_window=args.in_window, out_window=args.out_window, sv_name=args.sv_name, bcs=args.bc_select, out=args.outfile)
	if args.tool=="get_phased_basic":
		pipeline = get_phased_basic(vcf=args.vcf, out=args.outfile, chrom=args.chrom)
	if args.tool=="get_phase_blocks":
		pipeline = get_phase_blocks(infile_basic=args.infile, out=args.outfile)
	if args.tool=="get_phased_bcs":
		pipeline = get_phased_bcs(infile_basic=args.infile, ps=args.phase_block, out=args.outfile)
	#if args.tool=="select_bcs":
	#	pipeline = select_bcs(region_in=args.region_in, region_out=args.region_out, bam=args.bam, out=args.outfile)
	if args.tool=="count_bcs_list":
		pipeline = count_bcs_list(region=args.region_in, in_window=args.in_window, bam=args.bam, bcs=args.bcs, out=args.outfile)
	if args.tool=="get_bcs_in_region":
		pipeline = get_bcs_in_region(region=args.region_in,bam=args.bam, out=args.outfile)
	if args.tool=="plot_hmw":
		pipeline = plot_hmw(in_windows=args.infile, out=args.outfile)
	if args.tool=="extract_reads_interleaved":
		pipeline = extract_reads_interleaved(fqdir=args.fqdir, s_bcs=args.s_bcs, lanes=args.lanes, bcs=args.bcs, fq_outdir=args.outdir)
	if args.tool=="extract_reads_separate":
		pipeline = extract_reads_separate(bcs=args.bcs, fq_outdir=args.outdir, read1=args.read1, read2=args.read2, index1=args.index1)
	if args.tool=="refine_bcs":
		pipeline = refine_bcs(bed_in=args.infile, out=args.outfile, bam=args.bam, shrd_file = args.shrd_file)
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

##########################################################################################		
	if args.tool=="bedpe2window":
		#print "gemtools -T bedpe2window -i [LR_input.bedpe] -w [window_size] -o [out.bedpe]"
		if args.help:
			print """
Tool:	gemtools -T bedpe2window
Summary: Generate windows around SV breakpoints for SV analysis\n
Usage:   gemtools -T bedpe2window [OPTIONS] -i <LR_input.bedpe> -o <out.bedpe> -m <mol_size>
Input:
	-i  bedpe file of SV breakpoints (ex: sv_call.bedpe from Long Ranger)
	NOTE: must also supply either -m or -w (see Options below)
Output:
	-o  output file: bedpe file with windows around breakpoints
Options:
	-m  mean HMW size (Long Ranger estimates this); if supplied, gemtools will determine optimal window size for each SV breakpoint
	-w  size of window to generate around the breakpoints; if supplied, windows of this specified size will be generated for each SV breakpoint
			"""
			sys.exit(1)
		if not (args.infile or args.outfile):
			parser.error('Missing required input')
		if (not args.window_size and not args.mol_size):
			parser.error('Missing required input')
		if (args.window_size and args.mol_size):
			parser.error('Please only supply one of -m or -w (-m is recommended)')
		
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
	-i  output file from 'bedpe2window' tool
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
Output:
	-o  output file: List of breakpoints with phase id and number of barcodes supporting assignment to each haplotype
Options:
	-w  window size around breakpoint in which to look for phased SNVs; if not specified, will use windows already generated in input file (-i)
	-c  vcf file generated by Long Ranger for control sample (ex: normal sample); only specify if you want to use a different vcf to define phase blocks
			"""
			sys.exit(1)
		if not (args.infile or args.outfile or args.vcf or args.bcs):
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
		if not (args.infile or args.outfile or args.sv_name or args.bam):
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
Output:
	-o  output file: plot of barcode mapping locations in a given region
			"""
			sys.exit(1)		
		if not (args.infile or args.outfile):
			parser.error('Missing required input')
		
		if not os.path.isfile(args.infile):
			parser.error(str(args.infile) + " does not exist")

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
	if args.tool=="refine_bcs":
		#print "gemtools -T select_bcs [OPTIONS] -b [LR.bam] -f [region_in] -g [region_out] -o [out.bcs]"
		if args.help:
			print """
Tool:	gemtools -T refine_bcs
Summary: Get barcodes shared by ALL region_in's and present in NONE of the region_out's\n
Usage:   gemtools -T refine_bcs -i in.bed -b <LR.bam> -o <output.bcs.txt>
Input:
	-i  bed file of regions; order of columns must be: ['chrom','start','stop','name','status']; header line must be commented
	-b  bam file generated by Long Ranger
Output:
	-o  output file: barcode info summary for each event (specified by 'name')
Options:
	-e  name of file generated by 'get_shared_bcs'
			"""
			sys.exit(1)
		if not (args.bam or args.outfile or args.infile):
			parser.error('Missing required input')
		
		if not str(args.bam).endswith(".bam"):
			parser.error(str(args.bam) + " does not appear to be a bam file")

##########################################################################################
	pipeline = pipeline_from_parsed_args(args)
	runner = pipeline

if __name__ == '__main__':
	main()
