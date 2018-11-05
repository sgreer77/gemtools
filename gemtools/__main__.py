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
#from optparse import OptionParser, OptionGroup
import argparse

from gemtools import __version__
from gemtools.bedpe2window_f import bedpe2window
from gemtools.get_shared_bcs_f import get_shared_bcs
from gemtools.assign_sv_haps_f import assign_sv_haps
from gemtools.count_bcs_f import count_bcs
from gemtools.get_phased_basic_f import get_phased_basic
from gemtools.get_phase_blocks_f import get_phase_blocks
from gemtools.get_bcs_in_region_f import get_bcs_in_region
from gemtools.get_phased_bcs_f import get_phased_bcs
from gemtools.select_bcs_f import select_bcs
from gemtools.count_bcs_list_f import count_bcs_list
from gemtools.plot_hmw_f import plot_hmw
from gemtools.extract_reads_interleaved_f import extract_reads_interleaved
from gemtools.extract_reads_separate_f import extract_reads_separate

from gemtools.get_hmw_summary_f import get_hmw_summary


#class GemtoolsOptionParser(OptionParser):
#	def get_usage(self):
#		return self.usage.lstrip().replace('%version', __version__)
#		print("test1")
#	def error(self, msg):
#		print('Run "gemtools --help" to see command-line args.')
#		print("test2")
#		self.exit(2, "%s: error: %s\n" % (self.get_prog_name(), msg))

class CommandLineError(Exception):
	pass

def msg(name=None):                                                            
    return '''\tgemtools -T <tool> [options]
        '''

help_msg = """\tgemtools -T <tool> [options]\n 
The gemtools sub-commands include:\n 
[ Phase blocks ] 
    get_phased_basic   Obtain phasing information for all SNVs in the vcf file 
    get_phase_blocks   Summarize phase blocks -- coordinates, size, SNVs per phase block etc.\n
[ SV analysis tools ]
    bedpe2window	Generate windows around SV breakpoints for SV analysis
    get_shared_bcs	Determine barcodes shared between SV breakpoints
    assign_sv_haps	Assign SV barcodes to existing haplotypes (SNVs)
    count_bcs		Determine presence and quantity of given barcodes across a given region surrounding the SV breakpoints
    plot_hmw		Generate a plot of the mapping locations of reads with each barcode (SAME AS ABOVE) \n
[ Subset reads by barcode ]
    extract_reads_separate	Obtain reads with particular barcodes from Long Ranger fastq files (R1,R2,I1)
    extract_reads_interleaved	Obtain reads with particular barcodes from Long Ranger fastq files (RA,I1,I2)\n
[ General tools ]
    get_phased_bcs	For a particular phase block, return the haplotype 1 and haplotype 2 barcodes
    select_bcs		Get barcodes shared by ALL region_in's and present in NONE of the region_out's
    get_bcs_in_region	Get all the barcodes that exist in a given region of the genome
    count_bcs_list	Determine presence and quantity of given barcodes across a given region
    plot_hmw		Generate a plot of the mapping locations of reads with each barcode
        """

def get_option_parser():
	parser = argparse.ArgumentParser(description='Process some integers.', add_help=False, usage=msg())
	
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
	parser.add_argument("-t","--vcf_test", metavar="FILE",
		dest="vcf_test",
		help="Test vcf file")
	#parser.add_argument_group(group)	
	
	#group = OptionGroup(parser, "Windows")
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
	#parser.add_argument_group(group)

	#group = OptionGroup(parser, "Regions")
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
	#parser.add_argument_group(group)

	#group = OptionGroup(parser, "Barcodes")
	parser.add_argument("-l","--bc_list", metavar="FILE",
		dest="bcs",
		help="File with list of barcodes")			
	parser.add_argument("-q","--bc_select",metavar='(all|shared)',
		dest="bc_select",choices=('all', 'shared'), default="shared",
		help="BCs to consider: all bcs or shared bcs               "
			"default: shared")
	#parser.add_argument_group(group)

	#group = OptionGroup(parser, "Specifics")
	parser.add_argument("-n","--chrom", metavar="CHR",
		dest="chrom", default="None",
		help="Chromosome number; ex: 'chr22','22'")
	parser.add_argument("-p","--phase_block",
		dest="phase_block", metavar="PHASE_ID",
		help="Phase block id (from vcf)")
	parser.add_argument("-s","--sv_name",
		dest="sv_name", metavar="SV",
		help="Name of SV; ex: 'call_144', '144'")
	#parser.add_argument_group(group)

	#group = OptionGroup(parser, "Fastq")
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
	parser.add_argument("-m","--read1",
		dest="read1", metavar="READ1",
		help="read1 fastq file")
	parser.add_argument("-u","--read2",
		dest="read2", metavar="READ2",
		help="read2 fastq file")
	parser.add_argument("-r","--index1",
		dest="index1", metavar="INDEX1",
		help="index1 fastq file")
	#parser.add_argument_group(group)

	return parser

def pipeline_from_parsed_args(options):
	if args.tool=="bedpe2window":
		pipeline = bedpe2window(bedpe=args.infile, window=args.window_size, out=args.outfile)
	if args.tool=="get_shared_bcs":
		pipeline = get_shared_bcs(sv=args.infile, bam=args.bam, out=args.outfile)
	if args.tool=="assign_sv_haps":
		pipeline = assign_sv_haps(sv=args.infile, window=args.window_size, vcf_control=args.vcf_control, vcf_test=args.vcf_test, out=args.outfile)
	if args.tool=="count_bcs":
		pipeline = count_bcs(bam=args.bam, sv=args.infile, in_window=args.in_window, out_window=args.out_window, sv_name=args.sv_name, bcs=args.bc_select, out=args.outfile)
	if args.tool=="get_phased_basic":
		pipeline = get_phased_basic(vcf=args.vcf, out=args.outfile, chrom=args.chrom)
	if args.tool=="get_phase_blocks":
		pipeline = get_phase_blocks(infile_basic=args.infile, out=args.outfile)
	if args.tool=="get_phased_bcs":
		pipeline = get_phased_bcs(infile_basic=args.infile, ps=args.phase_block, out=args.outfile)
	if args.tool=="select_bcs":
		pipeline = select_bcs(region_in=args.region_in, region_out=args.region_out, bam=args.bam, out=args.outfile)
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
	return pipeline

def main(cmdlineargs=None):
	parser = get_option_parser()
	
	if cmdlineargs is None:
		cmdlineargs = sys.argv[1:]
		
	#options, args = parser.parse_args(args=cmdlineargs)
	args = parser.parse_args()
	
	if args.help:
		print help_msg
		sys.exit()
	
	if not args.tool:
		parser.error("Must provide a tool to run with -T")
		
	if args.tool=="bedpe2window":
		print "gemtools -T bedpe2window -i [LR_input.bedpe] -w [window_size] -o [out.bedpe]"
		if not (args.infile or args.outfile):
			parser.error('Missing required input')
		
		if not os.path.isfile(args.infile):
			parser.error(str(args.infile) + " does not exist")
		
		#if str(args.window_size).isdigit() and int(args.window_size)>0:
		#	print "window_size: " + str(args.window_size)
		#else:
		#	parser.error(str(args.window_size) + " must be an integer >0")
	
	if args.tool=="get_shared_bcs":
		print "gemtools -T get_shared_bcs -i [out.bedpe] -b [LR_bam_file] -o [out.shared]"
		if not (args.infile or args.outfile or args.bam):
			parser.error('Missing required input')

		if not os.path.isfile(args.infile):
			parser.error(str(args.infile) + " does not exist")
		if os.path.isfile(args.bam):
			print "input file: " + str(args.bam)
		else:
			parser.error(str(args.bam) + " does not exist")
		if str(args.bam).endswith(".bam"):
			print "Bam file: " + str(args.bam)
		else:
			parser.error(str(args.bam) + " does not appear to be a bam file")

	if args.tool=="assign_sv_haps":
		print "gemtools -T assign_sv_haps -i [out.shared] -c [LR_control.vcf.gz] -t [LR_test.vcf.gz] -o [out.haps]"
		if not (args.infile or args.outfile or args.vcf_control or args.vcf_test):
			parser.error('Missing required input')

		if not str(args.vcf_control).endswith(".vcf.gz"):
			parser.error(str(args.vcf_control) + " does not appear to be a gzipped vcf file")		
		if not str(args.vcf_test).endswith(".vcf.gz"):
			parser.error(str(args.vcf_test) + " does not appear to be a gzipped vcf file")

	if args.tool=="count_bcs":
		print "gemtools -T count_bcs -i [out.shared] -b [LR.bam] -x [in_window] -y [out_window] -s [sv_name] -q [all|shared] -o [out.bc_count]"
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

	if args.tool=="get_phased_basic":
		print "gemtools -T get_phased_basic -v [LR.vcf.gz] -o [output.phased_basic] -n [chr_num]"
		if not (args.outfile or args.vcf):
			parser.error('Missing required input')

		if not str(args.vcf).endswith(".vcf.gz"):
			parser.error(str(args.vcf) + " does not appear to be a gzipped vcf file")
	
	if args.tool=="get_phase_blocks":
		print "gemtools -T get_phase_blocks -i [out.phased_basic] -o [out.phase_blocks]"
		if not (args.infile or args.outfile):
			parser.error('Missing required input')

		if not os.path.isfile(args.infile):
			parser.error(str(args.infile) + " does not exist")

	if args.tool=="get_bcs_in_region":
		print "gemtools -T get_bcs_in_region -b [LR.bam] -f [region_in] -o [out.bcs]"
		if not (args.outfile or args.bam or args.region_in):
			parser.error('Missing required input')
		
		if not str(args.bam).endswith(".bam"):
			parser.error(str(args.bam) + " does not appear to be a bam file")
	
	if args.tool=="get_phased_bcs":
		print "gemtools -T get_phased_bcs -i [out.phased_basic] -p [phase_block_id] -o [out.phased_bcs]"
		if not (args.infile or args.outfile or args.phase_block):
			parser.error('Missing required input')

		if not os.path.isfile(args.infile):
			parser.error(str(args.infile) + " does not exist")

	if args.tool=="select_bcs":
		print "gemtools -T select_bcs -b [LR.bam] -f [region_in] -g [region_out] -o [out.bcs]"
		if not (args.bam or args.outfile or args.region_in or args.region_out):
			parser.error('Missing required input')
			
		if not str(args.bam).endswith(".bam"):
			parser.error(str(args.bam) + " does not appear to be a bam file")

	if args.tool=="count_bcs_list":
		print "gemtools -T count_bcs_list -b [LR.bam] -f [region_in] -x [in_window] -l [bc_list] -o [out.bc_count]"
		if not (args.bam or args.region_in or args.outfile or args.bcs):
			parser.error('Missing required input')

		if not str(args.bam).endswith(".bam"):
			parser.error(str(args.bam) + " does not appear to be a bam file")
		if not (str(args.in_window).isdigit() and int(args.in_window)>0):
			parser.error(str(args.in_window) + " must be an integer >0")
		if not os.path.isfile(args.bcs):
			parser.error(str(args.bcs) + " does not exist")

	if args.tool=="plot_hmw":
		print "gemtools -T plot_hmw -i [out.bc_count] -o [out.pdf]"
		if not (args.infile or args.outfile):
			parser.error('Missing required input')
		
		if not os.path.isfile(args.infile):
			parser.error(str(args.infile) + " does not exist")

	if args.tool=="extract_reads_interleaved":
		print "gemtools -T extract_reads_interleaved -l [bc_list] -z [fastq_output_dir] -d [LR_fastq_dir] -j [sample_barcodes] -k [sample_lanes]"
		if not (args.fqdir or args.s_bcs or args.lanes or args.bcs or args.outdir):
			parser.error('Missing required input')
		
		if not os.path.isfile(args.bcs):
			parser.error(str(args.bcs) + " does not exist")
		if os.path.isdir(args.outdir):
			parser.error(str(args.outdir) + " already exists")
	
	if args.tool=="extract_reads_separate":
		print "gemtools -T extract_reads_separate -l [bc_list] -z [fastq_output_dir] --read1 [LR_R1.fastq.gz] --read2 [LR_R2.fastq.gz] --index1 [LR_I1.fastq.gz]"
		if not (args.bcs or args.read1 or args.read2 or args.index1 or args.outdir):
			parser.error('Missing required input')
			
		if not os.path.isfile(args.bcs):
			parser.error(str(args.bcs) + " does not exist")
		
	
	pipeline = pipeline_from_parsed_args(options)
	runner = pipeline

if __name__ == '__main__':
	main()
