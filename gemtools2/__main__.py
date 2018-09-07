#!/usr/bin/env python


"""
gemtools2 version %version
Copyright (C) 2018 Stephanie Greer <sgreer2@stanford.edu>
gemtools2 is a suite of tools to work with linked read sequencing data.
Usage:
    gemtools2 -T TOOL [options] [-o output]
"""


import sys
from optparse import OptionParser, OptionGroup

from gemtools2 import __version__
from gemtools2.bedpe2window_f import bedpe2window
from gemtools2.get_shared_bcs_f import get_shared_bcs
from gemtools2.assign_sv_haps_f import assign_sv_haps
from gemtools2.count_bcs_f import count_bcs
from gemtools2.get_phased_basic_f import get_phased_basic
from gemtools2.get_phased_basic_chr_f import get_phased_basic_chr
from gemtools2.get_phase_blocks_f import get_phase_blocks
from gemtools2.get_bcs_in_region_f import get_bcs_in_region
from gemtools2.get_phased_bcs_f import get_phased_bcs
from gemtools2.select_bcs_f import select_bcs
from gemtools2.count_bcs_list_f import count_bcs_list
from gemtools2.extract_reads_v2_0_f import extract_reads_v2_0
from gemtools2.extract_reads_v2_1_f import extract_reads_v2_1
from gemtools2.get_hmw_summary_f import get_hmw_summary


class GemtoolsOptionParser(OptionParser):
	def get_usage(self):
		return self.usage.lstrip().replace('%version', __version__)

class CommandlineError(Exception):
	pass

def get_option_parser():
	parser = GemtoolsOptionParser(usage=__doc__, version=__version__)
	
	parser.add_option("-T", "--tool", default=None, dest="tool",
		help="Name of tool to use "
		"gemtools2 -T bedpe2window -i INPUT_FILE -w WINDOW_SIZE -o OUTPUT_FILE"
		"gemtools2 -T get_shared_bcs -i INPUT_FILE -b BAM_FILE -o OUTPUT_FILE"
		"gemtools2 -T assign_sv_haps -i INPUT_FILE -c VCF_CONTROL_FILE -t VCF_TEST_FILE -o OUTPUT_FILE"
		"gemtools2 -T count_bcs -i INPUT_FILE -x SMALL_WINDOW -y BIG_WINDOW -s SV_NAME -q (all|shared) -o OUTPUT_FILE
		"gemtools2 -T get_phased_basic -v VCF_FILE -o OUTPUT_FILE"
		"gemtools2 -T get_phased_basic_chr -v VCF_FILE -n CHROM -o OUTPUT_FILE"
		"gemtools2 -T get_phase_blocks -i INPUT_FILE -o OUTPUT_FILE"
		"gemtools2 -T get_bcs_in_region -b BAM_FILE -f REGION -o OUTPUT_FILE"
		"gemtools2 -T get_phased_bcs -i INPUT_FILE -p PHASE_ID -o OUTPUT_FILE"
		"gemtools2 -T select_bcs -b BAM_FILE -f REGION_IN -g REGION_OUT -o OUTPUT_FILE"
		"gemtools2 -T count_bcs_list -b BAM_FILE -f REGION -x WINDOW -l BARCODES_FILE -o OUTPUT_FILE")
	parser.add_option("-o", "--output", metavar="FILE",
		dest="outfile",
		help="Name of output file")

	group = OptionGroup(parser, "Input files")
	group.add_option("-i", "--input", metavar="FILE",
		dest="infile",
		help="Name of input file ")
	group.add_option("-b","--bam", metavar="FILE",
		dest="bam",
		help="bam file")
	group.add_option("-v","--vcf", metavar="FILE",
		dest="vcf",
		help="vcf file")
	group.add_option("-c","--vcf_control", metavar="FILE",
		dest="vcf_control",
		help="Control vcf file")
	group.add_option("-t","--vcf_test", metavar="FILE",
		dest="vcf_test",
		help="Test vcf file")
	parser.add_option_group(group)	
	
	group = OptionGroup(parser, "Windows")
	group.add_option("-w", "--window", type=int, default=50000,
		dest="window_size", metavar="WINDOW",
		help="BEDPE file of SVs (longranger output)")
	group.add_option("-x","--in_window", type=int, default=1000,
		dest="in_window",metavar="WINDOW",
		help="Size of small windows")
	group.add_option("-y","--out_window", type=int, default=50000,
		dest="out_window",metavar="WINDOW",
		help="Size of large window")
	parser.add_option_group(group)

	group = OptionGroup(parser, "Regions")
	group.add_option("-f","--region_in",
		dest="region_in", metavar='REGION',
		help="In regions")
	group.add_option("-g","--region_out",
		dest="region_out", metavar='REGION',
		help="Out regions")		
	parser.add_option_group(group)

	group = OptionGroup(parser, "Barcodes")
	group.add_option("-l","--bc_list", metavar="FILE",
		dest="bcs",
		help="File with list of barcodes")			
	group.add_option("-q","--bc_select",metavar='(all|shared)',
		dest="bc_select",choices=('all', 'shared'),
		help="BCs to consider: all bcs or shared bcs")
	parser.add_option_group(group)

	group = OptionGroup(parser, "Specifics")
	group.add_option("-n","--chrom", metavar="CHR",
		dest="chrom",
		help="Chromosome number; ex: 'chr22','22'")
	group.add_option("-p","--phase_block",
		dest="phase_block", metavar="PHASE_ID",
		help="Phase block id (from vcf)")
	group.add_option("-s","--sv_name",
		dest="sv_name", metavar="SV",
		help="Name of SV; ex: 'call_144', '144'")
	parser.add_option_group(group)

	return parser

def pipeline_from_parsed_args(options):
	if options.tool=="bedpe2window":
		pipeline = bedpe2window(bedpe=options.infile, window=options.window_size, out=options.outfile)
	if options.tool=="get_shared_bcs":
		pipeline = get_shared_bcs(sv=options.infile, bam=options.bam, out=options.outfile)
	if options.tool=="assign_sv_haps":
		pipeline = assign_sv_haps(sv=options.infile, vcf_control=options.vcf_control, vcf_test=options.vcf_test, out=options.outfile)
	if options.tool=="count_bcs":
		pipeline = count_bcs(sv=options.infile, in_window=options.in_window, out_window=options.out_window, sv_name=options.sv_name, bcs=options.bc_select, out=options.outfile)
	if options.tool=="get_phased_basic":
		pipeline = get_phased_basic(vcf=options.vcf, out=options.outfile)
	if options.tool=="get_phased_basic_chr":
		pipeline = get_phased_basic_chr(vcf=options.vcf, chr=options.chrom, out=options.outfile)
	if options.tool=="get_phase_blocks":
		pipeline = get_phase_blocks(infile_basic=options.infile, out=options.outfile)
	if options.tool=="get_phased_bcs":
		pipeline = get_phased_bcs(infile_basic=options.infile, ps=options.phase_block, out=options.outfile)
	if options.tool=="select_bcs":
		pipeline = select_bcs(region_in=options.region_in, region_out=options.region_out, bam=options.bam, out=options.outfile)
	if options.tool=="count_bcs_list":
		pipeline = count_bcs_list(region=options.region_in, in_window=options.in_window, bam=options.bam, bcs=options.bcs, out=options.outfile)
	if options.tool=="get_bcs_in_region":
		pipeline = get_bcs_in_region(region=options.region_in,bam=options.bam, out=options.outfile)
	return pipeline


def main(cmdlineargs=None):
	parser = get_option_parser()
	if cmdlineargs is None:
		cmdlineargs = sys.argv[1:]
	options, args = parser.parse_args(args=cmdlineargs)
	pipeline = pipeline_from_parsed_args(options)
	runner = pipeline

if __name__ == '__main__':
	main()