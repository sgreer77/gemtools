#!/usr/bin/env python

import sys
from optparse import OptionParser, OptionGroup

from gemtools2.bedpe2window_f import bedpe2window

class GemtoolsOptionParser(OptionParser):
	def get_usage(self):
		return self.usage.lstrip().replace('%version', '0.1')

class CommandlineError(Exception):
	pass

def get_option_parser():
	parser = GemtoolsOptionParser(usage=__doc__, version='0.1')
	
	parser.add_option("-t", "--tool", default=None,
		help="Name of tool to use", dest="tool")
	parser.add_option("-o", "--output", metavar="FILE",
		dest="outpre",
		help="Name of output file")
	
	group = OptionGroup(parser, "Window breakpoints",
		description="Parameters -b, -w specify file and windows")
	group.add_option("-b", "--bedpe", metavar="FILE",
		dest="sv_input",
		help="BEDPE file of SVs (longranger output)")
	group.add_option("-w", "--window", type=int, default=50000,
		dest="window_size",
		help="BEDPE file of SVs (longranger output)")
	parser.add_option_group(group)
	
	return parser

def pipeline_from_parsed_args(options):
	if options.tool=="bedpe2window":
		pipeline = bedpe2window()
	return pipeline


def main(cmdlineargs=None):
	parser = get_option_parser()
	print cmdlineargs
	if cmdlineargs is None:
		cmdlineargs = sys.argv[1:]
		print cmdlineargs
	options, args = parser.parse_args(args=cmdlineargs)
	print options
	print args
	pipeline = pipeline_from_parsed_args(options)
	runner = pipeline

if __name__ == '__main__':
	main()