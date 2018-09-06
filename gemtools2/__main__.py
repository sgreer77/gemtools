#!/usr/bin/env python

import sys
from optparse import OptionParser, OptionGroup

from gemtools.bedpe2window_f import bedpe2window

class GemtoolsOptionParser(OptionParser):
	def get_usage(self):
		return self.usage.lstrip().replace('%version', __version__)

class CommandlineError(Exception):
	pass

def get_option_parser():
	parser = GemtoolsOptionParser(usage=__doc__, version=__version__)
	
	parser.add_option("-t", "--tool", default=None,
		help="Name of tool to use")
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

def main(cmdlineargs=None):
	if cmdlineargs is None:
		cmdlineargs = sys.argv[1:]
	options, args = parser.parse_args(args=cmdlineargs)
	print options
	print args

if __name__ == '__main__':
	main()