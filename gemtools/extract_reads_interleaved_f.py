import os, sys, argparse, __main__ as main
import csv
from itertools import izip_longest, islice
import gzip
import time


def extract_reads_interleaved(out_dir='fastq_bc',**kwargs):

	print out_dir

	if 'fqdir' in kwargs:
		fq_path = kwargs['fqdir']
	if 's_bcs' in kwargs:
		bcs = kwargs['s_bcs']
	if 'lanes' in kwargs:
		lanes = kwargs['lanes']	
	if 'bcs' in kwargs:
		bc_file = kwargs['bcs']
	if 'fq_outdir' in kwargs:
		out_dir = kwargs['fq_outdir']

	print out_dir

	if os.path.isdir(out_dir):
		print str(out_dir) + " already exists -- exiting"
		sys.exit()

	if not os.path.isdir(out_dir):
		os.makedirs(out_dir)
        cur_version = 1.0
        
        bc_list = bcs.split(",")

        lane_list = lanes.split(",")
        lane_list = ["lane-" + "0"*(3-len(l)) + l for l in lane_list] 
             
        file_list = os.listdir(fq_path) #list of things in directory

        for bc in bc_list: #for each barcode
                bc_files = [f for f in file_list if bc in f]

                for lane in lane_list: #for each lane
                        lane_files = [g for g in bc_files if lane in g]
                        ra_file = [r for r in lane_files if 'read-RA' in r][0]
                        i_file = [i for i in lane_files if 'read-I1' in i][0]

                        if not out_dir.endswith("/"):
                                out_dir = out_dir + "/"

                        if not fq_path.endswith("/"):
                                fq_path = fq_path + "/"

                        cmd_list = [fq_path + ra_file, fq_path + i_file, bc_file, out_dir + ra_file, out_dir + i_file]

                        extract_reads(cmd_list)


def extract_reads(args_fq):

        bcs = [] 
        out_file = gzip.open(args_fq[3],'w')  
        out_si_file = gzip.open(args_fq[4],'w')  

        with open(args_fq[2],'r') as f:
                for line in csv.reader(f,delimiter='\t'): 
                        bcs.append(line[0])

        bcs = set(bcs)

        n = 0
        i = 0

        with gzip.open(args_fq[0], 'r') as f, gzip.open(args_fq[1],'r') as ind:   
                cur_time = time.time()

                while True:
                        lines = list(islice(f,8)) #slice after index 8 of fastq
                        lines_index = list(islice(ind,4)) #slice after index 4 of si_fastq (file of indices)

                        if not lines:
                                break

                        n = n + 1
                        if (n % 1000000 == 0): 
                                print >>sys.stderr, "%d reads processed, %d records matched bcs in a %d second chunk" % (n, i, time.time() - cur_time)
                                cur_time = time.time()

                        if (lines[1][0:16] in bcs): #if barcodes match
                                i = i + 1

                                for line in lines:
                                        out_file.write(line)
                                for line in lines_index:
                                        out_si_file.write(line)

        


                       
                   
