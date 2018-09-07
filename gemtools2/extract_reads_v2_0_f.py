#!/usr/bin/env python #is this necessary within package?

### LOAD THE NECESSARY PACKAGES ###
##This is for version 2.0(?)

import os, sys, argparse, __main__ as main
import pandas as pd

import csv
#import pysam
#import distance
from itertools import izip_longest, islice
import gzip
#import sys
import time

#def extract_reads_v2_0(fq_path, bcs, lanes, bc_file, out_dir): #fastq, si_fastq
def extract_reads_v2_0(out_dir='fastq_bc',**kwargs):

	if 'fq' in kwargs:
		fq_path = kwargs['fq']
	if 'sample_bcs' in kwargs:
		bcs = kwargs['sample_bcs']
	if 'l' in kwargs:
		lanes = kwargs['l']	
	if 'mol_bcs' in kwargs:
		bc_file = kwargs['mol_bcs']
	if 'out' in kwargs:
		out_dir = kwargs['out']

	if not os.path.isdir(out_dir):
		os.makedirs(out_dir)
		        
        #"/mnt/ix1/Seq_Runs/20160122_CC3_0317/Analysis/fastq"
        #"ACGACGCT,CGCCATTC,GTAGTCAG,TATTGAGA"
        #"1,5"
        #"bcs_metr.txt" "/mnt/ix2/avitko/170621_SV_phasing/A02_bcl_to_fastq/bcs_metr.txt"
        #"/mnt/ix2/avitko/170621_SV_phasing/A02_bcl_to_fastq/out_metr"
        cur_version = 1.0
        
        bc_list = bcs.split(",")

        lane_list = lanes.split(",")
        lane_list = ["lane-" + "0"*(3-len(l)) + l for l in lane_list] 
             
        file_list = os.listdir(fq_path) #list of things in directory

        counter = 0

        for bc in bc_list: #for each barcode
                bc_files = [f for f in file_list if bc in f]

                for lane in lane_list: #for each lane
                        counter = counter + 1
                        lane_files = [g for g in bc_files if lane in g]
                        ra_file = [r for r in lane_files if 'read-RA' in r][0]
                        i_file = [i for i in lane_files if 'read-I1' in i][0]

                        if not out_dir.endswith("/"):
                                out_dir = out_dir + "/"

                        if not fq_path.endswith("/"):
                                fq_path = fq_path + "/"

                        cmd_list = [fq_path + ra_file, fq_path + i_file, bc_file, out_dir + ra_file, out_dir + i_file]


                        extract_reads(cmd_list)

#if __name__ == '__main__':
                #args = parse_args()
                #if(args.usage):
                        #usage()
                #if(not fq_path or not bcs or not lanes or not bc_file or not out_dir): #or not args.s_in):
                        #print os.path.basename(main.__file__) + " missing a required input\n"
                        #usage()
                        #sys.exit(1)

def extract_reads(args_fq): #from (args_fq) ##fastq, si_fastq  #this used to be main

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

                        if (lines[1][0:16] in bcs): #if barcodes match?
                                i = i + 1

                                for line in lines:
                                        out_file.write(line)
                                for line in lines_index:
                                        out_si_file.write(line)

        


                       
                   
