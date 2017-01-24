#!/usr/bin/env python
## ----------------------
# Computes barcode overlaps between the before->left,before->right, before->after , left->right , left->after , right->after
## ----------------------
# Srividya Ramakrishnan (srividya.ramki@gmail.com)
# John Hopkins University
## ----------------------
import sys, os 
import re , hashlib
import time,logging,traceback,datetime
import shutil
import pysam
import subprocess
import argparse
import numpy
import pandas as pd

def count_barcode_overlaps(logging,bed_file,output_bed):
 """
   Count barcode overlap between the network paths and write it to a txt file 
 """
 before = time.time()
 ###
 barcode_pattern = re.compile(r'BX:(A-Z)?:[ACGT]*-[0-9]')
 output_bed = open(output_bed,"wb")
 data = pd.read_table(bed_file, header = None)
 #print data.head(10)
 v_infos = set(data[3])
 u_infos= set([ v.split(",")[0] for v in list(v_infos)])
 #print len(v_infos)
 #print len(u_infos)
 #print u_infos
 output_bed.write("Variant_Info\tbefore->left\tbefore->right\tbefore->after\tleft->right\tleft->after\tright->after\n")
 for u in u_infos:
	#print u+",before"
	#print u+',after'
 	before = data.loc[data[3] == u+",before"]
	after = data.loc[data[3] == u+',after']
	left =  data.loc[data[3] == u+',left']
	right = data.loc[data[3] == u+',right']
	b_to_a = set(before[4]) & set(after[4]) 
	b_to_l = set(before[4]) & set(left[4]) 
	b_to_r = set(before[4]) & set(right[4]) 
	l_to_r = set(left[4]) & set(right[4]) 
	l_to_a = set(left[4]) & set(after[4]) 
	r_to_a = set(right[4]) & set(after[4]) 
	output_bed.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\n".format(u,str(len(b_to_a)),str(len(b_to_l)),str(len(b_to_r)),str(len(l_to_r)),str(len(l_to_a)),str(len(r_to_a))))
 output_bed.close
 after = time.time()
 #print 'Counting barcode overlaps in the vcf regions completed in {0} seconds '.format((after - before))
 pass

if __name__ == "__main__":
    logger = logging.getLogger('compute_barcode_overlaps')
    streamHandler = logging.StreamHandler(sys.stdout)
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s - %(filename)s - %(lineno)d - %(levelname)s - %(message)s")
    formatter.converter = time.gmtime
    streamHandler.setFormatter(formatter)
    logger.addHandler(streamHandler)
    parser = argparse.ArgumentParser()
    parser.add_argument("-bed",dest="bed_file",help="Bed file for the regions near variants")
    parser.add_argument("-o",dest="output",help="output file name")
    args = parser.parse_args() 
    count_barcode_overlaps(logging,args.bed_file,args.output)
