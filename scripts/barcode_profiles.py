#!/usr/bin/env python
## ----------------------
# 10x Genomics Barcode profiling from bam file and Bed file
## ----------------------
# Srividya Ramakrishnan (srividya.ramki@gmail.com)
# John Hopkins University
## ----------------------
import sys, os 
import re , hashlib
import time,logging,traceback
import shutil
import pysam
import subprocess
import argparse
import pandas as pd
import numpy

def runProgram(logger=None,
	       progName=None,
	       argStr=None,
	       script_dir=None,
	       working_dir=None):
        """
        Convenience func to handle calling and monitoring output of external programs.
    
        :param progName: name of system program command
        :param argStr: string containing command line options for ``progName``
    
        :returns: subprocess.communicate object
        """
        # Ensure program is callable.
        if script_dir is not None:
                progPath= os.path.join(script_dir,progName)
        else:
		progPath = progName
        # Construct shell command
        cmdStr = "%s %s" % (progPath,argStr)
	logger.info("Processing Command : {0}".format(cmdStr))
	start = time.time()

        # Set up process obj
        process = subprocess.Popen(cmdStr,
                               shell=True,
                               stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE,
                               cwd=working_dir)
        # Get results
        result,stderr  = process.communicate()
        # keep this until your code is stable for easier debugging
        if logger is not None and result is not None and len(result) > 0:
            logger.info(result)
        if logger is not None and stderr is not None and len(stderr) > 0:
            logger.info(stderr)
	elasped=time.time()
        # Check returncode for success/failure
        if process.returncode != 0:
                raise RuntimeError('Return Code : {0} , result {1} , progName {2}'.format(process.returncode,result,progName))
        # Return result
        return { "result" : result , "stderr" :stderr, "exec_time" :' %.2fs' % (elasped - start) }

def getReadsinRegion(logger,bam_file,bed_file,output_file):
 """
  Get the Reads in the regions of bed file
 """ 
 cmd = "view -o {0} -hb -L {1} {2}".format(output_file,bed_file,bam_file)
 
 try:
 	res = runProgram(logger,"samtools",cmd,None,os.getcwd())
	logger.info("Extracting VCF regions from Alignment completed in {0}".format(res['exec_time']))
 except Exception, e:
	raise Exception("".join(traceback.format_exc()))


def buildSamIndex(logger,bam_file):
 """
  Build samtools index for the bam file
 """ 
 cmd = "index {0}".format(bam_file)
 print cmd
 try:
	dir_name=os.path.dirname(bam_file)
	res = runProgram(logger,"samtools",cmd,None,os.getcwd())
	logger.info("Building Samtools index for overlapping Alignments completed in {0}".format(res['exec_time']))
 except Exception, e:
	raise Exception("".join(traceback.format_exc()))

def buildcov(logging,bam_file,bed_file,output_bed):
 """
   Builds Barcode counts of the reads aligned to the genome
   input : bam_file and bed_file
   output : Returns bed file with Chr, start, end, barcode, no of occurences, strand
 """
 before = time.time()
 barcode_pattern = re.compile(r'BX:(A-Z)?:[ACGT]*-[0-9]')
 logger.info("WARNING: Dropping Reads Alignments that are not barcoded")
 output_bed = open(output_bed,"wb") 
 with open(bed_file) as regions:
    with pysam.Samfile(bam_file, "rb") as samfile:
	for line in regions:
		chr, start, end , info = line.rstrip('\n').rstrip('\r').split("\t")
 		b_dict = {}
		for read in samfile.fetch(chr,int(start),int(end)):
		    tags = dict(read.tags)
		    if 'BX' in tags.keys(): 
			if tags['BX'] in b_dict: b_dict[tags['BX']] += 1
			else: b_dict[tags['BX']] = 1 
	        for k, v in b_dict.items():
		    output_bed.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(chr,start,end,info,k,v)) 
	        
 regions.close()
 samfile.close()
 output_bed.close()
 after = time.time()
 logger.info('Computed Barcode profiles for the vcf regions completed in %.2fs' % (after - before))

def count_barcode_overlaps(logging,bed_file,tenkb_bed):
 """
   Count barcode overlap between the network paths and write it to a txt file 
 """
 before = time.time()
 ###
 tenkb_bed = open(output_bed,"wb")
 tenkb_bed.close()
 after = time.time()
 logging,info('Counting barcode overlaps in the vcf regions completed in %.2fs' % (after - before))

def count_barcode_overlaps(logging,bed_file,output_bed):
 """
   Count barcode overlap between the network paths and write it to a txt file 
 """
 before = time.time()
 ###
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
 output_bed.close()
 after = time.time()
 #logging.info('Counting barcode overlaps in the vcf regions completed in {0}'.format( str(after - before)))

if __name__ == "__main__":
    logger = logging.getLogger('barcode_profiling')
    streamHandler = logging.StreamHandler(sys.stdout)
    logger.setLevel(logging.INFO)
    formatter = logging.Formatter("%(asctime)s - %(filename)s - %(lineno)d - %(levelname)s - %(message)s")
    formatter.converter = time.gmtime
    streamHandler.setFormatter(formatter)
    logger.addHandler(streamHandler)
    parser = argparse.ArgumentParser()
    parser.add_argument("-bam",dest="bam_file",help="Barcoded Bam alignment file from 10x data")
    parser.add_argument("-bed",dest="bed_file",help="Bed file for the regions near variants")
    parser.add_argument("-o",dest="output",help="output dir name")
    args = parser.parse_args() 
    bam_file = args.bam_file
    bed_file = args.bed_file
    if not os.path.exists(args.output) : os.mkdir(args.output)
    output_file = os.path.join(args.output,"overlapping_reads.bam")
    getReadsinRegion(logger,args.bam_file,args.bed_file,output_file)
    buildSamIndex(logger,output_file)
    #if not os.path.exists(output_file):
    #	raise Exception("{0} Not Found".format(output_file))
    barcode_file = os.path.join(args.output,"reads_barcode_profile.txt")
    buildcov(logger,bam_file,bed_file,barcode_file)
    #tenkbfile="/seq/schatz/sramakri/NA24385_GRCh37/NA24385_GRCh37_chr20_shit_new_barcodes.bed" 
    overlap_file=os.path.join(args.output,"barcode_overlaps_between_regions.txt")
    count_barcode_overlaps(logger,barcode_file,overlap_file)
