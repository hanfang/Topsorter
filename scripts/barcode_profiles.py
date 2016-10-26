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
import pandas
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
        # Check returncode for success/failure
        if process.returncode != 0:
                raise RuntimeError('Return Code : {0} , result {1} , progName {2}'.format(process.returncode,result,progName))
        # Return result
        return { "result" : result , "stderr" :stderr }

def getReadsinRegion(logger,bam_file,bed_file,output_file):
 """
  Get the Reads in the regions of bed file
 """ 
 cmd = "view -o {0} -L {1} {2}".format(output_file,bam_file,bed_file)
 try:
 	runProgram(logger,"samtools",cmd,None,os.getcwd())
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
 print 'completed in %.2fs' % (after - before)
 pass

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
    parser.add_argument("-o",dest="output",help="output file name")
    args = parser.parse_args() 
    bam_file = args.bam_file
    bed_file = args.bed_file
    #if not os.path.exists(args.output_dir) : os.mkdir(args.output_dir)
    #output_file = os.path.join(args.output_dir,"overlapping_reads.sam")
    #getReadsinRegion(logger,args.bam_file,args.bed_file,output_file)
    #if not os.path.exists(output_file):
    #	raise Exception("{0} Not Found".format(output_file))
    buildcov(logger,bam_file,bed_file,args.output) 
