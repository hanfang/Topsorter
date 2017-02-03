#!/usr/bin/env python
## ----------------------
# 10x Genomics Barcode profiling from bam file and Bed file
## ----------------------
# Srividya Ramakrishnan (srividya.ramki@gmail.com)
# John Hopkins University
## ----------------------
from __future__ import print_function
import sys, os 
import re , hashlib
import time,logging,traceback
import shutil
import pysam
import subprocess
import argparse
import pandas as pd
import numpy
from datetime import datetime
#import multiprocessing as mp

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
	start = datetime.now()

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
	elasped=datetime.now()
        # Check returncode for success/failure
        if process.returncode != 0:
                raise RuntimeError('Return Code : {0} , result {1} , progName {2}'.format(process.returncode,result,progName))
        # Return result
        return { "result" : result , "stderr" :stderr, "exec_time" :'{0}'.format(str(elasped - start)) }


class BarcodeProfiler:
     '''
     Module Name:
     BarcodeProfiler
     Module Description:
     '''
        
     def __init__(self,bam_file,bed_file,out_dir):
        self.bam_file=bam_file
	self.bed_file=bed_file
	###TODO : Check the bed file format and throw error if needed"
        self.out_dir=out_dir
	self.obam_file=out_dir+"/overlapping_reads.bam"
	self.barcode_bed=out_dir+"/reads_barcode_profile.txt"
	self.overlap_file=out_dir+"/barcode_overlaps_between_regions.txt"
	# logging
        self.__LOGGER = logging.getLogger('BarcodeProfiler')
        self.__LOGGER.setLevel(logging.INFO)
        streamHandler = logging.StreamHandler(sys.stdout)
        formatter = logging.Formatter("%(asctime)s - %(filename)s - %(lineno)d - %(levelname)s - %(message)s")
        formatter.converter = time.gmtime
        streamHandler.setFormatter(formatter)
        self.__LOGGER.addHandler(streamHandler)
        self.__LOGGER.info("Logger was set")

     def getReadsinRegion(self):
        """
          Get the Reads in the regions of bed file
        """ 
	
        cmd = "view -o {0} -hb -L {1} {2}".format(self.obam_file,self.bed_file,self.bam_file)
 
        try:
             res = runProgram(self.__LOGGER,"samtools",cmd,None,os.getcwd())
             self.__LOGGER.info("Extracting VCF regions from Alignment completed in {0}".format(res['exec_time']))
        except Exception, e:
             raise Exception("".join(traceback.format_exc()))


     def buildSamIndex(self):
        """
           Build samtools index for the bam file
        """ 
        cmd = "index {0}".format(self.obam_file)
        print(cmd)
        try:
	    dir_name=os.path.dirname(self.obam_file)
	    res = runProgram(self.__LOGGER,"samtools",cmd,None,os.getcwd())
	    self.__LOGGER.info("Building Samtools index for overlapping Alignments completed in {0}".format(res['exec_time']))
        except Exception, e:
	    raise Exception("".join(traceback.format_exc()))

     def buildcov(self):
       """
          Builds Barcode counts of the reads aligned to the genome
          input : bam_file and bed_file
          output : Returns bed file with Chr, start, end, barcode, no of occurences, strand
       """
       before = datetime.now()
       barcode_pattern = re.compile(r'BX:(A-Z)?:[ACGT]*-[0-9]')
       self.__LOGGER.info("WARNING: Dropping Reads Alignments that are not barcoded")
       output_bed = open(self.barcode_bed,"wb") 
       with open(self.bed_file) as regions:
             with pysam.Samfile(self.obam_file, "rb") as samfile:
	         for line in regions:
		    chr, start, end , info = line.rstrip('\n').rstrip('\r').split("\t")
 		    b_dict = {}
		    try:
                        readsinregion=samfile.fetch(chr,int(start),int(end))
                        for read in readsinregion:
		            tags = dict(read.tags)
		            if 'BX' in tags.keys(): 
			        if tags['BX'] in b_dict: b_dict[tags['BX']] += 1
		 	        else: b_dict[tags['BX']] = 1
                    except:
                        self.__LOGGER.info("Skipping region : {0}:{1}:{2}".format(chr,start,end))
                        continue 
	            for k, v in b_dict.items():
                        output_bed.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\n".format(chr,start,end,info,k,v)) 
	        
       regions.close()
       samfile.close()
       output_bed.close()
       after = datetime.now()
       self.__LOGGER.info('Computed Barcode profiles for the vcf regions completed in {0}'.format(str(after - before)))


     def count_barcode_overlaps(self):
        """
           Count barcode overlap between the network paths and write it to a txt file 
        """
        before = datetime.now()
        ###
        output_bed = open(self.overlap_file,"wb")
        data = pd.read_table(self.barcode_bed, header = None)
        #print data.head(10)
        v_infos = set(data[3])
        #u_infos= set([ v.split(",")[-3]+","+v.split(",")[-2]+","+v.split(",")[-1] for v in list(v_infos)])
        u_infos= set([ ",".join(v.split(",")[0:3]) for v in list(v_infos)])
        #print len(v_infos)
        #print len(u_infos)
	output_bed.write("#chr\tstart\tend\tvariant_id\tbefore->left\tbefore->right\tbefore->after\tleft->right\tleft->after\tright->after\tunique_BA\tunique_LR\n")
        for u in u_infos:
	    v_id,start,end=u.split(",")
            before = data.loc[data[3] == u+",before"]
            after = data.loc[data[3] == u+',after']
            left =  data.loc[data[3] == u+',left']
            right = data.loc[data[3] == u+',right']
	    chr =  "".join(set(before[0]) | set(left[0]) | set(right[0]) |  set(after[0]))
            b_to_a = set(before[4]) & set(after[4])
            b_to_l = set(before[4]) & set(left[4])
            b_to_r = set(before[4]) & set(right[4])
            l_to_r = set(left[4]) & set(right[4])
            l_to_a = set(left[4]) & set(after[4])
            r_to_a = set(right[4]) & set(after[4])
            unique_BA = set(before[4])| set(after[4]) - set(after[4]) &  set(before[4])
            unique_LR = set(left[4]) |  set(right[4]) -  set(right[4]) & set(left[4])
            output_bed.write("{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}\t{8}\t{9}\t{10}\t{11}\n".format(chr,start,end,v_id,str(len(b_to_a)),str(len(b_to_l)),str(len(b_to_r)),str(len(l_to_r)),str(len(l_to_a)),str(len(r_to_a)),str(len(unique_BA)),str(len(unique_LR))))
        output_bed.close()
        after = datetime.now()
        #self.__LOGGER.info('Counting barcode overlaps in the vcf regions completed in {0}'.format(str(after - before)))

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
     worker = BarcodeProfiler(args.bam_file,args.bed_file,args.output)
     time = str(datetime.now())
     sys.stderr.write("[status]\tBarcode profiling Started at " + str(time) + "\n")
     #bam_file = args.bam_file
     #bed_file = args.bed_file
     if not os.path.exists(args.output) : os.mkdir(args.output)
     #output_file = os.path.join(args.output,"overlapping_reads.bam")
     worker.getReadsinRegion()
     worker.buildSamIndex()
     #barcode_file = os.path.join(args.output,"reads_barcode_profile.txt")
     worker.buildcov()
     #tenkbfile="/seq/schatz/sramakri/NA24385_GRCh37/NA24385_GRCh37_chr20_shit_new_barcodes.bed" 
     #overlap_file=os.path.join(args.output,"barcode_overlaps_between_regions.txt")
     worker.count_barcode_overlaps()
     time = str(datetime.now())
     sys.stderr.write("[status]\tBarcode profiling finished at " + str(time) + "\n")#, flush=True)

