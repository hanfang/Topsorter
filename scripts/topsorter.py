#!/usr/bin/env python

## ----------------------
# Topsorter
## ----------------------
# Han Fang (hanfang.cshl@gmail.com)
# Cold Spring Harbor Laboratory
## ----------------------

from __future__ import print_function
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from io import open
import os
import sys
import argparse
import vcf
import collections
import networkx as nx
import uuid
from networkx.drawing.nx_agraph import graphviz_layout
from datetime import datetime


# the class for the data structure
class TreeNode(object):
    def __init__(self, chr):
        self.chr = chr
        self.start = None
        self.end = None
        self.next = []


class topSorter:
    'parse & filter vcf files, construct & traverse DAGs'
    def __init__(self, inVcf):
        self.inVcf  = inVcf
        # self.sub_vcf = []
        self.largeVariants = []
        self.chrSizes = {}
        self.chrs = []

    # function for reading and filtering vcf files
    def readVcf(self):
        # read a vcf file
        self.vcfReader = vcf.Reader(open(self.inVcf, 'r'))
        self.chrs = [i for i in self.vcfReader.contigs.keys() if "random" not in i and "Un" not in i and "EBV" not in i]
        #samples = self.vcfReader.samples
        # [temp] filter out variants without END
        hashset = set()
        for variant in self.vcfReader:
            #call = variant.genotype(samples[0])
            if 'END' in variant.INFO and 'PAIR_COUNT' in variant.INFO and variant.INFO['SVTYPE'] in ("DEL", "DUP", "INV"):
                if abs(variant.INFO['END'] - variant.POS) > 10000:
                    if variant.INFO['PAIR_COUNT'] >= 10: #  or variant.samples['1']['SU'] >= 10:
                        id = (variant.CHROM, variant.POS, variant.INFO['END'])
                        if id not in hashset:
                            self.largeVariants.append(variant)
                            hashset.add(id)
        # create a dict of chromosome: size
        for k, v in self.vcfReader.contigs.items():
            self.chrSizes[k] = v[1]

    def createFolder(self):
        self.prefix = os.path.splitext(os.path.basename(self.inVcf))[0]
        cmd = 'mkdir -p ' + self.prefix + ';' + 'mkdir -p ' + self.prefix + '/DAGs'
        os.system(cmd)

    def exportVcfBed(self):
        # export the subset of large SVs to a vcf
        if int(sys.version_info.major) >= 3:
            vcf_writer = vcf.Writer(open(self.prefix + "/" + self.prefix + ".large_sv.vcf", 'w'), self.vcfReader)
        for record in self.largeVariants:
            vcf_writer.write_record(record)
        # export the flanking regions as a bed file
        bedStr = ""
        for var in self.largeVariants:
            chr = var.CHROM
            before = chr +"\t"+ str(var.POS-10001) +"\t"+ str(var.POS-1) +"\t"+ var.ID + ",before\n"
            left   = chr +"\t"+ str(var.POS-1) +"\t"+ str(var.POS+9999) +"\t"+ var.ID + ",left\n"
            right = chr +"\t"+ str(var.INFO['END']-10001) +"\t"+ str(var.INFO['END']-1) +"\t"+ var.ID + ",right\n"
            after = chr +"\t"+ str(var.INFO['END']-1) +"\t"+ str(var.INFO['END']+9999) +"\t"+ var.ID + ",after\n"
            bedStr += before + left + right + after
        # write bed file to local
        bedOut = open(self.prefix + "/" + self.prefix + ".bed", "w")
        bedOut.write(bedStr)

    def parseBed(self):
        pass

    def createNodes(self):
        # sort the variants by start coordinates
        self.largeVariants = sorted(self.largeVariants, key=lambda x: x.POS)

        # scan the chromosomes and split by variants
        self.graph = collections.defaultdict(list)
        lastPos = collections.defaultdict(int)
        for chr in self.chrSizes.keys():
            i = 0
            for variant in self.largeVariants:
                if variant.CHROM == chr:
                    # add ref node
                    refNode = (variant.CHROM, i, variant.POS-1, "REF")
                    self.graph[chr].append(refNode)
                    # add variant node
                    varNode = (variant.CHROM, variant.POS-1, variant.INFO['END']-1, variant.INFO['SVTYPE'])
                    self.graph[chr].append(varNode)
                    i = variant.INFO['END'] - 1
                    lastPos[chr] = i # keep track of the last pos
        # add the last node
        for chr in self.chrSizes.keys():
            last_node = (chr, lastPos[chr], self.chrSizes[chr]-1, "REF")
            self.graph[chr].append(last_node)

    def createWeightedDAG(self, chr):
        DAG = nx.DiGraph()
        # add nodes
        curr = self.graph[chr]
        ## weights = # placeholder
        # initialize the graph with equal weight
        for i in range(1, len(curr)-1):
            DAG.add_node(curr[i-1])
            DAG.add_node(curr[i])
            DAG.add_node(curr[i+1])
            DAG.add_edge(curr[i-1], curr[i], weight=1)
            DAG.add_edge(curr[i], curr[i+1], weight=1)
            if curr[i][3] == "DEL":
                DAG.add_edge(curr[i-1], curr[i+1], weight=2)
            elif curr[i][3] == "DUP":
                nodeCopy = (curr[i][0], curr[i][1], curr[i][2], "DUP_COPY")
                DAG.add_node(nodeCopy)
                DAG.add_edge(curr[i], nodeCopy, weight=0.5)
                DAG.add_edge(nodeCopy, curr[i+1], weight=0.5)
            elif curr[i][3] == "INV":
                nodeInv = (curr[i][0], curr[i][2], curr[i][1], "INV_FLIP")
                DAG.add_edge(curr[i-1], nodeInv, weight=1)
                DAG.add_edge(nodeInv, curr[i+1], weight=1)
            elif curr[i][3] == "BND":
                DAG.add_edge(curr[i-1], curr[i+1], weight=2) # currently as DELs
            else:
                pass
        return DAG

    def findLongestPath(self, chr):
        # topological sorting
        dag = self.DAGs[chr]
        sys.stderr.write("[execute]\tPerforming topological sorting for " + str(chr) + "\n")#, flush=True)
        order = nx.topological_sort(dag)
        print("[results]\t Topologically sorted order of nodes in ", chr, "\n", order, flush=True)
        # find the longest path
        longest_weighted_path = self.longestWeightedPath(dag)
        print("[results]\t Longest weighted path in ", chr, "\n", longest_weighted_path)
        return order, longest_weighted_path
        # temp:
        # longest path without weights:
        # longest = nx.dag_longest_path(self.DAG)
        # print ("[results]\t Longest\n", longest, flush=True)
        # longestDis = nx.dag_longest_path_length(self.DAG)
        # print ("[results]\t Longest length: ", longestDis, flush=True)

    def longestWeightedPath(self, G):
        dist = {} # stores [node, distance] pair
        # start with topological sorted order
        for node in nx.topological_sort(G):
            # pairs of dist,node for all incoming edges
            pairs = [(dist[v][0] + G[v][node]['weight'], v) for v in G.pred[node]] # incoming pairs
            if pairs:
                dist[node] = max(pairs)
            else:
                dist[node] = (0, node)
        node,(length,_)  = max(dist.items(), key=lambda x:x[1])
        path = []
        while length > 0:
            path.append(node)
            length,node = dist[node]
        return list(reversed(path))

    def allDAGs(self):
        # create a list of chromosome to analyze
        chroms = list(set(self.chrs).intersection(list(self.graph.keys()))) #chroms.remove("chr21") #chroms.remove("chrY")
        # self.graph.keys()
        self.DAGs, self.orders, self.paths = {}, {}, {}
        # create dag
        for chr in chroms:
            dag = self.createWeightedDAG(chr)
            self.DAGs[chr] = dag
        # find longest path
        for chr in chroms:
            order, longest_weighted_path = self.findLongestPath(chr)
            self.orders[chr] = order
            self.paths[chr] = longest_weighted_path
        # draw dags
        for chr in chroms:
            self.drawDAG(chr)

    def drawDAG(self, chr):
        labels = {}
        i = 0
        for node in self.orders[chr]:
            labels[node] = (i, node[3])
            i += 1
        # edge_labels = {(n1,n2): self.DAG[n1][n2]['weight'] for (n1,n2) in self.DAG.edges()}
        # pos = graphviz_layout(self.DAG)
        dag = self.DAGs[chr]
        pos = nx.spring_layout(dag)
        pos_higher = {}
        y_off = 1  # offset on the y axis
        # change label positions
        for k, v in pos.items():
            pos_higher[k] = (v[0], v[1]+y_off)
        # start plotting
        plt.figure()
        nx.draw(dag, pos_higher, labels=labels)
        # nx.draw_networkx(self.DAG, pos, labels=labels)
        # nx.draw_networkx_edge_labels(self.DAG, pos, edge_labels=edge_labels)
        plt.gcf()
        plt.savefig(self.prefix + "/DAGs/dag." + chr + ".pdf")
        plt.clf()
        plt.close()

    # exporting vcf files
    def exportVcf(self):
        ## TODO: output in sorted order
        if int(sys.version_info.major) >= 3:
            vcf_writer = vcf.Writer(open(self.prefix + "/" + self.prefix + ".topsorter.vcf", 'w'), self.vcfReader)
        elif int(sys.version_info.major) == 2 :
            vcf_writer = vcf.Writer(open(self.prefix + "/" + self.prefix + ".topsorter.vcf", 'wb'), self.vcfReader)
        # write to file
        for record in self.largeVariants: # sub_vcf
            record.INFO["Topsorter"] = "T" # placeholder
            vcf_writer.write_record(record)

# the main process
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='createGraph.py - a script for constructing a graph from a vcf file')
    parser.add_argument("-i", help="input vcf file [REQUIRED]",required=True)

    # check if there is any argument
    if len(sys.argv) <= 1:
        parser.print_usage()
        sys.stderr.write("[help] please use [-h] for details")
        sys.exit(1)
    else:
        args = parser.parse_args()

    # process the file if the input files exist
    if (args.i!=None):
        vcfFn = args.i
        # main
        worker = topSorter(vcfFn)
        time = str(datetime.now())
        sys.stderr.write("[status]\tTopsorter started at " + str(time) + "\n")#, flush=True)
        sys.stderr.write("[status]\tReading the vcf file: " + str(args.i) + "\n")#, flush=True)
        worker.readVcf()
        worker.createFolder()
        sys.stderr.write("[execute]\tExporting the large SVs to vcf and flanking regions to a bed file\n")#, flush=True)
        worker.exportVcfBed()
        sys.stderr.write("[execute]\tCreating the nodes\n")#, flush=True)
        worker.createNodes()
        sys.stderr.write("[execute]\tConstructing the graph and finding the longest paths\n")#, flush=True)
        worker.allDAGs()
        sys.stderr.write("[execute]\tDrew graphs for each chromosome\n")#, flush=True)
        time = str(datetime.now())
        sys.stderr.write("[execute]\tExporting the vcf file\n")#, flush=True)
        worker.exportVcf()
        sys.stderr.write("[status]\tTopsorter finished " + str(time) + "\n")#, flush=True)
    # print usage message if any argument is missing
    else:
        sys.stderr.write("[error]\tmissing argument")
        parser.print_usage()