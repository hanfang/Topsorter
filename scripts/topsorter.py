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
        self.vcfReader = vcf.Reader(open(self.inVcf, 'r'))
        self.prefix = os.path.splitext(os.path.basename(self.inVcf))[0]
        # self.sub_vcf = []
        self.largeVariants = []
        self.chrSizes = {}
        self.chrs = []
        self.graph = collections.defaultdict(list)
        self.DAGs, self.orders, self.paths, self.passedVars = {}, {}, {}, {}
        self.barcodes = collections.defaultdict(dict)

    # function for reading and filtering vcf files
    def readVcf(self):
        # remove alt contigs
        #self.chrs = [i for i in self.vcfReader.contigs.keys() if "random" not in i and "Un" not in i and "EBV" not in i]
        self.chrs = ["chr20"]
        #samples = self.vcfReader.samples
        # [temp] filter out variants without END
        hashset = set()
        for variant in self.vcfReader:
            #call = variant.genotype(samples[0])
            if 'END' in variant.INFO and 'PAIR_COUNT' in variant.INFO and variant.INFO['SVTYPE'] in ("DEL", "DUP", "INV"):
                if abs(variant.INFO['END'] - variant.POS) > 10000:
                    # if variant.INFO['PAIR_COUNT'] >= 10: #  or variant.samples['1']['SU'] >= 10:
                    id = (variant.CHROM, variant.POS, variant.INFO['END'])
                    if id not in hashset:
                        self.largeVariants.append(variant)
                        hashset.add(id)
        # create a dict of chromosome: size
        for k, v in self.vcfReader.contigs.items():
            self.chrSizes[k] = v[1]

    def createFolder(self):
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
            coordinate = str(var.POS) + "," + str(var.INFO['END'])
            before = chr +"\t"+ str(var.POS-1001) +"\t"+ str(var.POS-1) +"\t"+ var.ID +","+ coordinate + ",before" + "\n"
            left   = chr +"\t"+ str(var.POS-1) +"\t"+ str(var.POS+999) +"\t"+ var.ID +","+ coordinate + ",left" + "\n"
            right = chr +"\t"+ str(var.INFO['END']-1001) +"\t"+ str(var.INFO['END']-1) +"\t"+ var.ID +","+ coordinate + ",right" + "\n"
            after = chr +"\t"+ str(var.INFO['END']-1) +"\t"+ str(var.INFO['END']+999) +"\t"+ var.ID +","+ coordinate + ",after" +"\n"
            bedStr += before + left + right + after
        # write bed file to local
        bedOut = open(self.prefix + "/" + self.prefix + ".bed", "w")
        bedOut.write(bedStr)

    def parseBed(self):
        with open(self.prefix + "/" + self.prefix + ".barcodes.bed", 'r') as fl:
            next(fl)
            for line in fl:
                row = line.split("\t")
                [chr, start, end, varId, b2l, b2r, b2a, l2r, l2a, r2a]= row
                self.barcodes[chr][varId] = [int(b2l), int(b2r), int(b2a), int(l2r), int(l2a), int(r2a)]
        fl.close()

    def createNodes(self):
        # sort the variants by start coordinates
        self.largeVariants = sorted(self.largeVariants, key=lambda x: x.POS)
        # scan the chromosomes and split by variants
        lastPos = collections.defaultdict(int)
        for chr in self.chrs:
            i = 0
            for variant in self.largeVariants:
                if variant.CHROM == chr:
                    # add ref node
                    refNode = (variant.CHROM, i, variant.POS-1, "REF", variant.ID)
                    self.graph[chr].append(refNode)
                    # add variant node
                    varNode = (variant.CHROM, variant.POS-1, variant.INFO['END']-1, variant.INFO['SVTYPE'], variant.ID)
                    self.graph[chr].append(varNode)
                    i = variant.INFO['END'] - 1
                    lastPos[chr] = i # keep track of the last pos
        # add the last node
        for chr in self.chrs:
            lastNode = (chr, lastPos[chr], self.chrSizes[chr]-1, "REF")
            self.graph[chr].append(lastNode)

    def createWeightedDAG(self, chr):
        DAG = nx.DiGraph()
        # add nodes
        curr = self.graph[chr]
        weights = self.barcodes[chr]
        # initialize the graph with equal weight
        for i in range(1, len(curr)-1):
            varId = curr[i][4]
            [b2l, b2r, b2a, l2r, l2a, r2a] = weights[varId]
            DAG.add_node(curr[i-1])
            DAG.add_node(curr[i])
            DAG.add_node(curr[i+1])
            if curr[i][3] == "DEL":
                DAG.add_edge(curr[i-1], curr[i], weight=b2l)
                DAG.add_edge(curr[i-1], curr[i+1], weight=b2a)
                DAG.add_edge(curr[i], curr[i+1], weight=r2a)
            elif curr[i][3] == "DUP":
                nodeCopy = (curr[i][0], curr[i][1], curr[i][2], "DUP_COPY", varId)
                DAG.add_node(nodeCopy)
                DAG.add_edge(curr[i-1], curr[i], weight=b2l)
                DAG.add_edge(curr[i], nodeCopy, weight=l2r)
                DAG.add_edge(curr[i], curr[i+1], weight=l2a)
                DAG.add_edge(nodeCopy, curr[i+1], weight=r2a)
            elif curr[i][3] == "INV":
                nodeInv = (curr[i][0], curr[i][2], curr[i][1], "INV_FLIP", varId)
                DAG.add_edge(curr[i-1], curr[i], weight=b2l)
                DAG.add_edge(curr[i], curr[i+1], weight=r2a)
                DAG.add_edge(curr[i-1], nodeInv, weight=b2r)
                DAG.add_edge(nodeInv, curr[i+1], weight=l2a)
            elif curr[i][3] == "BND": # currently as DELs, TBD
                DAG.add_edge(curr[i-1], curr[i], weight=b2l)
                DAG.add_edge(curr[i], curr[i+1], weight=r2a)
                DAG.add_edge(curr[i-1], curr[i+1], weight=b2a) #
            else:
                pass
        return DAG

    def findLongestPath(self, chr):
        # topological sorting
        dag = self.DAGs[chr]
        # sys.stderr.write("[execute]\tPerforming topological sorting for " + str(chr) + "\n")
        order = nx.topological_sort(dag)
        # print("[results]\t Topologically sorted order of nodes in ", chr, "\n", order, flush=True)
        # find the longest path
        sys.stderr.write("[execute]\tFinding longest paths in " + str(chr) + "\n")
        longest_weighted_path = self.longestWeightedPath(dag)
        print("[results]\tLongest weighted path in", chr, "\n", longest_weighted_path)
        return order, longest_weighted_path

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
        node, (length,_)  = max(dist.items(), key=lambda x:x[1])
        path = []
        while length > 0:
            path.append(node)
            length, node = dist[node]
        return list(reversed(path))

    def drawDAG(self, chr):
        labels = {}
        i = 0
        for node in self.orders[chr]:
            labels[node] = node[3]
            i += 1
        # pos = graphviz_layout(self.DAG)
        dag = self.DAGs[chr]
        edge_labels = {(n1,n2): dag[n1][n2]['weight'] for (n1,n2) in dag.edges()}
        pos = nx.spring_layout(dag, k=0.5, scale=0.5)
        pos_higher = {}
        y_off = 1  # offset on the y axis
        # change label positions
        #for k, v in pos.items():
        #    pos_higher[k] = (v[0], v[1]+y_off)
        # start plotting
        plt.figure(figsize=(12, 12))
        #nx.draw(dag, pos_higher, labels=labels)
        # nx.draw_networkx(dag, pos, labels=labels)
        #nx.draw_networkx(dag, pos)
        #nx.draw_networkx_edge_labels(dag, pos, edge_labels=edge_labels)
        nx.draw_networkx_nodes(dag, pos)
        nx.draw_networkx_edges(dag, pos)
        #nx.draw_networkx_labels(dag, pos, labels=labels)
        nx.draw_networkx_edge_labels(dag, pos, edge_labels=edge_labels)
        plt.title(chr)
        plt.axis('off')
        plt.gcf()
        plt.savefig(self.prefix + "/DAGs/dag." + chr + ".pdf")
        plt.clf()
        plt.close()

    def allDAGs(self):
        # create a list of chromosome to analyze
        chroms = list(set(self.chrs).intersection(list(self.graph.keys())))
        if "chrM" in chroms:
            chroms.remove("chrM")
        # create dag
        for chr in chroms:
            dag = self.createWeightedDAG(chr)
            self.DAGs[chr] = dag
        # find longest path
        for chr in chroms:
            order, longest_weighted_path = self.findLongestPath(chr)
            self.orders[chr] = order
            self.paths[chr] = longest_weighted_path
            self.passedVars[chr] = set(
                [i[4] for i in longest_weighted_path[:-1] if i[3] in ('DEL', "DUP_COPY", "INV_FLIP", "BND")])
        # draw dags
        for chr in chroms:
            self.drawDAG(chr)

    # exporting vcf files
    def exportVcf(self):
        ## TODO: output in sorted order
        if int(sys.version_info.major) >= 3:
            vcf_writer = vcf.Writer(open(self.prefix + "/" + self.prefix + ".topsorter.vcf", 'w'), self.vcfReader)
        elif int(sys.version_info.major) == 2:
            vcf_writer = vcf.Writer(open(self.prefix + "/" + self.prefix + ".topsorter.vcf", 'wb'), self.vcfReader)
        # write to file
        for var in self.largeVariants:  # sub_vcf
            if var.CHROM in self.chrs:
                if self.passedVars[var.CHROM] and var.ID in self.passedVars[var.CHROM]:
                    var.INFO["Topsorter"] = "T"  # placeholder
                else:
                    var.INFO["Topsorter"] = "F"
                vcf_writer.write_record(var)

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
        worker.parseBed()
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