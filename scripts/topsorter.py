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
import re
import argparse
import vcf
import collections
import networkx as nx
import numpy as np
import uuid
from networkx.drawing.nx_agraph import graphviz_layout
from datetime import datetime
from intervaltree import Interval, IntervalTree


class topSorter:
    """ add me
    """
    def __init__(self, inVcf, windowSize):
        """
        :param inVcf: input vcf file
        :param windowSize: window size for barcode counting
        """
        self.inVcf  = inVcf
        self.windowSize = windowSize
        # self.sub_vcf = []
        self.largeVariants = []
        self.chrSizes = {}
        self.chrs = []
        self.graph = collections.defaultdict(list)
        self.DAGs, self.orders, self.paths, self.passedVars = {}, {}, {}, {}
        self.barcodes = collections.defaultdict(dict)
        self.treeDic = {}

    def readVcf(self):
        """
        reading and filtering vcf files
        :return: None
        """
        # import vcf
        self.vcfReader = vcf.Reader(open(self.inVcf, 'r'))
        self.prefix = os.path.splitext(os.path.basename(self.inVcf))[0]
        # remove alt contigs
        contigsToExclude = "^((?!random|Un|EBV|MT|M).)*$"
        if self.vcfReader.contigs.keys():
            self.chrs = [i for i in self.vcfReader.contigs.keys() if re.match(contigsToExclude, i)]
        else:
            self.chrs = sorted(list(set([var.CHROM for var in self.vcfReader])))
            self.chrs = [i for i in self.chrs if re.match(contigsToExclude, i)]
        chrsHash = set(self.chrs)
        #samples = self.vcfReader.samples
        # [temp] filter out variants without END
        hashset = set()
        self.vcfReader = vcf.Reader(open(self.inVcf, 'r'))
        for var in self.vcfReader:
            #call = variant.genotype(samples[0])
            if 'END' in var.INFO and var.INFO['SVTYPE'] in {"DEL", "DUP", "INV", "BND"}: # and 'PAIR_COUNT' in var.INFO
                if abs(var.INFO['END'] - var.POS) > 10000:
                    id = (var.CHROM, var.POS, var.INFO['END'])
                    if id not in hashset and var.CHROM in chrsHash:
                        self.largeVariants.append(var)
                        hashset.add(id)
        # create a dict of chromosome: size
        if not self.vcfReader.contigs:
            thisDir, thisFn = os.path.split(__file__)
            thisDir = "/seq/schatz/hfang/Develop/Topsorter"
            fn = os.path.join(thisDir, "data", "grch38.chrom.sizes.txt")
            with open(fn, 'r') as f:
                lines = f.readlines()
                for line in lines:
                    k, v = line.split("\t")
                    self.chrSizes[k] = int(v)
        else:
            for k, v in self.vcfReader.contigs.items():
                self.chrSizes[k] = v[1]

    def createFolder(self):
        """
        create folders
        :return: on disk
        """
        cmd = 'mkdir -p ' + self.prefix + ';' + 'mkdir -p ' + self.prefix + '/DAGs'
        os.system(cmd)

    def exportVcfBed(self):
        """
        export bed file to disk
        :return: on disk
        """
        # export the subset of large SVs to a vcf
        if int(sys.version_info.major) >= 3:
            vcf_writer = vcf.Writer(open(self.prefix + "/" + self.prefix + ".large_sv.vcf", 'w'), self.vcfReader)
        for record in self.largeVariants:
            vcf_writer.write_record(record)
        # export the flanking regions as a bed file
        bedStr = ""
        for var in self.largeVariants:
            chr = var.CHROM
            chrEnd = self.chrSizes[chr]
            if var.POS - self.windowSize - 1 >= 0 and var.INFO['END'] + self.windowSize - 1 <= chrEnd:
                windowSize = self.windowSize
            else:
                windowSize = min(var.POS, chrEnd - var.INFO['END'])
                print("[warning]\t" + str(var.ID) + " - distance to chromosome end smaller than specified: " + str(windowSize) +"\n")
            if windowSize <= 500:
                print("[warning]\t" + str(var.ID) + " - distance to chromosome end is less than 500bp\n")
            leftMost, rightMost = max(0, var.POS-windowSize-1), min(chrEnd, var.INFO['END']+windowSize-1)
            coordinate = str(var.POS) + "," + str(var.INFO['END'])
            before = chr +"\t"+ str(leftMost) +"\t"+ str(var.POS-1) +"\t"+ var.ID +","+ coordinate + ",before" + "\n"
            left = chr +"\t"+ str(var.POS-1) +"\t"+ str(var.POS+windowSize-1) +"\t"+ var.ID +","+ coordinate + ",left" + "\n"
            right = chr +"\t"+ str(var.INFO['END']-windowSize-1) +"\t"+ str(var.INFO['END']-1) +"\t"+ var.ID +","+ coordinate + ",right" + "\n"
            after = chr +"\t"+ str(var.INFO['END']-1) +"\t"+ str(rightMost) +"\t"+ var.ID +","+ coordinate + ",after" +"\n"
            bedStr += before + left + right + after
        # write bed file to local
        bedOut = open(self.prefix + "/" + self.prefix + "." + str(self.windowSize) +".bed", "w")
        bedOut.write(bedStr)

    def parseBed(self):
        """
        parse the input barcode bed file
        :return: None
        """
        with open(self.prefix + "/" + self.prefix + "." + str(self.windowSize) + ".barcodes.bed", 'r') as fl:
            next(fl)
            for line in fl:
                row = line.split("\t")
                [chr, start, end, varId, b2l, b2r, b2a, l2r, l2a, r2a, uniqBA, uniqLR]= row
                self.barcodes[chr][varId] = [int(b2l), int(b2r), int(b2a), int(l2r), int(l2a), int(r2a), int(uniqBA), int(uniqLR)]
        fl.close()

    def buildIntervalTree(self):
        """
        build interval tree for variants
        :return: None
        """
        # sort the variants by start coordinates
        self.largeVariants = sorted(self.largeVariants, key=lambda x: (x.CHROM, x.POS))
        # build interval tree for each chr
        for chr in self.chrs:
            variants = [var for var in self.largeVariants if var.CHROM == chr]
            intervals = [(var.POS, var.INFO['END'], (var.ID, var.INFO['SVTYPE'])) for var in variants] # (start, end, data)
            tree = IntervalTree(Interval(*iv) for iv in intervals)
            self.treeDic[chr] = tree

    def overlapFraction(self, varA, varB):
        """
        calculate reciprocal overlapping fraction of two variants
        :param varA: variant A as Interval object
        :param varB: variant B as Interval object
        :return: fraction (float)
        """
        Asize, Bsize = varA.end - varA.begin, varB.end - varB.begin
        overlap = min(varA.end, varB.end) - max(varA.begin, varB.begin)
        fraction = float(overlap) / float(min(Asize, Bsize))
        return fraction

    def buildGraph(self, chr):
        """
        build a graph for a chromosome
        :param chr: chromosome name (str)
        :return: dag (networkx graph object)
        """
        """
        :param chr:
        :return:
        """
        tree = self.treeDic[chr]
        if not tree:
            return
        treeLst = list(sorted(tree)) # sort by start/begin
        weights = self.barcodes[chr]
        existedNodes = set([])
        DAG = nx.DiGraph()
        leftNode, rightNode = None, None
        for i in range(len(treeLst)):
            currVar = treeLst[i]
            overlaps, envelops = set(tree.search(currVar)), set(tree.search(currVar, strict=True))
            if len(overlaps) == 1:
                # no overlap
                if currVar not in existedNodes:
                    DAG, leftNode, rightNode = self.nonOverlapGraph(DAG, rightNode, None, currVar, weights)
                    existedNodes.add(currVar)
            elif len(overlaps) > 1: #  and len(envelops) == 1:
                # partial overlap
                if currVar not in existedNodes:
                    overlaps.remove(currVar)
                    childVars = []
                    # define left/right ref nodes
                    leftSuperNode = Interval(currVar.begin - 2, currVar.begin - 1, (currVar.data[0], "REF"))
                    rightSuperNode = Interval(currVar.end + 1, currVar.end + 2, (currVar.data[0], "REF"))
                    if not rightNode:
                        rightNode = Interval(0, currVar.begin - 2, (currVar.data[0], "REF"))
                    DAG.add_edge(rightNode, leftSuperNode, weight=1)
                    # overlap variants
                    for childVar in overlaps:
                        faction = self.overlapFraction(currVar, childVar)
                        if faction > 0.1 and childVar not in existedNodes:
                            childVars.append(childVar)
                    if not childVars:
                        DAG, leftNode, rightNode = self.nonOverlapGraph(DAG, rightNode, None, currVar, weights)
                    else:
                        for childVar in sorted(childVars):
                            DAG, leftNode, rightNode = self.overlapGraph(DAG, leftSuperNode, rightSuperNode, currVar, childVar, weights)
                            existedNodes.add(childVar)
                    existedNodes.add(currVar)
        # plotting
        plt.figure()
        labels = {}
        for node in DAG.nodes():
            labels[node] = (node.data[0].split(":")[0], node.data[1])
        nx.draw_networkx(DAG, labels=labels)
        plt.title("Example")
        plt.axis('off')
        plt.gcf()
        plt.savefig("example.pdf")
        plt.clf()
        plt.close()
        sys.stderr.write("[status]\tFinished creating the graph for " + str(chr) + "\n")#, flush=True)
        return DAG
        #lst = [i for i in list(DAG.nodes()) if i]
        #print(sorted(lst))

    def nonOverlapGraph(self, DAG=None, leftNode=None, rightNode=None, currVar=None, weights=None):
        """
        build non overlapping graph locally
        :param DAG: directed acyclic graph (networkx graph object)
        :param leftNode: node on the left, Interval object
        :param rightNode: node on the right, Interval object
        :param currVar: current variant, Interval object
        :param weights: dictionary of weights for each variant, dict
        :return: DAG, node1, node3
        """
        svId, svType = currVar.data
        [b2l, b2r, b2a, l2r, l2a, r2a, uniqBA, uniqLR] = weights[svId]
        node1 = leftNode if leftNode else Interval(0, currVar.begin - 1, (svId, "REF"))
        node2 = currVar
        node3 = rightNode if rightNode else Interval(currVar.end, currVar.end + 1, (svId, "REF"))
        DAG.add_node(node1)
        DAG.add_node(node2)
        DAG.add_node(node3)
        DAG.add_edge(node1, node2, weight=b2l)
        if svType == "DEL":
            DAG.add_edge(node1, node3, weight=b2a)
            DAG.add_edge(node2, node3, weight=r2a)
        elif svType == "DUP":
            nodeCopy = Interval(currVar.begin, currVar.end, (svId, "DupCopy"))
            DAG.add_edge(node2, node3, weight=uniqLR)
            DAG.add_edge(node2, nodeCopy, weight=l2r)
            DAG.add_edge(nodeCopy, node3, weight=r2a) # undecided
        elif svType == "INV":
            nodeInv = Interval(currVar.end, currVar.begin, (svId, "InvFlip"))
            DAG.add_edge(node2, node3, weight=r2a)
            DAG.add_edge(node1, nodeInv, weight=b2r)
            DAG.add_edge(nodeInv, node3, weight=l2a)
        elif svType == "BND":
            DAG.add_edge(node2, node3, weight=r2a)
            DAG.add_edge(node1, node3, weight=uniqBA)
        return DAG, node1, node3

    def overlapGraph(self, DAG=None, leftNode=None, rightNode=None, currVar=None, nextVar=None, weights=None):
        """
        build overlapping graph locally when two or more variants are overlapping
        :param DAG: directed acyclic graph (networkx graph object)
        :param leftNode: node on the left (Interval object)
        :param rightNode: node on the right (Interval object)
        :param currVar: current variant (Interval object)
        :param nextVar: next variant (Interval object)
        :param weights: dictionary of weights for each variant (dict)
        :return: DAG, leftNode, rightNode
        """
        # curr var
        nodeUL = Interval(currVar.begin - 1, currVar.begin, (currVar.data[0], "REF"))
        nodeUR = Interval(currVar.end, currVar.end + 1, (currVar.data[0], "REF"))
        DAG.add_node(nodeUL)
        DAG.add_node(nodeUR)
        DAG.add_edge(leftNode, nodeUL, weight=1)
        DAG.add_edge(nodeUR, rightNode, weight=1)
        DAG, tmp1, tmp2 = self.nonOverlapGraph(DAG, nodeUL, nodeUR, currVar, weights)
        # next var
        nodeLL = Interval(nextVar.begin - 1, nextVar.begin, (nextVar.data[0], "REF"))
        nodeLR = Interval(nextVar.end, nextVar.end + 1, (nextVar.data[0], "REF"))
        DAG.add_node(nodeLL)
        DAG.add_node(nodeLR)
        DAG.add_edge(leftNode, nodeLL, weight=1)
        DAG.add_edge(nodeLR, rightNode, weight=1)
        DAG, tmp1, tmp2 = self.nonOverlapGraph(DAG, nodeLL, nodeLR, nextVar, weights)
        return DAG, leftNode, rightNode

    def findLongestPath(self, chr):
        """
        find the longest path for a chr
        :param chr: chromosome name (str)
        :return: order (list), longest_weighted_path (list)
        """
        # topological sorting and find the longest path
        dag = self.DAGs[chr]
        order = nx.topological_sort(dag)
        longest_weighted_path = self.longestWeightedPath(dag)
        # sys.stderr.write("[execute]\tPerforming topological sorting for " + str(chr) + "\n")
        # sys.stderr.write("[execute]\tFinding longest paths in " + str(chr) + "\n")
        print("[results]\tLongest weighted path in", chr, "\n", longest_weighted_path)
        return order, longest_weighted_path

    def longestWeightedPath(self, G):
        """
        find longest path in a weighted DAG
        :param G: dag (networkx graph object)
        :return: longest_weighted_path (list)
        """
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

    def hierarchy_pos(self, G, root, width=1., vert_gap=0.2, vert_loc=0, xcenter=0.5):
        """
        calculate hierarchy positions for nodes in a dag
        :param G: the graph
        :param root: the root node of current branch
        :param width: horizontal space allocated for this branch - avoids overlap with other branches
        :param vert_gap: gap between levels of hierarchy
        :param vert_loc: vertical location of root
        :param xcenter: horizontal location of root
        """
        def h_recur(G, root, width=4., vert_gap=8, vert_loc=2, xcenter=0.5, pos=None, parent=None, parsed=[]):
            if (root not in parsed):
                parsed.append(root)
                if pos == None:
                    pos = {root: (xcenter, vert_loc)}
                else:
                    pos[root] = (xcenter, vert_loc)
                neighbors = G.neighbors(root)
                if len(neighbors) != 0:
                    dx = width
                    nextx = xcenter - dx / 2 - width / 2
                    for neighbor in neighbors:
                        nextx += dx
                        pos = h_recur(G, neighbor, width=dx, vert_gap=vert_gap, vert_loc=vert_loc - vert_gap,
                                      xcenter=nextx, pos=pos, parent=root, parsed=parsed)
            return pos
        return h_recur(G, root, width=4., vert_gap=8, vert_loc=2, xcenter=0.5)

    def drawDAG(self, chr):
        """
        draw a dag
        :param chr: chromosome name (str)
        :return: on disk
        """
        labels = {}
        for node in self.orders[chr]:
            #labels[node] = node[3]
            if node.data[1] != "REF":
                if node.data[1] == "DEL":
                    labels[node] = "-"
                elif node.data[1] == "DUP" or node.data[1] == "DUP_COPY":
                    labels[node] = "+"
                elif node.data[1] == "INV":
                    labels[node] = "~"
                else:
                    labels[node] = ">"
            else:
                labels[node] = " "
        sourceNode = self.orders[chr][0]
        sinkNode = self.orders[chr][-1]
        labels[sourceNode] = "Source"
        labels[sinkNode] = "Sink"
        # extract dag for each chr
        dag = self.DAGs[chr]
        pos = self.hierarchy_pos(dag, sourceNode)
        plt.figure()
        nx.draw_networkx(dag, pos=pos, node_size=30, labels=labels, font_size = 4)#, with_labels=False, arrows=False)
        plt.title(chr)
        plt.axis('off')
        plt.gcf()
        plt.savefig(self.prefix + "/DAGs/dag." + chr + ".pdf")
        plt.clf()
        plt.close()

    def allDAGs(self):
        """
        execute the graphical analysis for all the chromosomes
        :return: None
        """
        # create a list of chromosome to analyze
        #chroms = list(set(self.chrs).intersection(list(self.graph.keys())))
        chroms = sorted(list(set(self.chrs))) #.intersection(list(self.graph.keys())))
        if "chrM" in chroms:
            chroms.remove("chrM")
        # create dag
        for chr in chroms:
            # sys.stderr.write("[execute]\tAnalyzing "+ str(chr) + "\n")
            dag = self.buildGraph(chr)
            if dag:
                self.DAGs[chr] = dag
        # find longest path
        for chr in chroms:
            if chr in self.DAGs:
                order, longest_weighted_path = self.findLongestPath(chr)
                self.orders[chr] = order
                self.paths[chr] = longest_weighted_path
                self.passedVars[chr] = set([])
                allDelBnds = set([i.data[0] for i in order if i.data[1] == "DEL" or i.data[1] == "BND"])
                falseDelBnds = set([i.data[0] for i in longest_weighted_path if i.data[1] == "DEL" or i.data[1] == "BND"])
                trueDelBnds = allDelBnds.difference(falseDelBnds)
                dups = set([i.data[0] for i in longest_weighted_path if i.data[1] == "DupCopy"])
                invs = set([i.data[0] for i in longest_weighted_path if i.data[1] == "InvFlip"])
                self.passedVars[chr] = trueDelBnds | dups | invs
        # draw dags
        sys.stderr.write("[execute]\tDrew graphs for each chromosome\n")
        for chr in chroms:
            if chr in self.DAGs:
                self.drawDAG(chr)

    def exportVcf(self):
        """
        export vcf file with topsorter tags
        :return: on disk
        """
        ## TODO: output in sorted order
        if int(sys.version_info.major) >= 3:
            vcf_writer = vcf.Writer(open(self.prefix + "/" + self.prefix + ".topsorter.vcf", 'w'), self.vcfReader)
        elif int(sys.version_info.major) == 2:
            vcf_writer = vcf.Writer(open(self.prefix + "/" + self.prefix + ".topsorter.vcf", 'wb'), self.vcfReader)
        # write to file
        good, bad = [] ,[]
        for var in self.largeVariants:  # sub_vcf
            if var.CHROM in self.chrs:
                if self.passedVars[var.CHROM] and var.ID in self.passedVars[var.CHROM]:
                    var.INFO["Topsorter"] = "T"  # placeholder
                    # good.append(int(var.INFO["SU"][0]))
                else:
                    var.INFO["Topsorter"] = "F"
                    # bad.append(int(var.INFO["SU"][0]))
                vcf_writer.write_record(var)
        #print(np.mean(np.array(good)), np.mean(np.array(bad)))
        #good = [i for i in good if i < 200]
        #bad = [i for i in bad if i < 200]
        #good, bad = np.array(good), np.array(bad)
        #data_to_plot = [good, bad]
        #print(np.mean(good), np.mean(bad))
        #fig = plt.figure(1, figsize=(9, 6))
        #ax = fig.add_subplot(111)
        #bp = ax.boxplot(data_to_plot)
        #fig.savefig('boxplot.png', bbox_inches='tight')

# the main process
if __name__ == '__main__':
    parser = argparse.ArgumentParser(description='createGraph.py - a script for constructing a graph from a vcf file')
    parser.add_argument("-i", help="input vcf file [REQUIRED]",required=True)
    parser.add_argument("-w", help="window size of counting barcodes, default: 1000 [OPTIONAL]", default=1000, required=False, type=int)
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
        windowSize = args.w
        # main
        worker = topSorter(vcfFn, windowSize)
        time = str(datetime.now())
        sys.stderr.write("[status]\tTopsorter started at " + str(time) + "\n")
        sys.stderr.write("[status]\tReading the vcf file: " + str(args.i) + "\n")
        worker.readVcf()
        worker.createFolder()
        sys.stderr.write("[execute]\tExporting the large SVs to vcf and flanking regions to a bed file\n")
        worker.exportVcfBed()
        sys.stderr.write("[execute]\tBuilding interval trees for the variants\n")
        worker.buildIntervalTree()
        sys.stderr.write("[execute]\tParsing the barcode bed file\n")
        worker.parseBed()
        sys.stderr.write("[execute]\tConstructing the graph and finding the longest paths\n")
        worker.allDAGs()
        time = str(datetime.now())
        sys.stderr.write("[execute]\tExporting the vcf file\n")
        worker.exportVcf()
        sys.stderr.write("[status]\tTopsorter finished " + str(time) + "\n")
    # print usage message if any argument is missing
    else:
        sys.stderr.write("[error]\tmissing argument")
        parser.print_usage()