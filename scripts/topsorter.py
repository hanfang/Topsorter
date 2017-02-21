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

# the class for the data structure
class TreeNode(object):
    def __init__(self, chr):
        self.chr = chr
        self.start = None
        self.end = None
        self.next = []


class topSorter:
    'parse & filter vcf files, construct & traverse DAGs'
    def __init__(self, inVcf, windowSize):
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

    # function for reading and filtering vcf files
    def readVcf(self):
        # import vcf
        self.vcfReader = vcf.Reader(open(self.inVcf, 'r'))
        self.prefix = os.path.splitext(os.path.basename(self.inVcf))[0]
        # remove alt contigs
        contigsToExclude = "^((?!random|Un|EBV|MT|M).)*$"
        self.chrs = [i for i in self.vcfReader.contigs.keys() if re.match(contigsToExclude, i)]
        # self.chrs = ["chr20"]
        #samples = self.vcfReader.samples
        # [temp] filter out variants without END
        hashset = set()
        for var in self.vcfReader:
            #call = variant.genotype(samples[0])
            if 'END' in var.INFO and 'PAIR_COUNT' in var.INFO and var.INFO['SVTYPE'] in ("DEL", "DUP", "INV"):
                if abs(var.INFO['END'] - var.POS) > 10000:
                    # if var.INFO['PAIR_COUNT'] >= 10: #  or var.samples['1']['SU'] >= 10:
                    id = (var.CHROM, var.POS, var.INFO['END'])
                    if id not in hashset and var.CHROM in self.chrs:
                        self.largeVariants.append(var)
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
            chrEnd = self.chrSizes[chr]
            if var.POS-self.windowSize-1 >=0 and var.INFO['END']+self.windowSize-1 <= chrEnd:
                windowSize = self.windowSize
            else:
                windowSize = min(var.POS, chrEnd-var.INFO['END'])
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
        with open(self.prefix + "/" + self.prefix + ".barcodes.bed", 'r') as fl:
            next(fl)
            for line in fl:
                row = line.split("\t")
                [chr, start, end, varId, b2l, b2r, b2a, l2r, l2a, r2a, uniqBA, uniqLR]= row
                self.barcodes[chr][varId] = [int(b2l), int(b2r), int(b2a), int(l2r), int(l2a), int(r2a), int(uniqBA), int(uniqLR)]
        fl.close()

    def buildIntervalTree(self):
        # sort the variants by start coordinates
        self.largeVariants = sorted(self.largeVariants, key=lambda x: (x.CHROM, x.POS))
        # build interval tree for each chr
        for chr in self.chrs:
            variants = [var for var in self.largeVariants if var.CHROM == chr]
            intervals = [(var.POS, var.INFO['END'], (var.ID, var.INFO['SVTYPE'])) for var in variants] # (start, end, data)
            tree = IntervalTree(Interval(*iv) for iv in intervals)
            self.treeDic[chr] = tree

    def overlapFraction(self, varA, varB):
        Asize, Bsize = varA.end - varA.begin, varB.end - varB.begin
        overlap = min(varA.end, varB.end) - max(varA.begin, varB.begin)
        fraction = float(overlap) / float(min(Asize, Bsize))
        return fraction

    def splitNodes(self, chr):
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
        # last variant
        lastVar = treeLst[-1]
        if lastVar not in existedNodes:
            DAG, leftNode, rightNode = self.nonOverlapGraph(DAG, rightNode, None, lastVar)
        #for i in list(DAG.edges()):
        #    print(i)
        # plotting
        plt.figure()# figsize=(12, 12))
        labels = {}
        #print(DAG.nodes())
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
        svId, svType = currVar.data
        [b2l, b2r, b2a, l2r, l2a, r2a, uniqBA, uniqLR] = weights[svId]
        node1 = leftNode if leftNode else Interval(0, currVar.begin - 1, (svId, "REF"))
        node2 = currVar
        node3 = rightNode if rightNode else Interval(currVar.end, currVar.end + 1, (svId, "REF"))
        DAG.add_node(node1)
        DAG.add_node(node2)
        DAG.add_node(node3)
        DAG.add_edge(node1, node2, weight=b2l)
        DAG.add_edge(node2, node3, weight=r2a)
        if svType == "DEL":
            DAG.add_edge(node1, node3, weight=b2a)
        elif svType == "DUP":
            nodeCopy = Interval(currVar.begin, currVar.end, (svId, "DupCopy"))
            DAG.add_edge(node2, nodeCopy, weight=uniqLR)
            DAG.add_edge(nodeCopy, node3, weight=uniqBA)
        elif svType == "INV":
            nodeInv = Interval(currVar.end, currVar.begin, (svId, "InvFlip"))
            DAG.add_edge(node1, nodeInv, weight=b2r)
            DAG.add_edge(nodeInv, node3, weight=l2a)
        elif svType == "BND":
            DAG.add_edge(node1, node3, weight=b2a)
        return DAG, node1, node3

    def overlapGraph(self, DAG=None, leftNode=None, rightNode=None, currVar=None, nextVar=None, weights=None):
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

    ## TODO: Change me
    def createNodes(self):
        # scan the chromosomes and split by variants
        lastPos = collections.defaultdict(int)
        for chr in self.chrs:
            i = 0
            variants = [var for var in self.largeVariants if var.CHROM == chr]
            for var in variants:
                # add ref node
                refNode = (var.CHROM, i, var.POS - 1, "REF", var.ID)
                self.graph[chr].append(refNode)
                # add variant node
                varNode = (var.CHROM, var.POS - 1, var.INFO['END'] - 1, var.INFO['SVTYPE'], var.ID)
                self.graph[chr].append(varNode)
                i = var.INFO['END'] - 1
                lastPos[chr] = i  # keep track of the last pos
        # add the last node
        for chr in self.chrs:
            lastNode = (chr, lastPos[chr], self.chrSizes[chr]-1, "REF")
            self.graph[chr].append(lastNode)

    def createWeightedDAG(self, chr):
        DAG = nx.DiGraph()
        # add nodes
        currChr = self.graph[chr]
        weights = self.barcodes[chr]
        # initialize the graph with equal weight
        for i in range(1, len(currChr)-1):
            varId = currChr[i][4]
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
            elif curr[i][3] == "BND": # currently as DELs, TBD B2A, unique
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

    def hierarchy_pos(self, G, root, width=1., vert_gap=0.2, vert_loc=0, xcenter=0.5):
        '''If there is a cycle that is reachable from root, then result will not be a hierarchy.

           G: the graph
           root: the root node of current branch
           width: horizontal space allocated for this branch - avoids overlap with other branches
           vert_gap: gap between levels of hierarchy
           vert_loc: vertical location of root
           xcenter: horizontal location of root
        '''

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
        #edge_labels = {(n1,n2): dag[n1][n2]['weight'] for (n1,n2) in dag.edges()}
        #pos_higher = {}
        #y_off = 1  # offset on the y axis
        # change label positions
        #for k, v in pos.items():
        #    pos_higher[k] = (v[0], v[1]+y_off)
        # start plotting
        plt.figure()#figsize=(12, 10))
        nx.draw_networkx(dag, pos=pos, node_size=30, labels=labels, font_size = 4)#, with_labels=False, arrows=False)
        plt.title(chr)
        plt.axis('off')
        plt.gcf()
        plt.savefig(self.prefix + "/DAGs/dag." + chr + ".pdf")
        plt.clf()
        plt.close()

    def allDAGs(self):
        # create a list of chromosome to analyze
        #chroms = list(set(self.chrs).intersection(list(self.graph.keys())))
        chroms = list(set(self.chrs)) #.intersection(list(self.graph.keys())))
        if "chrM" in chroms:
            chroms.remove("chrM")
        # create dag
        for chr in chroms:
            sys.stderr.write("[execute]\tAnalyzing "+ str(chr) + "\n")
            #dag = self.createWeightedDAG(chr)
            dag = self.splitNodes(chr)
            if dag:
                self.DAGs[chr] = dag
        # find longest path
        for chr in chroms:
            if chr in self.DAGs:
                order, longest_weighted_path = self.findLongestPath(chr)
                self.orders[chr] = order
                self.paths[chr] = longest_weighted_path
                self.passedVars[chr] = set([])
                longestVars = set([i.data[0] for i in longest_weighted_path])
                allDelBnds = set([i.data[0] for i in order if i.data[1] == "DEL" or i.data[1] == "BND"])
                falseDelBnds = set([i.data[0] for i in longest_weighted_path if i.data[1] == "DEL" or i.data[1] == "BND"])
                trueDelBnds = allDelBnds.difference(falseDelBnds)
                dups = set([i.data[0] for i in longest_weighted_path if i.data[1] == "DUP_COPY"])
                invs = set([i.data[0] for i in longest_weighted_path if i.data[1] == "INV_FLIP"])
                self.passedVars[chr] = trueDelBnds | dups | invs
                #self.passedVars[chr] = set(
                #    [i[4] for i in longest_weighted_path[:-1] if i.data[1] in ('DEL', "DUP_COPY", "INV_FLIP", "BND")])
        # draw dags
        sys.stderr.write("[execute]\tDrew graphs for each chromosome\n")  # , flush=True)
        for chr in chroms:
            if chr in self.DAGs:
                self.drawDAG(chr)

    # exporting vcf files
    def exportVcf(self):
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
                    good.append(int(var.INFO["PAIR_COUNT"]))
                else:
                    var.INFO["Topsorter"] = "F"
                    bad.append(int(var.INFO["PAIR_COUNT"]))
                vcf_writer.write_record(var)
        good, bad = np.array(good), np.array(bad)
        data_to_plot = [good, bad]
        print(np.mean(good), np.mean(bad))
        fig = plt.figure(1, figsize=(9, 6))
        ax = fig.add_subplot(111)
        bp = ax.boxplot(data_to_plot)
        fig.savefig('boxplot.png', bbox_inches='tight')

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
        sys.stderr.write("[status]\tTopsorter started at " + str(time) + "\n")#, flush=True)
        sys.stderr.write("[status]\tReading the vcf file: " + str(args.i) + "\n")#, flush=True)
        worker.readVcf()
        worker.createFolder()
        sys.stderr.write("[execute]\tExporting the large SVs to vcf and flanking regions to a bed file\n")#, flush=True)
        worker.exportVcfBed()
        sys.stderr.write("[execute]\tBuilding interval trees for the variants\n")#, flush=True)
        worker.buildIntervalTree()
        worker.parseBed()
        #worker.createNodes()
        sys.stderr.write("[execute]\tConstructing the graph and finding the longest paths\n")#, flush=True)
        worker.allDAGs()
        time = str(datetime.now())
        sys.stderr.write("[execute]\tExporting the vcf file\n")#, flush=True)
        worker.exportVcf()
        sys.stderr.write("[status]\tTopsorter finished " + str(time) + "\n")#, flush=True)
    # print usage message if any argument is missing
    else:
        sys.stderr.write("[error]\tmissing argument")
        parser.print_usage()