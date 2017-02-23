## Topsorter
### - Graphical assessment of structrial variants using 10x genomics data Edit
####  

![alt text](https://github.com/hanfang/Topsorter/blob/master/image/topsorter.png | width=100)

### Contact
- Han Fang (hanfang.cshl@gmail.com)
- Srividya Ramakrishnan (srividya.ramki@gmail.com)
- Fritz Sedlazeck (fritz.sedlazeck@gmail.com)

--------

### Topsorter

#### Traversing & finding the longest path (most confident haplotype) in a weighted directed acyclic graph (DAG).
- At the initial construction step, it splits regions of a chromosome by structural variants and creates a weighted DAG.
- Then it updates the weights of the edges according to barcodes information (and/or other quality metrics) from the script barcode_profiles.py.
- Finally Topsorter performs topological sorting of the graph and finds the longest path (the most confident haplotype).
- Input: vcf file
- Output: PDF files of graphs for each chromosme, longest paths
- Command: `python topsorter.py $vcf ` 

-------- 

### BarcodeProfiler

#### Building barcode profile for the alignments and count the overlapping barcodes for every split region of a chromosome 
-  Extract alignments from every split region in to a bam file and index them
-  Identify the barcodes in each of these split regions and count number of reads per barcode
-  Count the barcode overlaps between the  split regions of interest
-  Input: Phased bam file from 10x data, Constructed split regions in bed format using Topsorter class (func exportVCFBed )
-  Output: directory containing files
         overlapping_reads.bam
         overlapping_reads.bam.bai
         reads_barcode_profile.txt
         barcode_overlaps_between_regions.txt
   
-  Command: ./barcode_profiles.py -bam \<input phased bam\> -bed \<input constructed bed file\> -o \<output dir name\>


