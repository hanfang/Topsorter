## README
####  

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

### barcode_profiling.py

#### Given the list of all barcodes present in a regions left and right of each breakpoints for each SVs that are labeled as before, left, right and after , we generated a bed file that including all the 4 regions per SVs each of which 10kb long. 
- Next we used samtools to extract the reads overlapping these regions. 
- For each of these regions we collect the information of which barcodes occur and how many reads are supporting that barcode. 
- In the end a new bed file is generated holding the so obtained informations.

   


