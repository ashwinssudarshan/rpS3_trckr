# rpS3_trckr
for tracking of assembled rpS3 genes across metagenomes

The MIT License (MIT) Copyright (c) 2017 Alexander J Probst

Permission is hereby granted, free of charge, to any person obtaining a copy of this software and associated documentation files (the "Software"), to deal in the Software without restriction, including without limitation the rights to use, copy, modify, merge, publish, distribute, sublicense, and/or sell copies of the Software, and to permit persons to whom the Software is furnished to do so, subject to the following conditions:

The above copyright notice and this permission notice shall be included in all copies or substantial portions of the Software.

THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE SOFTWARE.

# Readme

Dependencies:
pullseq executable as "pullseq"
hmmsearch executable as "hmmsearch"
usearch64 executable as "usearch"
Ruby toolbox from https://github.com/AJProbst/tlbx
HMM databases provided in this repository
Set path for HMM databases and the toolbox in line 13 and 14 of the script before using it!

How it works:
The bash script finds rpS3 sequences of Archaea, Bacteria, and Eukaryotes in metagenome assemblies and clusters them at 99% identity, which corresponds to species level (Sharon et al., 2013 https://www.ncbi.nlm.nih.gov/pubmed/22936250).

Output files:
all_rpS3.fa contains protein sequences of all rpS3 genes detected in the samples.
OTU_table_from_assembly_coverage.csv is a comma-delimted table of all OTU clusters based on rpS3 sequences detected in the asesmblies. 0 indidcate that the rpS3 sequences was not present in the assembly.
final_rpS3_scaffolds.fna contains nucleotide sequences of the longest scaffolds per OTU cluster.
