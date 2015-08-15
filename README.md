Karect
======

KAUST Assembly Read Error Correction Tool  

Installation
============

Extract the source code from the tarball using:
tar -xzf {SOURCE}.tar.gz
  
Change your directory to the extracted code, then run the commands:  
make
make install  

Instructions
============

To get all instructions, run the program:  
karect  

Test Data and Running Example
=============================

Karect can accept as input any fasta/fastq file of assembly reads:  
Running example used in the paper of correcting Staphylococcus aureus Illumina reads:  

1) Download the files frag_1.fastq.gz and frag_2.fatstq.gz (and genome.fasta if you need to evaluate results) from:  
http://gage.cbcb.umd.edu/data/Staphylococcus_aureus/Data.original/  

2) Decompress frag_1.fastq.gz and frag_2.fatstq.gz by:  
gunzip frag_1.fastq.gz  
gunzip frag_2.fastq.gz  

3) Use Karect to correct the read sequences (modify file paths if needed):  
karect -correct -threads=12 -matchtype=hamming -celltype=haploid -inputfile=./frag_1.fastq -inputfile=./frag_2.fastq  
which produces the corrected read files: karect_frag_1.fastq and karect_frag_2.fastq  

4) If you need to evaluate correction accuracy using the reference genome (genome.fasta):
  
4a) First, align original reads to the reference genome, to produce the file ./align.txt  
karect -align -threads=12 -matchtype=hamming -inputfile=./frag_1.fastq -inputfile=./frag_2.fastq -refgenomefile=./genome.fasta -alignfile=./align.txt
  
4b) Second, evaluate the correction accuracy to produce the file ./eval.txt  
karect -eval -threads=12 -matchtype=hamming -inputfile=./frag_1.fastq -inputfile=./frag_2.fastq -resultfile=karect_frag_1.fastq -resultfile=karect_frag_2.fastq -refgenomefile=./genome.fasta -alignfile=./align.txt -evalfile=./eval.txt  

Author
======

Amin Allam  
amin.allam@kaust.edu.sa  

Reference
=========

Currently submitted to Bioinformatics, under the title:  
Karect: Accurate Correction of Substitution, Insertion and Deletion Errors for Next-generation Sequencing Data  
