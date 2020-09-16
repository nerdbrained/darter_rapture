# WGS Alignment and Analysis Pipeline
  
Pipeline for aligning WGS data to a given reference genome, extracting coverage data, and performing ANGSD/PCAngsd analyses

Note that many of these scripts are taken or modified from Ryan Peek's GitHub (https://ryanpeek.github.io)

Requirements: SAMtools, BWA, ANGSD, pcangsd

This pipeline assumes that raw forward (R1) and reverse (R2) shotgun sequences from a given Illumina are already demultiplexed and compiled in a folder, along with fasta file for reference genome.

## Align sequences to a reference

### Navigate to your sequence folder: 

    cd /<sequencefolder>

### Make list of all sequence files in folder:

    ls *R1* | sed "s/\.fq//g" > bamlist1

    ls *R2* | sed "s/\.fq//g" > bamlist2

    paste bamlist? > bamlist

### Index reference genome if needed

    samtools index <reference.fasta>

### Align sequences to reference using BWA, sort, filter, and marking + removing duplicates using SAMtools.

This can be done for an individual sequence files <sample1>.R1.fq.gz <sample1>.R2.fq.gz and a given reference <reference.fasta> using these commands:

    bwa mem <reference.fasta> <sample1>.R1.fq.gz <sample1>.R2.fq.gz > <sample1>.aln-pe.sam
    samtools view -Sb -o <sample1>.aln-pe.bam <sample1>.aln-pe.sam

    samtools sort -n -o <sample1>.sort.bam <sample1>.aln-pe.bam

    samtools view -f 0x2 -b <sample1>.sort.bam > <sample1>.sort.filt.bam

    samtools fixmate -m <sample1>.sort.filt.bam <sample1>.fixmate.bam

    samtools sort -o <sample1>.positionsort.bam <sample1>.fixmate.bam

    samtools markdup -r <sample1>.positionsort.bam <sample1>.rmdup.bam

This can also be implemented in parallel over a list of sequences <bamlist> using the run_WGAalign.sh script on a SLURM system:
  
    sh run_WGAalign.sh <bamlist> <reference.fasta>

