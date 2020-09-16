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

### Index reference genome if needed:

    samtools index <reference.fasta>

### Align sequences to reference using BWA, sort, filter, and marking + removing duplicates using BWA + SAMtools:

This can be done for an individual sequence files ('sample.R1.fastq.gz' and 'sample.R2.fastq.gz') and a given reference ('reference.fasta') using these commands:

    bwa mem <reference.fasta> <sample>.R1.fastq.gz <sample>.R2.fastq.gz > <sample>.aln-pe.sam
    samtools view -Sb -o <sample>.aln-pe.bam <sample>.aln-pe.sam
    samtools sort -n -o <sample>.sort.bam <sample>.aln-pe.bam
    samtools view -f 0x2 -b <sample>.sort.bam > <sample>.sort.filt.bam
    samtools fixmate -m <sample>.sort.filt.bam <sample>.fixmate.bam
    samtools sort -o <sample>.positionsort.bam <sample>.fixmate.bam
    samtools markdup -r <sample>.positionsort.bam <sample>.rmdup.bam

This can also be implemented in parallel over a list of sequences ('bamlist') using the run_WGAalign.sh script on a SLURM system:
  
    sh run_WGAalign.sh <bamlist> <reference.fasta>

### Merge bam files for samples split across multiple lanes if needed:

For two individual files ('file1.rmdup.bam and file2.rmdup.bam'):

    samtools merge mergebams/${str} <folder1>/${str} <folder2>/${str} > ${str}.sh
    
If you have two folders of files with identical filenames ('mergelist') to be merged you can use the merge.sh script (alter the foldernames first) on a SLURM system to merge bam files in parallel:

    sh merge.sh mergelist

## Calculate coverage and write to a file called 'covout' for each file ('name')

    echo "<name>" >> covout
    samtools depth <name>.rmdup.bam |  awk '{sum+=$3} END { print "Average (covered sites) = ",sum/NR}' >> covout
    samtools depth -a <name>.rmdup.bam |  awk '{sum+=$3} END { print "Average (whole genome) = ",sum/NR}' >> covout

## For a list of bam files ('bamlist') calculate genotype likelihoods in ANGSD

    angsd -b bamlist -GL 2 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -minInd 12 -minMaf 0.05 -nThreads 10 -out angsdgl_WGS

## Move beagle file to PCAngsd folder, call genotypes with PCAngsd

    python pcangsd.py -beagle angsdgl_WGS.beagle.gz -geno 0.9 -o arkgenosRaptEcr -threads 10
    
## Estimate covariance matrix and individual admixture proportions

    python pcangsd.py -beagle angsdgl_WGS.beagle.gz -admix -o WGS_admix -threads 10
    
## Estimate covariance matrix and perform selection scan

    python pcangsd.py -beagle angsdgl_WGS.beagle.gz -selection 1 -sites_save -o WGS_sel -threads 10
