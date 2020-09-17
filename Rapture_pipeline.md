# Rapture Alignment and Analysis Pipeline
  
Pipeline for aligning Rapture (RADseq + sequence capture) data to a given reference genome, extracting coverage data, and performing ANGSD/PCAngsd analyses

Note that many of these scripts are taken or modified from Ryan Peek's GitHub (https://ryanpeek.github.io)

Requirements: Must have Stacks, SAMtools, bedtools, BWA, ANGSD, Python with Numpy installed at the system level, R for plotting results.

PCAngsd comes from here:

    git clone https://github.com/Rosemeis/pcangsd.git
    cd pcangsd/
    python setup.py build_ext --inplace

## Pre-process and demultiplex sequences

First step: run Flip2BeRAD python script to remove reads that do not contain a barcode on either forward or reverse read and properly orient reads that do.

Need Flip2BeRAD script (Git Repository: https://github.com/tylerhether/Flip2BeRAD) and text file of barcode sequences ('barcodelist') in directory ('flip_dir') with raw forward and reverse sequencing files ('R1.fastq' & 'R2.fastq'). With current version of BestRAD barcodes need to have exact matches (-m 0):

    cd /<flip_dir>
    python ./Flip2BeRAD.py -f ./R1.fastq -r ./R2.fastq -b ./<barcodelist> -c TGCA -m 0 -o 2
    
After flipping, I found that I needed to trim the first two sites (GG) from all forward reads. This script removes the first two characters from sequences and from quality scores.
    
    sed '2~2s/^\(.\{2\}\)//' filtered_forward.fastq > filtered_fwd_trim2.fastq
    
Move filtered reads to clonefilter folder:
    
    mkdir /<flip_dir>/clonefilter
    cd /<flip_dir>
    mv filtered* /<flip_dir>/clonefilter
    cd /<flip_dir>/clonefilter
    mkdir clonefilter_out

Use Stacks to remove clones:

    clone_filter -1 ./filtered_fwd_trim2.fastq -2 ./filtered_reverse.fastq -i fastq -o ./clonefilter_out -D

Make directory for demultiplexed reads:

    cd <flip_dir>
    mkdir processRADtags
    cd <flip_dir>/processRADtags
    mkdir processed_samples

Demultiplexing with Stacks - also need a tab-delimited file ('barcode_key') with two columns (barcode sequence and sample ID);

    cd <flip_dir>/processRADtags
    process_radtags -1 ./filtered_fwd_trim2.1.fq -2 ./filtered_reverse.2.fq -i fastq -e sbfI -b ./<barcode_key>  -o ./processed_samples --barcode_dist_2 3 -r -q

## Align sequences to a reference

Move all demultiplexed .fq files (can pool from multiple Rapture runs) and reference genome sequence ('reference.fasta') to a new folder an navigate to that folder: 

    cd /<sequencefolder>

Make list of all sequence files in folder:

    ls *R1* | sed "s/\.fq//g" > bamlist1

    ls *R2* | sed "s/\.fq//g" > bamlist2

    paste bamlist? > bamlist

Index reference genome if needed:

    samtools index <reference.fasta>

## Align sequences to reference using BWA, sort and filter using SAMtools

I recommend running using the run_align.sh script (written for a SLURM system) to do this in parallel:
  
    sh run_WGAalign.sh <bamlist> <reference.fasta>
    
## Align Rapture baits on reference genome and create .bed files

I have separated these into baits for short loci only on one side of the restriction site and long loci spanning both sides of the restriction site
Align with BWA, sort and index with samtools:

    bwa mem <reference>.fasta <shortRapturebaits>.fasta > shortalign.sam
    samtools view -Sb -o shortalign.bam shortalign.sam
    samtools sort -o shortalign.sort.bam shortalign.bam
    samtools index shortalign.sort.bam

    bwa mem <reference>.fasta <longRapturebaits>.fasta > longalign.sam
    samtools view -Sb -o longalign.bam longalign.sam
    samtools sort -o longalign.sort.bam longalign.bam
    samtools index longalign.sort.bam
    
Create .bed file and merge all overlapping/adjacent baits:

    bedtools bamtobed -i shortalign.sort.bam > shortalign.bed
    bedtools merge -i shortalign.bed > shortmerge.bed

    bedtools bamtobed -i longalign.sort.bam > longalign.bed
    bedtools merge -i longalign.bed > longmerge.bed
    
Add a +500 bp buffer for short loci:

    samtools faidx <reference>.fasta
    awk -v OFS='\t' {'print $1,$2'} <reference>.fasta.fai > genomeFile.txt
    slopBed -i shortmerge.bed -g genomeFile.txt -r 500 -l 0 > shortbuffer.bed
    
Add +/- 500 bp buffer for long loci:

    slopBed -i longmerge.bed -g Ecr_genomeFile.txt -r 500 -l 500 > longbuffer.bed

Combining .bed files is very simple:

    cat shortbuffer.bed longbuffer.bed > allbuffer.bed


## Calculate coverage and filter .bam files according to .bed files

Make list of bam files for calculating coverage and filtering

    ls *.sort.flt1.bam > filterlist

Run bedtools and generate coverage report per locus in parallel on a SLURM system - covcomp.sh script (this generates one file per individual):

    sh covcomp.sh filterlist allbuffer.bed
    
Script for summarizing coverage files - covsum.sh:

    sh covsum.sh
    
Script for calculating per-base coverage over all loci in parallel on a SLURM system - perbase.sh: 

    sh perbase.sh filterlist allbuffer.bed

Script for filtering all .bam files in parallel on a SLURM system - filter.sh:

    sh filter.sh filterlist

## Generating individual sequence files

This can be done in Stacks - we generated a reduced subset of .bam files for phylogenetic/phylogeographic analyses. All files need to have a simple naming format ('samplename.bam') and you need to provide a tab-delimited popmap file with the samplename prefix and the identifier you want it to have in the sequence file. Stacks command: 

    ref_map.pl --samples bams/ --popmap bam_indmap -o stacks_ind -X "populations: --fasta_samples --phylip --phylip_var --phylip_var_all"

## Phylogenetic analyses and plotting

Running a simple concatenated maximum likelihood analyis in iqtree:

    bin/iqtree -s phyfile -B 1000

For time-calibrating tree and plotting trees see: [darterrapture_phyloscripts.R](darterrapture_phyloscripts.R)

## For a list of bam files ('bamlist') calculate genotype likelihoods in ANGSD
We did this locally and had to increase maximum # of open files. Change the -minInd parameter to 1/2 of total number of individuals.

    ulimit -n 5000
    angsd -b bamlist -GL 2 -doGlf 2 -doMajorMinor 1 -SNP_pval 1e-6 -doMaf 1 -minInd 890 -minMaf 0.05 -nThreads 10 -out angsdgl_Rapt

## Move beagle file to PCAngsd folder, call genotypes with PCAngsd

    python pcangsd.py -beagle angsdgl_Rapt.beagle.gz -geno 0.9 -o genosRapt -threads 10
    
## Estimate covariance matrix and individual admixture proportions

    python pcangsd.py -beagle angsdgl_Rapt.beagle.gz -admix -o Rapt_admix -threads 10
    
## Estimate covariance matrix and perform selection scan

    python pcangsd.py -beagle angsdgl_Rapt.beagle.gz -selection 1 -sites_save -o Rapt_sel -threads 10
    
## PCAngsd output is in numpy format - run these in python to read and output data in .csv form

To read genotypes and output a summary of per indiviual of # SNPs homozygous for the reference and alternate alleles, heterozygous, or missing data:

    import numpy
    data=numpy.load('genosRapt.geno.npy')
    numpy.savetxt("genosRapt.csv",data,delimiter=",")
    homalt = numpy.count_nonzero(data == 2, axis=1)
    homref = numpy.count_nonzero(data == 0, axis=1)
    homhet = numpy.count_nonzero(data == 1, axis=1)
    hommis = numpy.count_nonzero(data == -9, axis=1)
    genosum=numpy.column_stack((homref,homhet,homalt,hommis))
    numpy.savetxt("Rapt_genosum.csv",genosum,delimiter=",")

To read admixture results and output Q matrix as a CSV:

    data=numpy.load('Rapt_admix.admix.Q.npy')
    numpy.savetxt("Rapt_admix.csv",data,delimiter=",")

To read selection results and output selection scan statistics for each locus along each PC axis:

    data=numpy.load('Rapt_sel.selection.npy')
    numpy.savetxt("Rapt_sel.csv",data,delimiter=",")

## Plotting results

To plot admixture results (barplot and piecharts on a map), see [plot_admix.R](plot_admix.R)

To plot selection results (Manhattan plot), see [plot_selection.R](plot_selection.R)
