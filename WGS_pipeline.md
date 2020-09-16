<h1> WGS Alignment and Analysis Pipeline
  
Pipeline for aligning WGS data to a given reference genome, extracting coverage data, and performing ANGSD/PCAngsd analyses
Note that many of these scripts are taken or modified from Ryan Peek's GitHub (https://ryanpeek.github.io)
This pipeline assumes that raw forward (R1) and reverse (R2) shotgun sequences from a given Illumina are already demultiplexed and compiled in a folder

<h2> Align sequences to a reference

<h3> Navigate to your sequence folder: 

cd /<sequencefolder>

<h3> Make list of all sequence files in folder:

ls *R1* | sed "s/\.fq//g" > bamlist1
ls *R2* | sed "s/\.fq//g" > bamlist2
paste bamlist? > bamlist
