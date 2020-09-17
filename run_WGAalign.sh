#!/bin/bash -l
list=$1
ref=$2
wc=$(wc -l ${list} | awk '{print $1}')
x=1
while [ $x -le $wc ] 
do
	string="sed -n ${x}p ${list}" 
	str=$($string)
 	var=$(echo $str | awk -F"\t" '{print $1, $2}')   
	set -- $var
	c1=$1
	c2=$2
	c3=$(echo $c1 | awk -F'_' '{print $1}')
 	echo "#!/bin/bash -l
 	#SBATCH -o slurm_outs/01a_align_flt-%j.out
	#SBATCH -J alignflt
	module load OpenMPI/2.1.1
	module load GCC/6.4.0-2.28
	module load SAMtools
	module load BWA/0.7.17
 	bwa mem $ref ${c1} ${c2} > ${c3}.aln-pe.sam
	samtools view -Sb -o ${c3}.aln-pe.bam ${c3}.aln-pe.sam
	samtools sort -n -o ${c3}.sort.bam ${c3}.aln-pe.bam
	samtools view -f 0x2 -b ${c3}.sort.bam > ${c3}.sort.filt.bam
	samtools fixmate -m ${c3}.sort.filt.bam ${c3}.fixmate.bam
	samtools sort -o ${c3}.positionsort.bam ${c3}.fixmate.bam
	samtools markdup -r ${c3}.positionsort.bam ${c3}.rmdup.bam" > ${c3}.sh 
	sbatch -t 10:00:00 -p med --mem=4G ${c3}.sh
	sleep 2
	x=$(( $x + 1 ))
done
