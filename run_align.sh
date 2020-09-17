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
 	echo "#!/bin/bash -l
 	#SBATCH -o slurm_outs/01a_align_flt-%j.out
	#SBATCH -J alignflt
	module load OpenMPI/2.1.1
	module load GCC/6.4.0-2.28
	module load SAMtools
	module load BWA/0.7.17
 	bwa mem $ref ${c1}.fq ${c2}.fq > ${c1}.aln-pe.sam
	samtools view -Sb -o ${c1}.aln-pe.bam ${c1}.aln-pe.sam
	samtools sort -o ${c1}.sort.bam ${c1}.aln-pe.bam
	samtools view -f 0x2 -b ${c1}.sort.bam > ${c1}.sort.flt1.bam" > ${c1}.sh
	sbatch -t 24:00:00 -p med --mem=4G ${c1}.sh
	sleep 2
	x=$(( $x + 1 ))
done
