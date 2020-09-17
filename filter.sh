#!/bin/bash -l
list=$1
wc=$(wc -l ${list} | awk '{print $1}')
x=1
while [ $x -le $wc ] 
do
	string="sed -n ${x}p ${list}" 
	str=$($string)
 	var=$(echo $str | awk -F"\t" '{print $1}')   
	set -- $var
	c1=$1
	c2=$(echo $c1 | awk -F'.' '{print $1}')
 	echo "#!/bin/bash -l
 	#SBATCH -o slurm_outs/covcomp-%j.out
	#SBATCH -J covcomp
	module load OpenMPI/2.1.2
	module load GCC/6.4.0-2.28
	module load SAMtools
 	samtools view -b -h -L allbuffer.bed -o bamfilter/${c2}.filter.bam bams/${c1} " > ${c1}.sh
	sbatch -t 24:00:00 -p med --mem=10G ${c1}.NSsh
	sleep 2
	x=$(( $x + 1 ))
done
