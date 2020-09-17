#!/bin/bash -l
list=$1
ref=$2
wc=$(wc -l ${list} | awk '{print $1}')
x=1
while [ $x -le $wc ] 
do
	string="sed -n ${x}p ${list}" 
	str=$($string)
 	var=$(echo $str | awk -F"\t" '{print $1}')   
	set -- $var
	c1=$1
 	echo "#!/bin/bash -l
 	#SBATCH -o slurm_outs/covcomp-%j.out
	#SBATCH -J covcomp
	module load OpenMPI/2.1.2
	module load GCC/6.4.0-2.28
	module load SAMtools
	module load bedtools
 	bedtools coverage -d -a NSbuffer.bed -b ${c1} > ${c1}.perbase" > ${c1}NSperbase.sh
	sbatch -t 24:00:00 -p med --mem=10G ${c1}NSperbase.sh
	sleep 2
	x=$(( $x + 1 ))
done
