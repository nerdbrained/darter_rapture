#!/bin/bash -l
#you will need to change folder names to reflect your file structure
list=$1
wc=$(wc -l ${list} | awk '{print $1}')
x=1
while [ $x -le $wc ] 
do
	string="sed -n ${x}p ${list}" 
	str=$($string)
 	echo "#!/bin/bash -l
 	#SBATCH -o slurm_outs/merge-%j.out
	#SBATCH -J merge
	module load OpenMPI/2.1.1
	module load GCC/6.4.0-2.28
	module load SAMtools
 	samtools merge mergebams/${str} <folder1>/${str} <folder2>/${str}" > ${str}.sh 
	sbatch -t 3:00:00 -p med --mem=4G ${str}.sh
	sleep 2
	x=$(( $x + 1 ))
done
