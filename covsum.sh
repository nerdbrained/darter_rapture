#!/bin/bash -l
list=filterlist
echo "individual,5x,10x,20x,50x,mapped2rapt" >> indivcovsum.csv
seqno=$(wc -l ${list} | awk '{print $1}')
x=1
while [ $x -le $seqno ] 
do
	string="sed -n ${x}p ${list}" 
	seqname=$($string)
 	five=$(awk '$4>4' ${seqname}.cov | wc -l)
	ten=$(awk '$4>9' ${seqname}.cov | wc -l)
	twenty=$(awk '$4>19' ${seqname}.cov | wc -l)
	fifty=$(awk '$4>49' ${seqname}.cov | wc -l)
	total=$(awk 'NR > 3 { print $4 }' ${seqname}.cov | paste -sd+ - | bc)
	fq=$(cat ${seqname}.fq | wc -l)
	let filt=$fq/4
	echo $seqname $five $ten $twenty $fifty $total | tr ' ' , >> indivcovsum.csv
	x=$(( $x + 1 ))
done
