#!/bin/sh

##NOTE: This script works only on linux.
# To make function on MacOS, for the sed command, add an additional "" after -i

#list of scaffold names, or however naming scheme works
scaflist=("1085" "1086" "1087" "2151" "2446" "2531" "2532" "2533" "2684" "2685" "2686" "2687")

#bash for loop to convert fasta consensus files into proper format for gene tree estimation
for i in consensus*.fa;
do
	species=${i/.fa/}; species=${species/consensus_fullgenome_bws_/}
	
	#change from 50bp/break to no break
	awk '/^>/ { print (NR==1 ? "" : RS) $0; next } { printf "%s", $0 } END { printf RS }' $i > nobreaks.$i
	
	#split the nobreak fastas up by scaffold
	awk -v species=$species '/^>/ {out = substr($1, 2) "_" species ".fa"; print > out} !/^>/ {print >> out}' nobreaks.$i
	
	#for each of the scaffolds, replace the name to match the species name
	for j in scaffold*$species.fa;
	do
		sed -i "1s/^.*$/>$species/g" $j;
	done
	
	#clean up
	rm nobreaks.$i
done


#concatenate scaffold files
for i in "${scaflist[@]}"
do
	cat scaffold_$i* > allseqs_consensus_scaffold_$i.fa
	rm scaffold_$i*
done


