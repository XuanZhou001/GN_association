#!/usr/bin/bash

file_name=$(ls *.bed)
declare -i num=0
for i in $file_name
do
    let num+=1
    sed "s/^/$num\t&/g" $i >> all.bed
done
