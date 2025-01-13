#/usr/bin/bash
# Author Xuan Zhou

bed=$1
name=$2
number=$3 # the parameter is X

mkdir backup_dir

names=$(cat $name)
echo "===================="
echo "mark GN size < 2x+1"
echo "===================="

arr0=()
for line in $names
do
    grep -A $number -B $number  -w $line $bed > $line.txt # origin result is a txt file
    awk "{print \$1}" $line.txt | sort | uniq > tmp1 # the name of chromosome
    x=$(cat tmp1)
    declare -i SUM=0
    for n in $x # the number of chromosome
    do
        let SUM+=1
    done
    rm tmp1
    if [ $SUM -eq 1 ] # GN size < 2X+1 ?
    then
        cp $line.txt $line.bed
        rm $line.txt
    else # get the chromosome name of the scaffolding gene
        grep "$line" $line.txt | awk '{print $1}' > tmp.chr_name 
        cat tmp.chr_name | while read l
        do
            grep -w $l $line.txt > $line.bed
            mv $line.bed backup_dir/$line.bed # mv the GN to backup_dir
            rm $line.txt
            echo $line # mark the gene name in the log file
        done
        rm tmp.chr_name
    fi
done

echo "===================="
echo "neighborhood-overlap"
echo "===================="
f_name=$(ls *.bed)
for i in $f_name
do
    arr1=("${arr1[@]}" $i) # list of neighborhood-overlap
    f_name=( "${f_name[*]/$i}" )
    for j in $f_name
    do
        if [ $i!=$j ]
        then 
            cat $i $j | sort | uniq -d > $i-$j.tmp # i and j have overlap or not
            if [ -s $i-$j.tmp ] # 
            then
                arr1=(${arr1[@]} $j) # 
                rm $i-$j.tmp # 
                f_name=( "${f_name[*]/$j}" ) #
            else
                rm $i-$j.tmp #
            fi
        fi
    done
    if [ ${#arr1[@]} != 1 ] # combine the GN have overlap
    then
        a="C"
        echo ${arr1[@]} #
        for m in ${arr1[@]} #
        do
            cat $m >> $i.cross
            arr2=("${arr2[@]}" $m) # the file have been combined 
            n2=${m%.*}
            a=$a-$n2
        done
        sort -n -k2 $i.cross | uniq > $a.bed
        rm $i.cross # delete the temporary file 
    fi
    unset arr1
done

echo 'removing some file to backup_dir'

for p in ${arr2[@]}
do
    mv $p backup_dir/$p
done
unset arr2


# according to the file name to combine the GNs that have overlap
f2_name=$(ls C-*)
for b in $f2_name
do
    d="S"
    n4=${b%.*}
    n5=${n4#*-}
    d=$d-$n5
    f2_name=( "${f2_name[*]/$b}" )
    for r in $f2_name
    do
        cat $b $r | sort | uniq -d > tmp
        if [ -s tmp ]
        then
            arr3=("${arr3[@]}" $b $r)
            n3=${r%.*}
            n6=${n3#*-}
            cat $b $r >> S-$b
            d=$d-$n6
            f2_name=( "${f2_name[*]/$r}" )
        fi
        rm tmp
    done
    if [ -s S-$b ]
    then
        sort -n -k3 S-$b | uniq > $d.bed
        rm S-$b
    fi
done


if [ ${#arr3[@]} != 1 ] 
then
    for m in ${arr3[@]}
    do
        mv $m backup_dir/$m
    done
fi
unset arr3
