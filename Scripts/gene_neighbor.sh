#/usr/bin/bash

bed=$1
name=$2
number=$3

mkdir backup_dir

names=$(cat $name)
echo "===================="
echo "marker gene number < 2x+1"
echo "===================="

arr0=()
for line in $names
do
    grep -A $number -B $number  -w $line $bed > $line.txt # 原始的数据是txt
    awk "{print \$1}" $line.txt | sort | uniq > tmp1 # 统计染色体个数
    x=$(cat tmp1)
    declare -i SUM=0
    for n in $x # 计算染色体个数
    do
        let SUM+=1
    done
    rm tmp1
    if [ $SUM -eq 1 ] # 判断是否跨越了染色体或者scaffold
    then
        cp $line.txt $line.bed
        rm $line.txt
    else # 如果跨越了scaffold，提取基因所在的scaffold信息
        grep "$line" $line.txt | awk '{print $1}' > tmp.chr_name # 提取基因所在行的染色体信息
        cat tmp.chr_name | while read l
        do
            grep -w $l $line.txt > $line.bed
            mv $line.bed backup_dir/$line.bed # 移动
            rm $line.txt
            echo $line # 标记基因数目少于11的neighborhood信息
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
    arr1=("${arr1[@]}" $i) # 统计有overlap的列表
    f_name=( "${f_name[*]/$i}" ) # 向后比较，所以删除$i变量
    for j in $f_name
    do
        if [ $i!=$j ] # 这句起始没必要
        then 
            cat $i $j | sort | uniq -d > $i-$j.tmp # 判断是否重叠
            if [ -s $i-$j.tmp ] # $i-$j.tmp文件存在且不为空是，为True
            then
                arr1=(${arr1[@]} $j) # $i-$j.tmp文件不为空说明 有重叠，这是将$j添加到arr1列表中
                rm $i-$j.tmp # 删除无用数据
                f_name=( "${f_name[*]/$j}" ) #从f_name循环中删除该数据，已经进入了重叠的范围，后续会进行进一步整理
            else
                rm $i-$j.tmp # 为空的文件更需要删除
            fi
        fi
    done
    if [ ${#arr1[@]} != 1 ] # 判断重叠列表里是否为1个元素，不为1 则继续合并
    then
        a="C"
        echo ${arr1[@]} # 打印合并bed文件的数据
        for m in ${arr1[@]} # 逐个合并重叠的bed文件
        do
            cat $m >> $i.cross
            arr2=("${arr2[@]}" $m) # 对合并过的元素进行统计，后续移动数据到别的文件夹
            n2=${m%.*} # 获取文件的基因名
            a=$a-$n2
        done
        sort -n -k3 $i.cross | uniq > $a.bed
        rm $i.cross # 删除中间文件
    fi
    unset arr1 # 清空合并数据的列表
done

echo 'removing some file to backup_dir'

for p in ${arr2[@]}
do
    mv $p backup_dir/$p
done
unset arr2 # 清空arr2

# 根据文件名对上面的数据进行进一步合并
f2_name=$(ls C-*) # 统计合并过的文件
for b in $f2_name
do
    d="S"
    n4=${b%.*}
    n5=${n4#*-}
    d=$d-$n5
    f2_name=( "${f2_name[*]/$b}" ) # 向后比较，删除
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


if [ ${#arr3[@]} != 1 ] # 判断重叠列表里是否为1个元素，不为1 则继续合并
then
    for m in ${arr3[@]} # 逐个合并重叠的bed文件
    do
        mv $m backup_dir/$m
    done
fi
unset arr3 # 清空合并数据的列表
