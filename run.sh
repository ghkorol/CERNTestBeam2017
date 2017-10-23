#!/bin/bash

here=`pwd`

if [ ! -d "$here/runs" ]; then
  mkdir $here/runs
fi

while read line
do
    [[ $line = \#* ]] && continue
    lineArr=($line)
    if [ "${lineArr[0]}" = "$1" ]; then 
      runName=${lineArr[1]}
    fi
done < ./runlist

mkdir $here/runs/$1
if [ ! -e $here/runs/$1/$runName.list ]; then
  ls $here/data/$runName
  ls $here/data/$runName | grep \.bin > $here/runs/$1/$runName.list
fi
time $here/read $here/runs/$1/$runName.list $here/data/$runName/ $here/runs/$1/out.root


# #valgrind --trace-children=yes --tool=massif time ./../read ./$1/$1.list ./../data/$1/ ./$1/$1.root  
# #run for memory check
# #HEAPCHECK=normal LD_PRELOAD=/usr/lib/libtcmalloc.so ./../read ./$1/$1.list ./../data/$1/ ./$1/$1.root
