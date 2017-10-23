#!/bin/bash
#usage: ./run.sh runNr
#       where runNr is a run number from runlist

here=`pwd`
if [ ! -d "$here/runs" ]; then
  mkdir $here/runs
fi


runNr=$1
while read line
do
    [[ $line = \#* ]] && continue
    lineArr=($line)
    if [ "${lineArr[0]}" = "$runNr" ]; then 
      runName=${lineArr[1]}
    fi
done < ./runlist

mkdir $here/runs/$runNr
if [ ! -e $here/runs/$runNr/$runName.list ]; then
  ls $here/data/$runName | grep \.bin > $here/runs/$runNr/$runName.list
fi
time $here/read $here/runs/$runNr/$runName.list $here/data/$runName/ $here/runs/$runNr/out.root


# #valgrind --trace-children=yes --tool=massif time ./../read ./$1/$1.list ./../data/$1/ ./$1/$1.root  
# #run for memory check
# #HEAPCHECK=normal LD_PRELOAD=/usr/lib/libtcmalloc.so ./../read ./$1/$1.list ./../data/$1/ ./$1/$1.root
