#!/bin/bash
mkdir $1
if [ ! -e $1.list ]
then
  #ls ./../data/$1 | grep \.bin | grep $2 > $1.list
  ls ./../data/$1 | grep \.bin > ./$1/$1.list
fi

./../read ./$1/$1.list ./../data/$1/ ./$1/$1.root
#HEAPCHECK=normal LD_PRELOAD=/usr/lib/libtcmalloc.so /home/korol/WC/read ./$1/$1.list /home/korol/WC/data/$1/ ./$1/$1.root
