#!/bin/bash
#usage: ./runall.sh


if [ ! -e runlist ]; then
  echo "Error: runlist is not exists. make: ln -s YourRunlistFile runlist"
fi

here=`pwd`
if [ ! -d "$here/runs" ]; then
  mkdir $here/runs
fi


while read line
do
    [[ $line = \#* ]] && continue
    lineArr=($line)
    echo ${lineArr[1]}
    runNr=${lineArr[0]}
    runName=${lineArr[1]}
    mkdir $here/runs/$runNr
    if [ ! -e $here/runs/$runNr/$runName.list ]; then
      ls $here/data/$runName | grep \.bin > $here/runs/$runNr/$runName.list
    fi
    time $here/read $here/runs/$runNr/$runName.list $here/data/$runName/ $here/runs/$runNr/out.root ${lineArr[0]} ${lineArr[2]} ${lineArr[3]} ${lineArr[4]} ${lineArr[5]} ${lineArr[6]} ${lineArr[7]} ${lineArr[8]} ${lineArr[9]}
done < runlist



