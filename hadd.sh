#!/bin/bash
#usage: ./hadd.sh

if [ ! -e runlist ]; then
  echo "Error: runlist is not exists. make: ln -s YourRunlistFile runlist"
fi
here=`pwd`

#rm temp.root
rm hadd.root
while read line
do
    [[ $line = \#* ]] && continue
    lineArr=($line)
    runNr=${lineArr[0]}
    haddLine=$haddLine" ./runs/$runNr/out.root "
#     if [ ! -e temp.root ]
#     then
#       cp $here/runs/$runNr/out.root temp.root
#     else
#       hadd hadd.root temp.root $here/runs/$runNr/out.root
#     fi
#     mv hadd.root temp.root
done < runlist
hadd hadd.root $haddLine
cp hadd.root CERNTestBeam2017.root
