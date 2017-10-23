#!/bin/bash

rm read
g++ read.C `root-config --libs --cflags` -o read
#compile for memory test
#g++ read.C `root-config --libs --cflags` -o read  -ltcmalloc
