#!/bin/bash

rm read
g++ read.C `root-config --libs --cflags` -o read
#g++ read.C `root-config --libs --cflags` -o read  -ltcmalloc