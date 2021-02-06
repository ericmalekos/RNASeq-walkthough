#!/usr/bin/env bash

# This script runs ngmerge in adapter mode  on all fq.gz files in the
# current directory. Inputs must be paired-ends
# The script assumes that pairs are named with the same prefix except that
# they are numbered with either a 1 or 2.
# EX: A01_1.fq.gz, A01_2.fq.gz
#
# The option '-u' option increases maximum Qscore to 41 (from 40) which is seen
# in Illumina >=1.8 sequencing.
#
# https://github.com/jsh58/NGmerge

readDir=raw_reads/
minRead=31
threads=8
maxQ=41
outDir=adapter_trimmed_reads/
outName=trim
for i in ${readDir}*_1.fq.gz
do
        prefix=$(basename $i _1.fq.gz)
        echo $prefix

        ./NGmerge/NGmerge \
                -a \
                -z \
                -v \
                -e $minRead \
                -u $maxQ \
                -n $threads \
                -1 ${readDir}${prefix}_1.fq.gz \
                -2 ${readDir}${prefix}_2.fq.gz \
                -o ${outDir}${prefix}_${outName}
done
