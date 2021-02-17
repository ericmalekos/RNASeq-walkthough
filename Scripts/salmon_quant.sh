#!/usr/bin/env bash

# This script runs salmon quant in mapping based mode.
# The script assumes that pairs are named with the same prefix except that
# they are numbered with either a 1 or 2.
# EX: A01_1.fastq.gz, A01_2.fastq.gz
#
# Arguments:
#  libType: A autodetects library strandedness
#  index: points to the output of "salmon index ..."
#  threads: number of threads
#  softclip: allow soft clipping.
#
# There are MANY additional parameters that can be considered.
# Use "salmon quant --help-reads"
#
# Salmon docs: https://salmon.readthedocs.io/en/latest/salmon.html

threads=12
readDir=../adapter_trimmed_reads/
index=./salmon_index/
outDir=./

for read in ${readDir}*_1.fastq.gz
do
        base="$(basename "$read" _1.fastq.gz)"
        name1="$base"_1.fastq.gz
        name2="$base"_2.fastq.gz
        echo
        echo First Read:  $name1
        echo Second Read: $name2

        outDir="$base"/
        echo Output Directory: $outDir

        salmon/bin/salmon quant \
        --libType A \
        --index $index \
        --mates1 ${readDir}/${name1} \
        --mates2 ${readDir}/${name2} \
        --threads $threads \
        --softclip \
        --output ${outDir}/$"base"

done
