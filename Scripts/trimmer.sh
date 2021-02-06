#!/usr/bin/env bash

# This script runs trimmomatic in paired-end mode on all fq.gz files in the
# current directory
# The script assumes that pairs are named with the same prefix except that
# they are numbered with either a 1 or 2.
# EX: A01_1.fq.gz, A01_2.fq.gz
#
# In this example the ILLUMINACLIP line is assuming Nextera adatpers based
# on FastQC analysis
# More on parameters: http://www.usadellab.org/cms/?page=trimmomatic

inputDir=./raw_reads/
outputDir=./trimmed_reads/
threads=8
for i in ${inputDir}*.fq.gz
do
	prefix=$(basename $i _1.fq.gz)

	java -jar ./Trimmomatic-0.39/trimmomatic-0.39.jar PE \
		-threads $threads \
		-phred33 \
		${readDir}${prefix}_1.fq.gz \
		${readDir}${prefix}_2.fq.gz \
		${outputDir}${prefix}_paired1.fq.gz \
		${outputDir}${prefix}_unPaired1.fq.gz \
		${outputDir}${prefix}_paired2.fq.gz \
		${outputDir}${prefix}_unPaired2.fq.gz \
		ILLUMINACLIP:./Trimmomatic-0.39/adapters/NexteraPE-PE.fa:2:30:10:2:keepBothReads \
		LEADING:36 \
		TRAILING:41 \
		MINLEN:36
done
