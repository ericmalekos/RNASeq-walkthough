#!/usr/bin/env bash

# This script runs bwa mem alignment on reads which are then piped through
# samtools to convert to sorted bam files
# The script assumes that pairs are named with the same prefix except that
# they are numbered with either a 1 or 2.
# EX: A01_1.fq.gz, A01_2.fq.gz
#
#
# https://github.com/jsh58/NGmerge

readDir=raw_reads/
threads=12
outName=trim
outDir=raw_alignments/
for i in {readDir}*_1.fq.gz
do
	prefix=$(basename $i _1.fq.gz)

	bwa mem -t $threads \
	./Reference/hg38_index \
	${readDir}${prefix}_1.fq.gz \
	${readDir}${prefix}_2.fq.gz | \
	samtools sort \
	-@ $threads \
	-O bam \
	-o ${outDir}${prefix}.bam

done
