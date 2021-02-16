#!/usr/bin/env bash

# This script runs STAR alignment on paired end reads and counts genes
# The script assumes that pairs are named with the same prefix except that
# they are numbered with either a 1 or 2.
# EX: A01_1.fastq.gz, A01_2.fastq.gz
#
# Arguments:
#  readFilesCommand: set to "zcat" if files are in ".gz" format, otherwise remove this
#  outSAMtype: sorted BAM file
#  outSAMattributes: attributes to include in alignment file. NH HI AS nM are standard.
#               I include XS because I believe it is reequired for juncBASE
#  twopassMode: default is single pass, "Basic" indicates to 2-pass mapping. I understand
#               this to be preferred for considering splice junctions
#  outFilterMultimapNmax: maximum number of loci the read is allowed to map to. 10 is default.
#               reads that map to >10 loci are not included in BAM Alignment file.
#  quantMode: GeneCounts results in output file containing gene counts the same as htseq count with
#               default settings. All library combinations are accounted for (unstranded, first/second stranded)
#               each as a separate column, so you need to pick the correct column from this output.
#               See STAR manual/htseq manual for more info.
#  runThreadN: number of threads
#  alignEndsType: alignmnent. Local is default and performs soft clipping.
#  outFileNamePrefix: there will be many output files. It's probably best to send each to its own
#               directory. This assigns the prefix to each.
#
# There are MANY additional parameters that can and should be considered especially with respect
# to dowstream analysis. Given the current parameters we will have:
#  Sorted BAM files for each pair of fastq.gz
#  ReadsPerGene.out.tab which contain gene count files for each pair. This can be used as input
#               to DESeq2 directly after selecting the appropriate colum
#
# STAR manual: https://github.com/alexdobin/STAR/blob/master/doc/STARmanual.pdf

threads=12
readDir=../adapter_trimmed_reads/
index=./M25_index/
outDir=./
for read in ${readDir}*1.fastq.gz
threads=12
readDir=../adapter_trimmed_reads/
index=./M25_index/
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

        STAR-2.7.7a/source/STAR \
        --genomeDir $index \
        --readFilesIn \
        ${readDir}/${name1} \
        ${readDir}/${name2} \
        --readFilesCommand zcat \
        --outSAMtype BAM SortedByCoordinate \
        --outSAMattributes NH HI AS nM XS \
        --twopassMode Basic \
        --outFilterMultimapNmax 10 \
        --quantMode GeneCounts \
        --runThreadN $threads \
        --alignEndsType Local \
        --outFileNamePrefix ${outDir}/$"base"

done
