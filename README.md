rescueMisreadIndex
------------------

Small C program that from gzipped / non compressed paired-end/single-end fastq files basecalled with at least the Illumina 1.8+ pipeline, look for indexed reads where bases were misread.

Compilation using kseq.h written by Heng Li and zlib.h:

	gcc Wall -g -O2 -o rescueMisReadIndex rescueMisReadIndex.c -lz

or simply type:
	
	make

Inputs:
------

	rescueMisReadIndex - 1 R1.fq.gz [-2 R2.fq.gz] -i index -m mismaches -n Ns [-t] [-p] 

	where mandatory arguments are:

	    -1 <first read input fastq filename> (gziped or not)
	    -i <index> (6 nucleotides expected)
	    -m <# mismaches> how many mismatches are allowed
	    -n <# Ns> how many Ns are allowed

	and optional arguments are:

	    -2 <second read input fastq filename> (gziped or not)
	    -t do not ouput result files, only stats files to see how much you could recover with a set of options
	    -p prefix, useful is you don't have rights to write in the folder where the inputs files are. A slash will added if not supplied

An example of effective run with an output prefix:

	rescueMisreadIndex -1 lane4_Undetermined_L004_R1_001.fastq.gz -2 lane4_Undetermined_L004_R2_001.fastq.gz -i ACAGTG -m 0 -m 1 -p mydir/

Use no arguments or -h to displays a short usage

Outputs:
-------

compressed fastq files of reads in sync. Filenames are appended with the chosen parameters:

    read1.fq_0_1_ACAGTG.gz
    read2.fq_0_1_ACAGTG.gz

a stat file contains a summary:

    read1.fq_0_1_ACAGTG.stats

Appendix:
--------

You may used the following AWK line to sum up the statistics:

	awk '{if($1~/^read/){sumRead+=$2};if($1~/^rescued/){sumGood+=$2}} END {print "tot read\t"sumRead"\ntot rescued\t"sumGood"\npcent\t"sumGood/sumRead*100}' *0_2_*stats

You can also the kind of bash script to loop on the different 4-million-read paired-end files:

	#!/bin/bash
	RUN=/basecalls/hiseq/Undetermined_indices/
	OUT=/data/
	INDEX=ATCACG
	cd $RUN
	for folder in Sample_lane? ; do
	      cd $folder 
	      for file in *_R1*gz ; do 
	              # first rescue with 2N 
	              rescueMisReadIndex1.1 -1 $file -2 ${file/_R1_/_R2_} -i $INDEX -m 0 -n 2 -p $OUT/$folder/
	              # then dry run with 1 N
	              rescueMisReadIndex1.1 -1 $file -2 ${file/_R1_/_R2_} -i $INDEX -m 0 -n 1 -t -p $OUT/$folder/
	      done
	      cd ..
	done

