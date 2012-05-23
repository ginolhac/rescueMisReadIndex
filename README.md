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
	    -p prefix, useful is you don't have rights to write in the folder where the inputs files are. MUST finish by a slash

An example of effective run with an output prefix:

	rescueMisreadIndex -1 lane4_Undetermined_L004_R1_001.fastq.gz -2 lane4_Undetermined_L004_R2_001.fastq.gz -i ACAGTG -m 0 -m 1 -p mydir/

Use no arguments or -h to displays a short usage

Outputs:
-------

compressed fastq files of reads in sync. Filenames are appended with the chosen parameters:

    read1.fq.gz_0_1_ACAGTG.out.gz
    read2.fq.gz_0_1_ACAGTG.out.gz

stats file which sum up:

    read1.fq.gz_0_1_ACAGTG.stats

