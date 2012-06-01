#include <zlib.h>  
#include <stdio.h>  
#include<string.h>
#include <stdbool.h>
#include "kseq.h"  

/*
Aurelien Ginolhac 120125
From gzipped/non gzipped paired fastq files from Illumina 1.8+
look for indexed reads were the index was misread

compile with  gcc -o rescueMisReadIndex rescueMisReadIndex.c -lz
Heng Li's kseq.h must be present in the same folder as c code.
doc here http://lh3lh3.users.sourceforge.net/parsefastq.shtml */


#define VERSION ("1.1")

// STEP 1: declare the type of file handler and the read() function  
KSEQ_INIT(gzFile, gzread)  


void distance (char *id1, char *id2, int *score, int *n) {
	
	int i;
	
	for (i=0;i<strlen(id1);i++){

		if(id1[i] == id2[i]) {
			(*score)++; // match +1
		}
		else  if ((id1[i] != 'N' ) && (id2[i] != 'N')) {
			//(*score)--; // mismatch, do nothing
		}
		else {
			(*n)++; // mismatch with one N, +1
			(*score)++;
		}
	}
	//printf("score %d #Ns %d\n", *score, *n);
	//return (score);
}

void help ( char *prog_name ) {
	fprintf(stderr, "Usage: %s - 1 R1.fq.gz [-2 R2.fq.gz] -i index -m mismaches -n Ns [-t] [-p]\n", prog_name);
	fprintf(stderr, "NOTE: The output is always gziped compressed.\n");
	fprintf(stderr, "Required Arguments:\n" );
	fprintf(stderr, "\t-1 <first read input fastq filename> (gziped or not)\n" );
	fprintf(stderr, "\t-i <index> (6 nucleotides expected)\n" );
	fprintf(stderr, "\t-m <# mismaches> how many mismatches are allowed\n" );
	fprintf(stderr, "\t-n <# Ns> how many Ns are allowed\n" );
	fprintf(stderr, "Optional:\n" );
	fprintf(stderr, "\t-2 <second read input fastq filename> (gziped or not)\n" );
	fprintf(stderr, "\t-t do not ouput result files, only stats files\n" );
	fprintf(stderr, "\t-q quiet, do not print progression info\n" );
	fprintf(stderr, "\t-p prefix to be added to output files, MUST end by a '/'\n\n" );
	exit (2);  
}

  int main(int argc, char *argv[])  {  
	gzFile fp, fp2, w1, w2, stat, m;
	kseq_t *seq, *seq2;  
	w1 = w2 = m = seq2 = fp2 = NULL;  
	int l, opt; 
	int cpt = 0;
	int good = 0;
	int bad = 0;
	int lg = 0;
	int req_args = 0;
	char *fin1 = (char *) calloc(201, sizeof(char)); // for output names
	char *fin2 = (char *) calloc(201, sizeof(char));
	char *index = (char *) calloc(10, sizeof(char));
	char *prefix = (char *) calloc(201, sizeof(char));
	char *mism_allowed = (char *) calloc(10, sizeof(char));
	char *n_allowed = (char *) calloc(10, sizeof(char));
	int mism_threshold = 0;
	int nbn_threshold = 0;
	bool paired = false;
	bool dry = false;
	bool append = false;
	bool quiet = false;
	  
	char delims[] = ":"; // for splitting the index field such as 1:N:0:ACAGTG
	if ((argc == 1)) {  
		help(argv[0]);
	}  

	// example taken from https://github.com/jstjohn/SeqPrep/blob/master/SeqPrep.c
	while( (opt=getopt( argc, argv, "1:2:i:m:n:p:tqh" )) != -1 ) {
		switch( opt ) {
	
		//REQUIRED ARGUMENTS
		case '1' :
		req_args ++;
		strcpy( fin1, optarg );
		break;
		
		case 'i' :
		req_args ++;
		strcpy( index, optarg );
		lg = strlen(index) ;
		// check length of the index provided by the user
		if( lg != 6) {
			fprintf(stderr, "index %s is not 6 characters long (%d) \n", index, lg);  
			exit (1);
		}
		break;
	
		case 'm' :
		req_args ++;
		mism_allowed = optarg;
		mism_threshold = lg - atoi(optarg);
		break;
		
		case 'n' :
		req_args ++;
		n_allowed = optarg;
		nbn_threshold = atoi(optarg);
		break;
		
		// OPTIONAL
		case '2' :
		paired=true;
		strcpy( fin2, optarg );
		break;
		
		case 'p' :
		append = true;
		strcpy( prefix, optarg );
		break;
		
		case 't' :
		dry=true;
		break;
		
		case 'q' :
		quiet=true;
		break;
		
		case 'h' :
		help(argv[0]);
		break;
		
		default :
		help(argv[0]);
		}
	}
	if(req_args < 4){
		fprintf(stderr, "Missing a required argument! (%d < 4)\n", req_args);
		help(argv[0]);
	}
	fp = gzopen(fin1, "r"); // open the file handler  
	seq = kseq_init(fp); //  initialize seq  
 	if( seq == NULL) {
		fprintf(stderr, "Cannot open file %s \n", fin1);
		exit(1);
	}	
	char *fout1 = (char *) calloc(201, sizeof(char)); // for output names
	char *fstat  = (char *) calloc(201, sizeof(char));
	char *fout2 = (char *) calloc(201, sizeof(char));
	
	// single-end mode
	if( paired) {
		fp2 = gzopen(fin2, "r"); 
		seq2 = kseq_init(fp2); 
 		if( seq2 == NULL) {
			fprintf(stderr, "Cannot open file %s \n", fin2);
			exit(1);
		}
	}
	
	if( append ) {
		strcpy(fout1, prefix);
		strcat(fout1 , fin1);
		strcpy(fstat, prefix);
		strcat(fstat , fin1);
		// single-end mode
		if( paired) {
			strcpy(fout2, prefix);
			strcat(fout2 , fin2);
		}
		fprintf(stderr, "Add prefix %s to outfiles\n", prefix);
	}
	else {
		strcpy(fout1, fin1);
		// single-end mode
		if( paired) {
			strcpy(fout2, fin2);
		}
		strcpy(fstat, fin1);
	}
	strcat(fout1 , "_");	
	strcat(fout1 , mism_allowed);	
	strcat(fout1 , "_");	
	strcat(fout1 , n_allowed);
	strcat(fout1 , "_");
	strcat(fout1 , index);
	strcat(fout1 , ".out.gz");
	strcat(fstat , "_");
	strcat(fstat , mism_allowed);
	strcat(fstat , "_");
	strcat(fstat , n_allowed);
	strcat(fstat , "_");
	strcat(fstat , index);
	strcat(fstat , ".stats");
	
	// single-end mode
	if( paired) {
		strcat(fout2 , "_");
		strcat(fout2 , mism_allowed);	
		strcat(fout2 , "_");		
		strcat(fout2 , n_allowed);
		strcat(fout2 , "_");
		strcat(fout2 , index);
		strcat(fout2 , ".out.gz");
	}
	
	if(dry) {
		printf("dry run only, writing only stats file %s\n", fstat);
	}
	else{
		w1 = gzopen(fout1, "w"); //  open the file handlers to write
 		if( w1 == NULL) {
			fprintf(stderr, "Cannot write file %s \n", fout1);
			exit(1);
		}
		if( paired ) {
			w2 = gzopen(fout2, "w");
		 	if( w2 == NULL) {
				fprintf(stderr, "Cannot write file %s \n", fout2);
				exit(1);
			}	
		}
	} 
	stat = fopen(fstat, "wb");  
	if( stat == NULL) {
		fprintf(stderr, "Cannot write file %s \n", fstat);
                return 1;
	}	
	// single-end mode
	if( paired ) {
		printf("Paired-end mode\n");
	}
	else {
		printf("Single-end mode\n");
	}
	printf("Looking for index %s in %s, #mismatch allowed %s, #Ns allowed %s\n", index, fin1, mism_allowed, n_allowed );
	
	while ((l = kseq_read(seq)) >= 0 ) { // STEP 4: read sequence  from pairs in parallel
		
		if( paired ) {
			kseq_read(seq2) ;
		}
		
		//printf("name: %s\n", seq->name.s);  
		cpt++;
		int flag = 0;
		int i=0;
		if (!seq->comment.l) fprintf(stderr, "Comment with index not defined! Check if produced by illumina 1.8+\n");  
		//printf("seq: %s\n", seq->seq.s);  
		//if (seq->qual.l) printf("qual: %s\n", seq->qual.s);  
		//printf("c: %s\n", seq->comment.s);  
		
		// start for split comment, ie index field	
		char *result = NULL;
		char *idx1 = NULL;
		// strtok modifies the string, we must copy it before it mess up
		char com1[20] ;
		strcpy(com1 , seq->comment.s);  
		result = strtok( seq->comment.s, delims );
		while( result != NULL ) {
			i++;
			//printf( "result is \"%s\" and i %d from %s\n", result, i, com1 );
			// 5th field is the index, record it
			if( i == 4) { idx1 = result; }
			result = strtok( NULL, delims );
		}
		char *idx2 = NULL;
		char com2[20] ;
		if( paired) {
			// idem R2
			char *result2 = NULL;
			strcpy(com2 , seq2->comment.s);  
			result2 = strtok( seq2->comment.s, delims );
			i=0;
			while( result2 != NULL ) {
				i++;
				if( i == 4) { idx2 = result2; }
				result2 = strtok( NULL, delims );
			}
			if(idx2 == NULL) {
			fprintf(stderr, "indexes were not read correctly for %s (%s)\n", seq2->name.s, com2);
			exit(1);
			}
		}
		if(idx1 == NULL) {
			fprintf(stderr, "indexes were not read correctly for %s\n", seq->name.s);
			exit(1);
		}
		// compare index 1 with the one expected
		int nbn=0;
		int score = 0;
		// compute distance and nb of Ns, use variable addresses to get both variables
		distance(idx1, index, &score, &nbn);
		
		//fprintf(stderr, "mm %d nbn allowed %d ", mism_threshold , nbn_threshold);
		
		if((score >= mism_threshold) && (nbn <= nbn_threshold)){
			flag++;
		}
		
		if( paired ) {
			// check if pair is ok ie, same index read in both paired files
			if( strcmp ( idx1, idx2) == 0 ){ flag++; }
			else { bad++;}
		}
		else { flag++;}
		
		//fprintf(stderr, "idx1 is %s, idx2 is %s expected %s flag %d score %d #Ns %d\n", idx1, idx2, index, flag, score, nbn);
		
		// pair can be rescued, print them
		if( flag == 2 ) {
			good++;
			if(!dry) {
				gzprintf(w1,"@%s %s\n%s\n+\n%s\n", seq->name.s, com1, seq->seq.s, seq->qual.s);
				// single-end mode
				if( paired) {
					gzprintf(w2,"@%s %s\n%s\n+\n%s\n", seq2->name.s, com2, seq2->seq.s, seq2->qual.s);
				}
			}
		}
		if((cpt % 2000000 == 0 ) && (!quiet)) {fprintf(stderr, "%8d reads processed, %8d rescued so far\n", cpt, good); }
	}  
	printf("%d seq read, %d rescued\n", cpt, good);  
	if(paired) {
		fprintf(stat,"Paired-end\nVersion %s\nfiles\t%s\t%s\nindex\t%s\n#mismatches\t%s\n#Ns\t%s\nread\t%d\nrescued\t%d\nunpaired\t%d\n", VERSION, fin1, fin2,index, mism_allowed, n_allowed, cpt, good, bad);
	}
	else {
		fprintf(stat,"Single-end\nVersion %s\nfile\t%s\nindex\t%s\n#mismatches\t%s\n#Ns\t%s\nread\t%d\nrescued\t%d\n", VERSION, fin1, index, mism_allowed, n_allowed, cpt, good);
	}
	kseq_destroy(seq); // STEP 5: destroy seq  
	kseq_destroy(seq2);
	gzclose(fp); // STEP 6: close the file handlers
	if(!dry) {	gzclose(w1); }
	if(paired) {
		gzclose(fp2);
		if(!dry) {	gzclose(w2); }   
	}
	fclose(stat);
	return 0;  
} 
