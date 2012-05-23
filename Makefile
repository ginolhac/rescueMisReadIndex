CC=gcc
CFLAGS=-g -Wall -O2

rescueMisReadIndex:rescueMisReadIndex.c kseq.h
		$(CC) $(CFLAGS) rescueMisReadIndex.c -o $@ -lz 

clean:
		rm -fr *.o a.out rescueMisReadIndex *~ *.a 

