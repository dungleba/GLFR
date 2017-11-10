CC=g++
LOP=-o
LOPT=-O3 -funroll-loops

MAIN=./benchm
TAG=benchmark


$(MAIN).o :
	$(CC) $(LOPT) $(LOP) $(TAG) $(MAIN).cpp


