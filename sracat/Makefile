OBJS = 

CC = c++
cc = cc

FLAGS = -O3 -Wall

PROFILE = #-pg
OPENMP = #-fopenmp

# Edit the SRA_INCLUDE_PATH to point to the SRA toolkit on your system.
# This makefile assumes that all of the SRA libraries for both "ngs" and "ncbi-vdb" have
# been installed in $SRA_PATH/lib64
SRA_PATH = ../SRA

INC = -I. -I$(SRA_PATH)/include

# 1) The order of the libraries in LIBS is important for g++, do not change!
# 2) g++ requires -lpthread, while clang does not
# 3) g++ requires -ldl to be specified (last)
LIBS = -lm -lz -lpthread -L$(SRA_PATH)/lib64 \
	-lncbi-ngs-c++ -lncbi-ngs -lngs-c++ -lncbi-vdb-static -ldl
	
.SUFFIXES : .o .cpp .c
.cpp.o:
	$(CC) $(FLAGS) $(PROFILE) $(OPENMP) $(INC) -c $<
 
.c.o:
	$(cc) $(FLAGS) $(PROFILE) $(OPENMP) $(INC) -c $<

all: sracat 
	 	
sracat: $(OBJS) sracat.o 
	$(CC) $(PROFILE) -o sracat $(OBJS) $(OPENMP) sracat.o $(LIBS)

clean:
	-rm -f *.o


