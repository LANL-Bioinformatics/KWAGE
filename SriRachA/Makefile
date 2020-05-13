OBJS = date.o mpi_util.o options.o parse_sequence.o file_util.o ifind.o word.o \
	sra_stream.o search_by_kmer.o

CC = mpic++

PROFILE = #-pg
OPENMP = #-Xpreprocessor -fopenmp
FLAGS = $(PROFILE) -O3 -Wall $(OPENMP) -std=c++11

SRA=$(HOME)/src/BIGSI/SRA

INC = -I. \
	-I$(SRA)/include \
	-I$(SRA)/ncbi-vdb/interfaces

LIBS = -lm  \
	-L$(SRA)/lib64 \
	-lncbi-ngs-c++       \
	-lncbi-vdb-static    \
	-ldl


.SUFFIXES : .o .cpp .c
.cpp.o:
	$(CC) $(FLAGS) $(INC) -c $<

all: sriracha
	
sriracha : $(OBJS) main.o
	$(CC) $(PROFILE) -o sriracha $(OBJS) main.o $(LIBS) $(OPENMP)

clean:
	-rm -f *.o


