OBJS = 
	
CC = c++
cc = cc

FLAGS = -O3 -Wall

PROFILE = #-pg
OPENMP = #-fopenmp

# Edit the SRA_INCLUDE_PATH to point to the SRA toolkit on your system
SRA_PATH = ../SRA

INC = -I. -I$(SRA_PATH)/include -I$(SRA_PATH)/ngs/ngs-sdk
LIBS = -lm -lz -L$(SRA_PATH)/lib64 -L$(SRA_PATH)/ngs/ngs-sdk/lib64 \
	-lncbi-ngs -lngs-c++-static -lncbi-ngs-c++-static -lncbi-vdb-static -ldl

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


