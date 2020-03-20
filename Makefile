OBJS = options.o ifind.o hash.o word.o mpi_util.o file_util.o \
	update.o bloom.o database.o binary_io.o parse_sequence.o

DOWNLOAD_OBJS = options.o ifind.o hash.o file_util.o parse_tar.o binary_io.o \
	download.o

CC = mpic++

PROFILE = #-pg
OPENMP = -Xpreprocessor -fopenmp
FLAGS = $(PROFILE) -O3 -Wall $(OPENMP) -std=c++11

INC = -I. -I$(HOME)/zlib/include \
	-I$(HOME)/src/BIGSI/SRA/include \
	-I$(HOME)/src/BIGSI/SRA/ncbi-vdb-master/interfaces \
	-I$(HOME)/llvm-project/build-openmp/runtime/src

LIBS = -lm $(HOME)/zlib/lib/libz.a \
	-L$(HOME)/llvm-project/build-openmp/runtime/src -lomp \
	-L$(HOME)/src/BIGSI/SRA/lib64 \
	-lncbi-ngs-c++       \
	-lngs-c++            \
	-lncbi-vdb-static    \
	-ldl


.SUFFIXES : .o .cpp .c
.cpp.o:
	$(CC) $(FLAGS) $(INC) -c $<

all: build_db bigsi++ sra_download bloomer checkbloom
	
build_db : $(OBJS) build_db.o
	$(CC) $(PROFILE) -o build_db $(OBJS) build_db.o $(LIBS) $(OPENMP)

bigsi++ : $(OBJS) bigsi++.o
	$(CC) $(PROFILE) -o bigsi++ $(OBJS) bigsi++.o $(LIBS) $(OPENMP)

bloomer : $(OBJS) bloomer.o
	$(CC) $(PROFILE) -o bloomer $(OBJS) bloomer.o $(LIBS) $(OPENMP)

checkbloom : bloom.o checkbloom.o file_util.o ifind.o binary_io.o
	$(CC) $(PROFILE) -o checkbloom bloom.o checkbloom.o file_util.o ifind.o binary_io.o $(LIBS) $(OPENMP)

sra_download: $(DOWNLOAD_OBJS) sra_download.o
	$(CC) $(PROFILE) -o sra_download $(DOWNLOAD_OBJS) sra_download.o $(LIBS) $(OPENMP)
	
clean:
	-rm -f *.o


