BUILD_OBJS = options.o ifind.o hash.o mpi_util.o file_util.o \
	update.o bloom.o database.o binary_io.o

BIGSI_OBJS = options.o ifind.o hash.o word.o mpi_util.o file_util.o \
	update.o bloom.o database.o binary_io.o parse_sequence.o

BLOOMER_OBJS = options.o ifind.o hash.o word.o mpi_util.o file_util.o \
	bloom.o binary_io.o
	
DOWNLOAD_OBJS = options.o ifind.o hash.o file_util.o parse_tar.o binary_io.o

CHECK_OBJS = bloom.o file_util.o ifind.o binary_io.o

CC = mpic++

# The profile flag is only needed for development
PROFILE = #-pg

OPENMP = -fopenmp
FLAGS = $(PROFILE) -O3 -Wall $(OPENMP) -std=c++0x

INC = -I. -I$(HOME)/zlib/include
	
LIBS = -lm $(HOME)/zlib/lib/libz.a

# The SRA libraries are *only* needed by the bloomer program to allow
# direct parsing of *.sra files
SRA_LIBS = -L$(HOME)/src/BIGSI/SRA/lib64 \
	-lncbi-ngs-c++       \
	-lngs-c++            \
	-lncbi-vdb-static    \
	-ldl
	
.SUFFIXES : .o .cpp
.cpp.o:
	$(CC) $(FLAGS) $(INC) -c $<

# A special rule to compile bloomer.cpp (which needs the SRA include files)
bloomer.o : bloomer.cpp
	$(CC) $(FLAGS) $(INC) -I$(HOME)/src/BIGSI/SRA/include -c $<
	
all: build_db bigsi++ sra_download bloomer checkbloom
	
build_db : $(BUILD_OBJS) build_db.o
	$(CC) $(PROFILE) -o build_db $(BUILD_OBJS) build_db.o $(LIBS) $(OPENMP)

bigsi++ : $(BIGSI_OBJS) bigsi++.o
	$(CC) $(PROFILE) -o bigsi++ $(BIGSI_OBJS) bigsi++.o $(LIBS) $(OPENMP)

# Note that the bloomer program requires the SRA libraries to enable the direct
# parsing of *.sra files
bloomer : $(BLOOMER_OBJS) bloomer.o
	$(CC) $(PROFILE) -o bloomer $(BLOOMER_OBJS) bloomer.o $(LIBS) $(SRA_LIBS) $(OPENMP)

checkbloom : $(CHECK_OBJS) checkbloom.o
	$(CC) $(PROFILE) -o checkbloom  $(CHECK_OBJS) checkbloom.o $(LIBS) $(OPENMP)
	
sra_download: $(DOWNLOAD_OBJS) sra_download.o
	$(CC) $(PROFILE) -o sra_download $(DOWNLOAD_OBJS) sra_download.o $(LIBS) $(OPENMP)

clean:
	-rm -f *.o
