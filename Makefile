BUILD_OBJS = options.o ifind.o hash.o mpi_util.o file_util.o \
	update.o bloom.o database.o binary_io.o date.o

BIGSI_OBJS = options.o ifind.o hash.o word.o mpi_util.o file_util.o \
	update.o bloom.o database.o binary_io.o parse_sequence.o date.o

BLOOMER_OBJS = options.o ifind.o hash.o word.o mpi_util.o file_util.o \
	bloom.o binary_io.o date.o
	
DOWNLOAD_OBJS = options.o ifind.o hash.o file_util.o parse_tar.o binary_io.o \
	date.o split.o

CHECK_OBJS = bloom.o file_util.o ifind.o binary_io.o date.o

DUMP_OBJS = binary_io.o date.o

CC = mpic++

# The profile flag is only needed for development
PROFILE = #-pg

OPENMP = #-fopenmp
FLAGS = $(PROFILE) -O3 -Wall $(OPENMP) -std=c++0x

INC = -I. -I$(HOME)/zlib/include
	
LIBS = -lm $(HOME)/zlib/lib/libz.a

# Edit the SRA_LIB_PATH variable to point to the directory on your system that contains
# the SRA library files
SRA_LIB_PATH = $(HOME)/src/BIGSI/SRA/lib64

# Edit the SRA_INCLUDE_PATH variable to point to the directory on your system that contains
# the SRA include files
SRA_INCLUDE_PATH = $(HOME)/src/BIGSI/SRA/include

# The SRA libraries are *only* needed by the bloomer program to allow
# direct parsing of *.sra files
SRA_LIBS = -L$(SRA_LIB_PATH) \
	-lncbi-ngs-c++       \
	-lngs-c++            \
	-lncbi-vdb-static    \
	-ldl
	
.SUFFIXES : .o .cpp
.cpp.o:
	$(CC) $(FLAGS) $(INC) -c $<
	
all: build_db bigsi++ sra_download bloomer checkbloom dump_db

build_db : $(BUILD_OBJS) build_db.o
	$(CC) $(PROFILE) -o build_db $(BUILD_OBJS) build_db.o $(LIBS) $(OPENMP)

bigsi++ : $(BIGSI_OBJS) bigsi++.o
	$(CC) $(PROFILE) -o bigsi++ $(BIGSI_OBJS) bigsi++.o $(LIBS) $(OPENMP)

# A special rule to compile bloomer.cpp (which needs the SRA include files)
bloomer.o : bloomer.cpp
	$(CC) $(FLAGS) $(INC) -I$(SRA_INCLUDE_PATH) -c $<

# Note that the bloomer program requires the SRA libraries to enable the direct
# parsing of *.sra files
bloomer : $(BLOOMER_OBJS) bloomer.o
	$(CC) $(PROFILE) -o bloomer $(BLOOMER_OBJS) bloomer.o $(LIBS) $(SRA_LIBS) $(OPENMP)

checkbloom : $(CHECK_OBJS) checkbloom.o
	$(CC) $(PROFILE) -o checkbloom  $(CHECK_OBJS) checkbloom.o $(LIBS) $(OPENMP)
	
sra_download: $(DOWNLOAD_OBJS) sra_download.o
	$(CC) $(PROFILE) -o sra_download $(DOWNLOAD_OBJS) sra_download.o $(LIBS) $(OPENMP)

dump_db: $(DUMP_OBJS) dump_db.o
	$(CC) $(PROFILE) -o dump_db $(DUMP_OBJS) dump_db.o $(LIBS) $(OPENMP)

clean:
	-rm -f *.o
