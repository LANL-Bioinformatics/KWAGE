OBJS = maestro_main.o worker_main.o mpi_util.o options.o string_conversion.o \
	file_util.o hash.o ifind.o binary_io.o sra_accession.o build_db.o \
	make_bloom.o word.o bloom.o date.o mem_usage.o file_io.o sra_meta.o

INVENTORY_OBJS = options.o ifind.o hash.o file_util.o parse_tar.o binary_io.o \
	split.o date.o string_conversion.o sra_accession.o bloom.o

SEARCH_OBJS = parse_sequence.o bloom.o word.o options.o ifind.o hash.o \
	file_util.o binary_io.o date.o string_conversion.o sra_accession.o

DATABASE_OBJS = file_io.o binary_io.o sra_accession.o file_util.o ifind.o

CC = mpic++

PROFILE = #-g
OPENMP = #-Xpreprocessor -fopenmp
FLAGS = $(PROFILE) -O3 -Wall $(OPENMP) -mavx2 -std=c++11
 
# Edit the SRA_LIB_PATH to point to the SRA toolkit libraries on your system
SRA_LIB_PATH = $(HOME)/src/KWAGE/SRA/lib64

# Edit the SRA_INCLUDE_PATH to point to the SRA toolkit header files on your system
SRA_INCLUDE_PATH = $(HOME)/src/KWAGE/SRA/include

# Edit the NCBI_VDB_INCLUDE_PATH to point to the NCBI-VDB header files on your system
NCBI_VDB_INCLUDE_PATH = $(HOME)/src/KWAGE/SRA/ncbi-vdb/interfaces

OMP_LIBS = #-L$(HOME)/llvm-project/build-openmp/runtime/src -lomp
OMP_INCLUDE = #-I$(HOME)/llvm-project/build-openmp/runtime/src

INC = -I. \
	-I$(SRA_INCLUDE_PATH) \
	-I$(NCBI_VDB_INCLUDE_PATH) \
	$(OMP_INCLUDE)

LIBS = -lm -lz \
	-L$(SRA_LIB_PATH) \
	-lncbi-ngs-c++       \
	-lngs-c++            \
	-lncbi-vdb-static    \
	-ldl

.SUFFIXES : .o .cpp .c
.cpp.o:
	$(CC) $(FLAGS) $(INC) -c $<

all: sra_inventory maestro sra_diff inventory_dump dump_db dump_bloom kwage sra_dump \
	bloom_test bff bloom_diff manual_db merge_db
	
maestro : $(OBJS) maestro.o
	$(CC) $(PROFILE) -o maestro $(OBJS) maestro.o $(LIBS) $(OMP_LIBS) $(OPENMP)

bff : $(OBJS) bff.o
	$(CC) $(PROFILE) -o bff $(OBJS) bff.o $(LIBS) $(OMP_LIBS) $(OPENMP)

kwage : $(SEARCH_OBJS) kwage.o
	$(CC) $(PROFILE) -o kwage $(SEARCH_OBJS) kwage.o $(LIBS) $(OMP_LIBS) $(OPENMP)

sra_inventory: $(INVENTORY_OBJS) sra_inventory.o
	$(CC) $(PROFILE) -o sra_inventory $(INVENTORY_OBJS) sra_inventory.o $(LIBS) $(OMP_LIBS) $(OPENMP)

inventory_dump: $(INVENTORY_OBJS) inventory_dump.o
	$(CC) $(PROFILE) -o inventory_dump $(INVENTORY_OBJS) inventory_dump.o $(LIBS) $(OMP_LIBS) $(OPENMP)

sra_dump: sra_dump.o mem_usage.o string_conversion.o
	$(CC) $(PROFILE) -o sra_dump sra_dump.o mem_usage.o string_conversion.o $(LIBS) $(OMP_LIBS) $(OPENMP)

bloom_test: $(SEARCH_OBJS) bloom_test.o
	$(CC) $(PROFILE) -o bloom_test $(SEARCH_OBJS) bloom_test.o $(LIBS) $(OMP_LIBS) $(OPENMP)

dump_db: $(INVENTORY_OBJS) dump_db.o
	$(CC) $(PROFILE) -o dump_db $(INVENTORY_OBJS) dump_db.o $(LIBS) $(OMP_LIBS) $(OPENMP)

merge_db: $(DATABASE_OBJS) merge_db.o
	$(CC) $(PROFILE) -o merge_db $(DATABASE_OBJS) merge_db.o $(LIBS) $(OMP_LIBS) $(OPENMP)

manual_db: $(DATABASE_OBJS) manual_db.o
	$(CC) $(PROFILE) -o manual_db $(DATABASE_OBJS) manual_db.o $(LIBS) $(OMP_LIBS) $(OPENMP)

db_debug: $(INVENTORY_OBJS) db_debug.o build_db.o
	$(CC) $(PROFILE) -o db_debug $(INVENTORY_OBJS) db_debug.o build_db.o $(LIBS) $(OMP_LIBS) $(OPENMP)

dump_bloom: $(INVENTORY_OBJS) dump_bloom.o
	$(CC) $(PROFILE) -o dump_bloom $(INVENTORY_OBJS) dump_bloom.o $(LIBS) $(OMP_LIBS) $(OPENMP)

sra_diff: $(INVENTORY_OBJS) sra_diff.o
	$(CC) $(PROFILE) -o sra_diff $(INVENTORY_OBJS) sra_diff.o $(LIBS) $(OMP_LIBS) $(OPENMP)

bloom_diff: $(INVENTORY_OBJS) bloom_diff.o
	$(CC) $(PROFILE) -o bloom_diff $(INVENTORY_OBJS) bloom_diff.o $(LIBS) $(OMP_LIBS) $(OPENMP)

clean:
	-rm -f *.o


