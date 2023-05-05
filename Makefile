CC 		= gcc
CFLAGS 	+= -Wall -g -Wextra
INCLUDE += -I./lib
LDLIB	+= -L./lib -lmpi_subset -lnetcdf -lm

LIB 	= ./lib
TEST 	= ./test
DATA	= ./test/data_FEM

.PHONY: clean run-hello run-pi run-FEM

LIST_TEST = $(TEST)/test_hello $(TEST)/test_pi $(TEST)/FEM_MPI $(TEST)/FEM $(TEST)/FEM_MPI_REF

all: $(LIB)/libmpi_subset.a $(LIST_TEST)

$(LIB)/kahn.o: $(LIB)/kahn.c $(LIB)/kahn.h
	$(CC) $< $(CFLAGS) -c -o $@

$(LIB)/mpi_subset.o: $(LIB)/mpi_subset.c $(LIB)/mpi_subset.h
	$(CC) $< $(CFLAGS) -c -o $@

$(LIB)/libmpi_subset.a: $(LIB)/kahn.o $(LIB)/mpi_subset.o
	ar rcu $@ $+
	ranlib $@

$(TEST)/FEM: $(TEST)/FEM.o
	$(CC) $< $(CFLAGS) -o $@ -lm -lnetcdf

$(TEST)/FEM.o: $(TEST)/FEM.c
	$(CC) $< $(CFLAGS) -c -o $@

$(TEST)/FEM_MPI_REF: $(TEST)/FEM_MPI_REF.o
	$(CC) $< $(CFLAGS) -o $@ -L/home/sylvain/.local/lib -lm -lnetcdf -lmpi

$(TEST)/FEM_MPI_REF.o: $(TEST)/FEM_MPI_REF.c 
	$(CC) $< $(CFLAGS) -c -o $@ -I/home/sylvain/.local/include

$(TEST)/%: $(TEST)/%.o $(LIB)/libmpi_subset.a
	$(CC) $< $(CFLAGS) -o $@  $(LDLIB)

$(TEST)/%.o: $(TEST)/%.c 
	$(CC) $< $(CFLAGS) -c -o $@ $(INCLUDE)

run-hello:
	$(TEST)/test_hello -np 4

run-pi:
	$(TEST)/test_pi -np 4

run-FEM-mpi:
	$(TEST)/FEM_MPI -np 4
#	ncview $(DATA)/FEM_MPI.nc

run-FEM:
	$(TEST)/FEM
#	ncview $(DATA)/FEM.nc

run-FEM-mpi-ref:
	mpiexec -np 4 $(TEST)/FEM_MPI_REF

clean: 
	rm $(TEST)/*.o $(LIB)/*.o $(LIST_TEST) $(LIB)/*.a