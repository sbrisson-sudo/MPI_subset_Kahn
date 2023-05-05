#include <stdlib.h>
#include <stdio.h>
#include <sys/types.h>
#include <unistd.h>

#include "mpi_subset.h"

int main(int argc, char** argv){
    int rank, comm_size;

    MPI_Init(argc, argv);
    MPI_Comm_size(MPI_COMM_WORLD, &comm_size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    MPI_Status status;


    printf("[%d]\tHello from process %d of %d\tPID [%d]\n", rank, rank, comm_size, (int)getpid());

    printf("[%d]\tHitting barrier line %d\n", rank, __LINE__);
    MPI_Barrier(MPI_COMM_WORLD);

    float i[3];

    if (rank == 1) {    
        i[0] = 10.1; i[1] = 20.2; i[2] = 30.3; 
        for (int k=0; k<5; k++) {
            i[0] += 5; i[1] += 5; i[2] += 5;
            MPI_Send(&i, 3, MPI_INT, 0, 0, MPI_COMM_WORLD);
            printf("[%d]\tSending\t\t[%f,%f,%f]\n", rank, *i,i[1],i[2]); fflush(stdout);
            MPI_Recv(&i, 3, MPI_INT, 0, 0, MPI_COMM_WORLD, &status);
            printf("[%d]\tReceiving\t[%f,%f,%f]\n", rank, *i,i[1],i[2]); fflush(stdout);

        }
    }

    if (rank == 0) {    
        i[0] = 70.1; i[1] = 80.2; i[2] = 90.3; 
        for (int k=0; k<5; k++) {
            i[0] += 5; i[1] += 5; i[2] += 5;
            MPI_Send(&i, 3, MPI_INT, 1, 0, MPI_COMM_WORLD);
            printf("[%d]\tSending\t\t[%f,%f,%f]\n", rank, *i,i[1],i[2]); fflush(stdout);
            MPI_Recv(&i, 3, MPI_INT, 1, 0, MPI_COMM_WORLD, &status);
            printf("[%d]\tReceiving\t[%f,%f,%f]\n", rank, *i,i[1],i[2]); fflush(stdout);
        }
    }


    printf("[%d]\tGood bye from process %d of %d\tPID [%d]\n", rank, rank, comm_size, (int)getpid());
    
    MPI_Finalize();
    return 0;
}



