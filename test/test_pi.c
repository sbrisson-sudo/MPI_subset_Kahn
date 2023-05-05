#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#include "mpi_subset.h"

/*      Programme calculant une valeur approchée de pi (formule de Madhava-Leibniz) de manière parallèle.
        Options:
            --nb-iter   // nombre d'itérations de la somme à réaliser
            --quiet     // le programme n'affiche rien sur stdout
        
*/


int main(int argc, char **argv)
{   

    int rank, size, istart, istop;

    MPI_Status status;
    MPI_Comm comm;
    comm = MPI_COMM_WORLD;

    double pi, receive_pi;

    MPI_Init(argc, argv);

    int N       = 1000000;
    int output  = 1;

    while(argv++,--argc) {
        if (strcmp(*argv, "--nb-iter") == 0) {
            if (argc > 1) {
                N = (int)strtol((const char *)*(++argv), NULL, 10);
                argc--;
            }
        }
        if (strcmp(*argv, "--quiet") == 0) {
            output = !output;
        }
    }

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if ((rank == 0) && output) {
        printf("[global] N = %d\n", N);
    }

    istart = N/size * rank + 1;
    istop = istart + N/size -1;

    if (output) { printf("[%d] istart = %d ; istop = %d\n", rank, istart, istop);}

    double local_bit_of_pi = 0.0;
    for (int i=istart; i<= istop; i++)
    {
        local_bit_of_pi += 1.0 / ( 1.0 + ( (i-0.5) / N)*((i-0.5) / N) );
    }

    if (output) { printf("[%d] local_bit_of_pi = %f\n", rank, local_bit_of_pi); }

    if (rank == 0)
    {
        pi = local_bit_of_pi;
        for (int source = 1; source < size; source++)
        { 
            int tag = 0;
            MPI_Recv(&receive_pi, 1, MPI_DOUBLE, source, tag, comm, &status);

            if (output) { printf("[0] receiving %f from %d\n", receive_pi, source); }

            pi += receive_pi;
        }
        
        pi *= 4.0/(long double)N;

        if (output) { printf("\n[0] pi = %.15f\n", pi); }

    }

    else 
    {
        int tag = 0;

        if (output) { printf("[%d] sending %f to 0\n", rank, local_bit_of_pi); }

        MPI_Send(&local_bit_of_pi, 1, MPI_DOUBLE, 0, tag, comm);
    }


    MPI_Finalize();

    if (output) { printf("[%d] Ended\n", rank); }

    return 0;

}

