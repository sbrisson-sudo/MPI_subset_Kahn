#ifndef MPI_INCLUDED
#define MPI_INCLUDED

/*      
        Interface implémentation Message Passing Interface
            - MPI_Datatype  // types mis à disposition pour les communications
            - MPI_Comm      // communicateurs (cette implémentation ne permet l'utilisation que d'un seul communicateur: MPI_Comm_World)
            - MPI_Init()    // initie l'environement MPI, lance les processus en parallèle et alloue les moyens de communication
            - MPI_Finalize()    // termine environement MPI
            - MPI_Comm_size // donne la taille d'un communicateur (nombre de processus)
            - MPI_Comm_rank // donne le rang du processus au sein du communicateur
            - MPI_Send()    // envoi de données
            - MPI_Receive() // reception données
            - MPI_Barrier() // synchronisation: appel bloquant tant que tous les processus n'ont pas appelé la fonction
*/

typedef enum { 
    MPI_BYTE,
    MPI_INT, 
    MPI_FLOAT,
    MPI_DOUBLE,
    MPI_LONG_DOUBLE
} MPI_Datatype;

typedef int MPI_Comm;
typedef int MPI_Status;

#define MPI_COMM_WORLD      0
#define MPI_SUCCESS         0

int MPI_Init(int argc, char **argv);
int MPI_Finalize(void);

int MPI_Comm_size(MPI_Comm comm, int *psize);
int MPI_Comm_rank(MPI_Comm comm, int *prank);

int MPI_Send(void *buf, int cnt, MPI_Datatype dtype, int dest, int tag, MPI_Comm comm);
int MPI_Recv(void *buf, int cnt, MPI_Datatype dtype, int src, int tag,
 MPI_Comm comm, MPI_Status *status);

int MPI_Barrier(MPI_Comm comm);

 #endif



