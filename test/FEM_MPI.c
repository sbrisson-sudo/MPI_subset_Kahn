#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <assert.h>
#include <netcdf.h>
#include <string.h>
#include <sys/time.h>   // bench-marking


#include <mpi_subset.h>
//#include <mpi.h>

/*      FEM_MPI.c

    Script de résolution de l'équation de d'Alembert pour des ondes sismiques en milieu non dispersif non absorbant homogène, dérivé en 3 variantes:
                - FEM.c         // script de base, calcul non parallèle
                - FEM_MPI.c     // script parallélisé appelant mpi_subset.h
                - FEM_MPI_REF.c // script parraléliséa appelant mpi.h (nécésitte d'avoi une implémentation de MPI installée)

    a parallélisation se fait sur 4 coeurs avec uniquemnt communications point-à-point.

    Ces 3 scripts écrivent leurs résultats sous format netCDF, pas d'implémentation d'écriture parallèle, les scripts parallélisés n'écrivent qu'un quart des données
    Les données peuvent être visualisées avec l'utilitaire ncview.

    Les variables de préprocesseurs DEBUG et EXPLICIT_ECHANGE permettent de controler le niveau de verbose.
    La variable de préprocesseur TIMER permet de mesurer le temps de calcul.
*/


#define N_SPACE     300
#define N_TIME      1000

#define N_DIMS      3

#define FILE_NAME   "./test/data_FEM/FEM_MPI.nc"

#define ERR(e) {printf("Error: %s\n", nc_strerror(e)); return 2;}


//#define DEBUG
//#define EXPLICIT_ECHANGE
#define TIMER


// int main()
int main(int argc, char** argv)
{
    // Initiation MPI

    int rank, size;
    MPI_Status status;

    //MPI_Init(NULL, NULL);
    MPI_Init(argc, argv);

    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    MPI_Comm_size(MPI_COMM_WORLD, &size);

    if (size != 4) {
        fprintf(stderr, "Error, this algorithm is made to run on at least 4 processes\n"); exit(1);
    }

    const int H      = 3000;         // dimension spatiale
    const float mu   = 40e9;         // module de cisaillement
    const float rho  = 3000;         // masse volumique
    const float c    = sqrt(mu/rho); // célérité ondes compression  

    const float eps  = 0.5;            // critère de stabilité de CFL
    const float dx   = (float)H / (float)N_SPACE;
    const float dy   = (float)H / (float)N_SPACE;
    const float dt   = dx * eps / c;

    const float T    = dt * N_TIME;   // durée de la simulation

    if (rank == 0) {
        printf("Discrétisation simulation:\n\tN_space\t=\t%d\n\tN_time\t=\t%d\n", N_SPACE, N_TIME);
        printf("Paramètres de la simulation:\n\tc\t=\t%g m/s\n\tdx\t=\t%g m\n\tdt\t=\t%g s\n",c, dx, dt);
        printf("Durée de la simulation\t%g s\n", T);
        printf("Nombre de processus MPI\t%d\n", size);
        printf("Fichier écriture données\t%s\n\n", FILE_NAME);
    }

    float d2_field_x, d2_field_y;

    float time[N_TIME], source[N_TIME], space_x[N_SPACE], space_y[N_SPACE];

    // Initiation vecteurs temps, espace
    time[0] = 0.f;
    for (int i =1; i<N_TIME; i++) {
        time[i] = time[i-1] + dt;
    }
    space_x[0] = 0.f;
    space_y[0] = 0.f;
    for (int i =1; i<N_SPACE; i++) {
        space_x[i] = space_x[i-1] + dx;
        space_y[i] = space_y[i-1] + dy;
    }

    // Allocation buffer données

    int N_SPACE_LOCAL = N_SPACE / 2 + 1;

    float field_old[N_SPACE_LOCAL][N_SPACE_LOCAL],
            field[N_SPACE_LOCAL][N_SPACE_LOCAL],
            field_new[N_SPACE_LOCAL][N_SPACE_LOCAL];

    for (int i=0; i<N_SPACE_LOCAL; i++) {
        for (int j=0; j<N_SPACE_LOCAL; j++) {
            field_old[i][j] = 0.f;
            field_new[i][j] = 0.f;
            field[i][j] = 0.f;
        }
    }

    // Terme de source

    int src_idx_x     = N_SPACE_LOCAL / 2;
    int src_idx_y     = N_SPACE_LOCAL / 2;

    if (rank == 0) {
        float freq     = 40.f;
        float t0       = 4 / freq;

        for (int i =0; i<N_TIME; i++) {
            source[i] = -2*(time[i] - t0)*freq*freq * \
            exp(-1*(time[i]-t0)*(time[i]-t0)*freq*freq);
        }
    }

    // Initiation fichier netCDF pour écriture

    int nc_id, x_dimid, y_dimid, time_dimid, dim_ids[N_DIMS];
    int field_varid, time_varid, x_varid, y_varid; 
    size_t start[N_DIMS], count[N_DIMS];
    int retval;
    // MPI_Info info=MPI_INFO_NULL;

    if (rank == 0) {

    if ((retval = nc_create(FILE_NAME, NC_CLOBBER, &nc_id))){ERR(retval);}

    if ((retval = nc_def_dim(nc_id, "x", N_SPACE, &x_dimid))){ERR(retval);}
    if ((retval = nc_def_dim(nc_id, "y", N_SPACE, &y_dimid))){ERR(retval);}
    if ((retval = nc_def_dim(nc_id, "time_dim", N_TIME, &time_dimid))){ERR(retval);}

    if ((retval = nc_def_var(nc_id, "space_x", NC_FLOAT, 1, &x_dimid, &x_varid))){ERR(retval);}
    if ((retval = nc_def_var(nc_id, "space_y", NC_FLOAT, 1, &y_dimid, &y_varid))){ERR(retval);}
    if ((retval = nc_def_var(nc_id, "time", NC_FLOAT, 1, &time_dimid, &time_varid))){ERR(retval);}

    dim_ids[0]  = time_dimid; 
    dim_ids[1]  = x_dimid;
    dim_ids[2]  = y_dimid;

    if ((retval = nc_def_var(nc_id, "u", NC_FLOAT, N_DIMS, dim_ids, &field_varid))){ERR(retval);}

    if ((retval = nc_enddef(nc_id))){ERR(retval);}

    if ((retval = nc_put_var_float(nc_id, x_varid, &space_x[0]))){ERR(retval);}
    if ((retval = nc_put_var_float(nc_id, y_varid, &space_y[0]))){ERR(retval);}
    if ((retval = nc_put_var_float(nc_id, time_varid, &time[0]))){ERR(retval);}

    start[1]    = 0;
    start[2]    = 0;

    count[1]    = N_SPACE_LOCAL;
    count[2]    = N_SPACE_LOCAL;
    count[0]    = 1;
    }

    // buffer échanges colonnes

    float buf[N_SPACE_LOCAL];
    int side1 = (rank+1)%size, side2 = (rank-1)%size;
    if (side2 < 0) { side2 = side2 + size;}

    //MPI_Barrier(MPI_COMM_WORLD);

#ifdef DEBUG
    printf("[%d] Communication: side1 = %d ; side2 = %d\n", rank, side1, side2);
#endif

#ifdef TIMER
    // Bench-marking
    struct timeval t1, t2;
    double elapsedTime;
    gettimeofday(&t1, NULL);
#endif

    for (int t=0; t<N_TIME; t++) {

#ifdef DEBUG
        printf("[%d] Entrée itération [%d]\n", rank, t);
#endif

        // Résolution approchée équation de D'Alembert

        for (int x=1; x < N_SPACE_LOCAL; x++) {
            for (int y=1; y < N_SPACE-1; y++) {
                d2_field_x = (field[x+1][y] + field[x-1][y] - 2*field[x][y]) / (dx*dx);
                d2_field_y = (field[x][y+1] + field[x][y-1] - 2*field[x][y]) / (dy*dy);

                field_new[x][y] = 2*field[x][y] - field_old[x][y] + dt*dt * c*c *(d2_field_x + d2_field_y);
            }
        }

        // Conditions aux limites (2 surfaces libres + 2 échanges)

        // Surfaces libres
        for (int x=1; x < N_SPACE_LOCAL-1; x ++) {
            field_new[0][x] = field_new[1][x];
            field_new[x][0] = field_new[x][1];
        }

        // Echanges
        MPI_Send(field_new[N_SPACE_LOCAL-2], N_SPACE_LOCAL, MPI_FLOAT,side1 , 0, MPI_COMM_WORLD);
        for (int i=0; i<N_SPACE_LOCAL; i++) buf[i] = field_new[i][N_SPACE_LOCAL-2];
        MPI_Send(buf, N_SPACE_LOCAL, MPI_FLOAT, side2, 0, MPI_COMM_WORLD);

#ifdef EXPLICIT_ECHANGE
        printf("[%d][i=%d]\tSending\n\t[%.6f,...,%.6f,...,%.6f] to %d\n\t[%.6f,...,%.6f,...,%.6f] to %d\n",
        rank, t, field_new[N_SPACE_LOCAL-2][0], field_new[N_SPACE_LOCAL-2][N_SPACE_LOCAL/2], field_new[N_SPACE_LOCAL-2][N_SPACE_LOCAL-1], side1, buf[0], buf[N_SPACE_LOCAL/2], buf[N_SPACE_LOCAL-1], side2);
        fflush(stdout);
#endif

        MPI_Recv(field_new[N_SPACE_LOCAL-1], N_SPACE_LOCAL, MPI_FLOAT, side1, 0, MPI_COMM_WORLD, &status);
        MPI_Recv(buf, N_SPACE_LOCAL, MPI_FLOAT, side2, 0, MPI_COMM_WORLD, &status);
        for (int i=0; i<N_SPACE_LOCAL; i++) field_new[i][N_SPACE_LOCAL-1] = buf[i];

#ifdef EXPLICIT_ECHANGE
        printf("[%d][i=%d]\tReceiving\n\t[%.6f,...,%.6f,...%.6f] from %d\n\t[%.6f,...,%.6f,...%.6f] from %d\n",
        rank, t, field_new[N_SPACE_LOCAL-1][0], field_new[N_SPACE_LOCAL-1][N_SPACE_LOCAL/2], field_new[N_SPACE_LOCAL-1][N_SPACE_LOCAL-1], side1, buf[0], buf[N_SPACE_LOCAL/2], buf[N_SPACE_LOCAL-1], side2);
        fflush(stdout);

#endif

        // Terme de source
        if (rank == 0){
            field_new[src_idx_x][src_idx_y] = field_new[src_idx_x][src_idx_y] + source[t];
#ifdef DEBUG
            printf("[%d] u_src = %g\n", rank, field_new[src_idx_x][src_idx_y]);
#endif
        }

        if (rank == 0) {
        // Ecriture
        start[0] = t;
        nc_put_vara_float(nc_id, field_varid, start, count, &field[0][0]);
        }

        // Avancement d'un pas de temps
        for (int x=0; x < N_SPACE_LOCAL; x++) {
            for (int y=0; y < N_SPACE_LOCAL; y++) {
                field_old[x][y] = field[x][y];
                field[x][y] = field_new[x][y];
            }
        }
    }

    MPI_Finalize();

#ifdef TIMER
    if (rank == 0) {
        gettimeofday(&t2, NULL);
        elapsedTime = (t2.tv_sec - t1.tv_sec) * 1000.0;      // sec to ms
        elapsedTime += (t2.tv_usec - t1.tv_usec) / 1000.0;   // us to ms
        printf("Temps total (hors initiation) %f ms.\n", elapsedTime);
    }   
#endif

    return 0;
}