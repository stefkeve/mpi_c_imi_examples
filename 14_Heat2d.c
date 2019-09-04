#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#define INNER_START(rank) (rank == 0 ? 2 : 1)
#define INNER_UPPER_BOUND(rank, numProc, npX) (rank == numProc - 1 ? npX - 1 : npX)

int main(int argc, char **argv) {
    int id;
    int numProc;
    int tag;
    double wallTime;

    int ncellsX = 480;
    int ncellsY = 480;
    int nSteps  = 500;
    int lengthX = 180;
    int lengthY = 180;
    double ao   = 1.0;
    double coeff = 0.01;
    double sigma = 3.0;

    double dx    = lengthX / ncellsX;
    double dy    = lengthY / ncellsY;

    int npX;
    int innerStartX;
    int innerUpperBoundX;
    double **temp, **tempNew;

    MPI_Status status;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &numProc);
    
    if(id == 0) {
	    npX = ncellsX / numProc;
    }

    MPI_Bcast(&npX, 1, MPI_INT, 0, MPI_COMM_WORLD);

    innerStartX = INNER_START(id);
    innerUpperBoundX = INNER_UPPER_BOUND(id, numProc, npX);

    temp = (double **) malloc((npX + 2) * sizeof(double*));
    tempNew = (double **) malloc((npX + 2) * sizeof(double*));

    for(int i = 0; i <= npX + 1; ++i) {
        temp[i]    = (double *) malloc(ncellsY * sizeof(double));
        tempNew[i] = (double *) malloc(ncellsY * sizeof(double));
    }

   if(id == 0) {
        wallTime = -MPI_Wtime();
    }


    for(int j = 0; j < ncellsY; ++j) {
        temp[0][j] = temp[npX + 1][j] = 0.0;
    }

    for(int i = 1; i <= npX; ++i) {
        for(int j = 0; j < ncellsY; ++j) {
            double x = (i-1 + id * npX) * dx;
            double y = j * dy;
            temp[i][j] = ao*exp(-x*x/(2.0*sigma*sigma)) + ao*exp(-y*y/(2.0*sigma*sigma));
        }
    }

    for(int step = 0; step < nSteps; ++step) {
       // communication part
       if(id > 0) {
         MPI_Send(&temp[1][0], ncellsY, MPI_DOUBLE, id-1, 1, MPI_COMM_WORLD);
       }

       if(id < numProc - 1) {
         MPI_Recv(&temp[npX+1][0], ncellsY, MPI_DOUBLE, id+1, 1, MPI_COMM_WORLD, &status);
       }

       if(id < numProc - 1) {
         MPI_Send(&temp[npX][0], ncellsY, MPI_DOUBLE, id+1, 2, MPI_COMM_WORLD);
       }

       if(id > 0) {
         MPI_Recv(&temp[0][0], ncellsY, MPI_DOUBLE, id-1, 2, MPI_COMM_WORLD, &status);
       }

       for(int i = innerStartX; i <= innerUpperBoundX; ++i) {
           for(int j = 1; j < ncellsY-1; ++j) {
               tempNew[i][j] = temp[i][j] + coeff*(temp[i-1][j] + temp[i+1][j] - 4.0*temp[i][j] +
                                                   temp[i][j-1] + temp[i][j+1]);
           }
       }

       for(int ii = innerStartX; ii <= innerUpperBoundX; ++ii) {
           for(int jj = 1; jj < ncellsY-1; ++jj) {
               temp[ii][jj] = tempNew[ii][jj];
           }
       }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    if(id == 0) {
	    wallTime += MPI_Wtime();
        printf("Walltime clock is %lf\n", wallTime);
    }

    // sync output
    for (int procId = 0; procId < numProc; ++procId) {
        if(procId == id) {
            for(int i = 1; i <= npX; ++i) {
                for(int j = 0; j < ncellsY; ++j) {
                    printf("%.5lf ", temp[i][j]);
                }
                printf("\n");
            }
        }

        MPI_Barrier(MPI_COMM_WORLD);
    }

    MPI_Finalize();
}
