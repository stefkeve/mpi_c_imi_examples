#include <stdio.h>
#include <stdlib.h>
#include <mpi.h>
#include <math.h>

#ifndef M_PI
    #define M_PI 3.14159265358979323846
#endif

int main(int argc, char **argv) {
    int id;
    int numProc;
    int tag;
    double wallTime;

    int nX         = 80;
    int nSteps     = 300;
    int npX;
    double lenghtX = 1.0;
    double alpha   = 0.005;
    double tMax    = 1.0;
    double dx = lenghtX / (nX-1);
    double dt = tMax / (nSteps-1);
    double r  = alpha * dt/pow(dx,2);
    double r2 = 1 - 2*r;

    double **temp;

    MPI_Status status;

    MPI_Init (&argc, &argv);
    MPI_Comm_rank (MPI_COMM_WORLD, &id);
    MPI_Comm_size (MPI_COMM_WORLD, &numProc);

    if(id == 0) {
	    npX = nX / numProc; 
    }
  
    MPI_Bcast(&npX, 1, MPI_INT, 0, MPI_COMM_WORLD); 

    temp = (double **) malloc((nSteps+1) * sizeof(double*));
    for(int i = 0; i <= nSteps; ++i) {
        temp[i] = (double *) malloc((npX + 2) * sizeof(double));
        temp[i][0] = temp[i][npX + 1] = 0.0;
    }   
 
    if(id == 0) {
        wallTime = -MPI_Wtime();
    }
 
    for(int i = 1; i <= npX; ++i) {
       temp[0][i] = sin(M_PI * (i-1 + id * npX)*dx);
    }

    for(int step = 1; step <= nSteps; ++step) {
       // communication part
       if(id > 0) {
         MPI_Send(&temp[step-1][1], 1, MPI_DOUBLE, id-1, 1, MPI_COMM_WORLD);
       }

       if(id < numProc - 1) {
         MPI_Recv(&temp[step-1][npX+1], 1, MPI_DOUBLE, id+1, 1, MPI_COMM_WORLD, &status);
       }

       if(id < numProc - 1) {
         MPI_Send(&temp[step-1][npX], 1, MPI_DOUBLE, id+1, 2, MPI_COMM_WORLD);
       }
              
       if(id > 0) {
         MPI_Recv(&temp[step-1][0], 1, MPI_DOUBLE, id-1, 2, MPI_COMM_WORLD, &status);
       }

       for(int j = 1; j <= npX; ++j) {
           temp[step][j] = r*temp[step-1][j-1] - r2*temp[step-1][j] + r*temp[step-1][j+1];
       }

       if(id == 0) {
           temp[step][1] = 0;
       }
       
       if(id == numProc - 1) {
           temp[step][npX] = 0;
       }
    }        
 
    if(id == 0) {
            wallTime += MPI_Wtime();
        printf("Walltime clock is %lf\n", wallTime);
    }
   
    for(int i = 1; i <= npX; ++i) {
        printf("temp[%d] = %lf\n", (i-1 + id * npX), temp[nSteps][i]);
    }

    MPI_Finalize();
}
