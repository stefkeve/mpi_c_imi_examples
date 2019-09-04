#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const long A = 1000000000;

int main(int argc, char *argv[])
{
  int i, j;
  int id; //Process rank
  int np; //Num of processes
  double x, y, PI; //coordinates
  long in, inI; // number of points inside the circle
  double elapsed_time;
  double PI25DT = 3.141592653589793238462643;
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Comm_size(MPI_COMM_WORLD, &np);

  if(id == 0)
    elapsed_time = - MPI_Wtime();

  srand ( time(NULL)*id );
  in = 0;

  for(i = 0; i < A/np; i++)
  {
    x = ((double)rand() / RAND_MAX);
    y = ((double)rand() / RAND_MAX);

    if (x*x + y*y <= 1) in++; 		
  }

  if(id != 0)
      MPI_Send(&in, 1, MPI_INT, 0, 20, MPI_COMM_WORLD);
  else  
  {
    for(i = 1; i < np; i++)
    {
      MPI_Recv(&inI, 1, MPI_INT, i, 20, MPI_COMM_WORLD, &status);
      in += inI;
    }
  }

  if(id == 0)
  {
    elapsed_time += MPI_Wtime();
    PI = (double)4*in/A;
    printf("PI = %26.25f,  and elapsed time is %f\n",PI, elapsed_time);
    printf("PI25 - PIour = %26.25f\n",PI25DT - PI);
  }

  MPI_Finalize();

  return 0;
}

