/*
 * Int.c
 * Izracunavanje odredjenog integrala date funkcije f(x) pomocu Simpsonovog pravila
 * Created on: Aug 20, 2010
 * Author: djordje
 */

#include "mpi.h"
#include <stdio.h>
#include <math.h>

const int A = 0;
const int B = 1;
double f( double a )
{
    return (4.0 / (1.0 + a*a));
}

int main( int argc, char *argv[])
{
  long n_intervals;
  long  i;
  int  myid, numprocs;
  double PI25DT = 3.141592653589793238462643;
  double mypi, pi, h, sum, sumI;
  double elapsed_time;
  MPI_Status status;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);

  if (myid == 0)
  {
    scanf("%ld", &n_intervals);
    elapsed_time = -MPI_Wtime();
  }
   
  if(myid == 0)
    for(i = 1; i < numprocs; i++)
      MPI_Send(&n_intervals, 1, MPI_LONG, i, 10, MPI_COMM_WORLD);
  else  
    MPI_Recv(&n_intervals, 1, MPI_LONG, 0, 10, MPI_COMM_WORLD, &status);

  h = 1.0 / n_intervals;
  sum = 0.0;

  for (i = myid + 1; i <= n_intervals/2; i += numprocs)
  {
      sum += 4 * f((2*i - 1) * h) + 2 * f((2*i) * h);
  }

  if(myid != 0)
      MPI_Send(&sum, 1, MPI_DOUBLE, 0, 20, MPI_COMM_WORLD);
  else  
  {
    for(i = 1; i < numprocs; i++)
    {
      MPI_Recv(&sumI, 1, MPI_DOUBLE, i, 20, MPI_COMM_WORLD, &status);
      sum += sumI;
    }
  }

  if (myid == 0)
  {
    elapsed_time += MPI_Wtime();
    mypi = (f(A) - f(B) + sum) / 3.0 / n_intervals;
    printf("Int[%d, %d](f(x)) = %.16f, Error is %.16f\n", A, B, mypi, fabs(mypi - PI25DT));
    printf("wall clock time = %f\n", elapsed_time);
  }

  MPI_Finalize();

  return 0;
}