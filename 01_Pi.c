/*
 * PI.c
 * Izracunavanje broja PI pomocu odredjenog integrala(od 0 do 1) funkcije 4.0 / (1.0 + x*x)
 *  Created on: Aug 20, 2010
 *      Author: djordje
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double f( double );
double f( double a )
{
    return (4.0 / (1.0 + a*a));
}

int main( int argc, char *argv[])
{
  long n_intervals;// = 100000000;
  long n;
  int myid, numprocs, i;
  double PI25DT = 3.141592653589793238462643;
  double mypi, pi, h, sum, x;
  double startwtime = 0.0, endwtime;
  int  namelen;
  char processor_name[MPI_MAX_PROCESSOR_NAME];
  MPI_Status status;
  double sumI;

  MPI_Init(&argc,&argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myid);
  MPI_Get_processor_name(processor_name,&namelen);

  printf("Process %d on %s: n=%d\n",myid, processor_name,numprocs);
  if( myid == 0 ) 
  {
    scanf("%ld", &n_intervals);
    printf("Using %ld intervals\n",n_intervals);
    startwtime = MPI_Wtime();
  }
  n = n_intervals;
  
  if(myid == 0)
    for(i = 1; i < numprocs; i++)
      MPI_Send(&n, 1, MPI_LONG, i, 10, MPI_COMM_WORLD);
  else  
    MPI_Recv(&n, 1, MPI_LONG, 0, 10, MPI_COMM_WORLD, &status);

  h   = 1.0 / (double) n;
  sum = 0.0;
  for (i = myid + 1; i <= n; i += numprocs)
  {
      x = h * ((double)i - 0.5);
      sum += f(x);
  }
  mypi = h * sum;

  if(myid != 0)
      MPI_Send(&mypi, 1, MPI_DOUBLE, 0, 20, MPI_COMM_WORLD);
  else  
  {
    for(i = 1; i < numprocs; i++)
    {
      MPI_Recv(&sumI, 1, MPI_DOUBLE, i, 20, MPI_COMM_WORLD, &status);
      mypi += sumI;
    }
  }
  
  if (myid == 0)
  {
      printf("pi is approximately %.16f, Error is %.16f\n", mypi, fabs(pi - PI25DT));
      endwtime = MPI_Wtime();
      printf("wall clock time = %f\n", endwtime-startwtime);
  }

  MPI_Finalize();

  return 0;
}

