/*
 * ArraySum.c
 *
 *  Created on: Aug 20, 2010
 *      Author: djordje
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

const long int A = 100000000;

int *allocate_1D(int n)
{
	int *a;
	a = (int *) malloc((unsigned)n*sizeof(int));
	return a;
}

int* postaviPoKoliko(int n, int brproc)
{
	int i, poproc, ostatak;
	int *po;
	
	po = allocate_1D(brproc);
	
	poproc = n / brproc;
	ostatak = n % brproc;

	for(i = 0; i < brproc; i++)
	  if(i >= brproc-ostatak)
	    po[i] = poproc + 1;
	  else po[i] = poproc;
  
	return po;
}

int* postaviOdakleKo(int n, int brproc, int *po)
{
	int i;
	int *od;
	
	od = allocate_1D(brproc);
	
	od[0] = 0;
	for(i = 1; i < brproc; i++)
	{
	    od[i] = od[i-1] + po[i-1];
	}
	
	return od;
}

int sumaNiza(int *a, int n)
{
	int i, sum = 0;
	for(i = 0; i < n; i++) sum += a[i];
	return sum;
}

int* array_input(int n)
{
	int i, *a;
	a = allocate_1D(n);
	for(i = 0; i < n; i++)
	  a[i] = ((double)rand() / RAND_MAX) * 10;
	return a;
}

int main(int argc, char *argv[])
{
	int global_sum;
	int i, n, *a, *pocTacke, *svakoPo, nn;
	int id; //Process rank
	int np; //Num of processes
	int sum;
	double elapsed_time;
	MPI_Status status;

	srand(0);

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	if(id == 0)
	{
	  a = array_input(A);

	  svakoPo = postaviPoKoliko(A,np);
	  pocTacke = postaviOdakleKo(A,np,svakoPo);
	  	    
	  for(i = 0; i < np; i++)
	  printf("proc %d, radi od %d, i radi %d elemenata!\n", i, pocTacke[i], svakoPo[i]);
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	elapsed_time = - MPI_Wtime();
	
	MPI_Scatter(svakoPo, 1, MPI_INT, &nn, 1, MPI_INT, 0, MPI_COMM_WORLD);

	if(id == 0)
	  for(i = 1; i < np; i++)
	    MPI_Send(&a[pocTacke[i]], svakoPo[i], MPI_INT, i, 17, MPI_COMM_WORLD);
	else
	{
	  a = allocate_1D(nn);
	  MPI_Recv(a, nn, MPI_INT, 0, 17, MPI_COMM_WORLD, &status);
	}

	sum = sumaNiza(a, nn);

	MPI_Reduce(&sum, &global_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	elapsed_time += MPI_Wtime();
	
	free(a);
	if (id == 0) 
	{
	  printf("Sum of array is %d. Elapsed time is %f\n", global_sum, elapsed_time);
	  free(pocTacke);
	}
	
	MPI_Finalize();

	return 0;
}




