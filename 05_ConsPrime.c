/*
 * Prime1.c
 * prost ciji je prvi neparni naslednik takodje prost
 *
 *  Created on: Aug 20, 2010
 *      Author: djordje
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

const int A = 1000000;

int prime(int n)
{
	int i;
	if((n == 0) || (n == 1)) return 0;            
        if(n == 2) return 1;
	for(i = 2; i <= n/2; i++) if((n % i) == 0) return 0;

	return 1;
}

int main(int argc, char *argv[])
{
	int i, solutions, globalSolution;
	int id; //Process rank
	int p; //Num of processes
	double elapsed_time;

	MPI_Init(&argc, &argv);
	MPI_Barrier(MPI_COMM_WORLD);
	elapsed_time = - MPI_Wtime();
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	solutions = 0;

	for(i = 2 * id + 1; i < A; i += 2*p)
		if(prime(i) && prime(i + 2))
		{
			solutions++;
			//printf("%d\n", i);
		}
	
	MPI_Reduce(&solutions, &globalSolution, 1 , MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	elapsed_time += MPI_Wtime();
	fflush(stdout);
	MPI_Finalize();
	if(id == 0)printf("%d, Elapsed time is %f\n", globalSolution, elapsed_time);
	return 0;
}

