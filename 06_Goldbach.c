#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

const int A = 1000000;

int *allocateI(int n)
{
	int *a;
	a = (int *) malloc((unsigned)n*sizeof(int));
	return a;
}

int prime(int n)
{
	int i;
	if((n == 0) || (n == 1)) return 0;
	for(i = 2; i <= n/2; i++) if((n % i) == 0) return 0;
	return 1;
}


int main(int argc, char *argv[])
{
	int i, j, theorem, nfind, t1;
	int id; //Process rank
	int np; //Num of processes
	int *primes, *globalPrimes;
	double elapsed_time;
	

	MPI_Init(&argc, &argv);
	MPI_Barrier(MPI_COMM_WORLD);
	elapsed_time = - MPI_Wtime();
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &np);
	
	
	primes = allocateI(A);
	globalPrimes = allocateI(A);
	
	if(id == 0)
		for(i = 0; i < A; i++) primes[i] = 0;
		
	MPI_Bcast(primes, A, MPI_INT, 0, MPI_COMM_WORLD);
	
	for(i = 2*id+1; i < A; i += 2*np)
		if(prime(i)) primes[i] = 1;
	
	//printf("Check point 01\n");
		
	MPI_Allreduce(primes, globalPrimes, A, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
	
	globalPrimes[2] = 1; // @ je prost broj, a nismo za parne ispitivali
	
	theorem = 1;
	
	for(i = 2*id+4; i <= A; i += 2*np)
	{
		nfind = 1;
		j = 1; 
		while((nfind) && (j < i/2))
		{
			if((globalPrimes[j]) && (globalPrimes[i-j])) nfind = 0;
			j++;			
		}
		if(nfind)
		{
			printf("Teorema ne vazi, pokazao proces br %d, za broj %d \n", id, i);
			theorem = 0;
		}		
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
	MPI_Reduce(&theorem, &t1, 1, MPI_INT, MPI_LAND, 0, MPI_COMM_WORLD);

	elapsed_time += MPI_Wtime();
	
	if(id ==0)
	{
		printf("Elapsed time is %f \n", elapsed_time);
		if(theorem) printf("Teorema je dokazana!!!!!!!\n");
	}
	MPI_Finalize();
	
	return 0;
}

