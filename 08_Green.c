#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

//const int A = 1000000;

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

int digitsSum(int n)
	{
		int s = 0;
		while(n > 0)
		{
			s += (n % 10);
			n = n /10;
		}
		return s;
	}

void swap(int *px, int *py)
{
  int tmp;
  tmp = *px;
  *px = *py;
  *py = tmp;
}
                      
/* sortira niz u rastuci */
void sort(int *a, int n)
{
  int i, j;
                        
  for(i = 0; i < n-1; i++)
      for(j = i+1; j < n; j++)
          if(a[i] > a[j]) swap(&a[i], &a[j]);
}


int main(int argc, char *argv[])
{
	FILE *fOut;
	int A;
	int i, j, nGreens;
	int size,mySize;
	int id; //Process rank
	int np; //Num of processes
	int *primes, *globalPrimes, *greens;
	double elapsed_time;
	MPI_Status status;
	
	
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	if(id == 0)
	{
		fOut = fopen("Greens.dat", "wt");
		printf("Unesite broj do kog zelite da trazite zelene brojeve!\n");
		//scanf("%d",&A);
		A = 1000000;
		primes = allocateI(A);
		for(i = 0; i < A; i++) primes[i] = 0;
		mySize = 0;
	}	

	MPI_Barrier(MPI_COMM_WORLD);
	if(id == 0)
	  elapsed_time = - MPI_Wtime();


	MPI_Bcast(&A, 1, MPI_INT, 0, MPI_COMM_WORLD);
	if(id != 0) primes = allocateI(A);
	MPI_Bcast(primes, A, MPI_INT, 0, MPI_COMM_WORLD);
	globalPrimes = allocateI(A);
	
	for(i = 2*id+3; i < A; i += 2*np)
		if(prime(i)) primes[i] = 1;
		
	MPI_Allreduce(primes, globalPrimes, A, MPI_INT, MPI_LOR, MPI_COMM_WORLD);
		
	globalPrimes[2] = 1; // @ je prost broj, a nismo za parne ispitivali

	greens = allocateI(A);

	nGreens = 0;
	
	if(id == 0)
	{
	  greens[nGreens++] = 2;
	}

	for(i = 2*id+3; i < A; i += 2*np)
		if((globalPrimes[digitsSum(i)] == 1) && (globalPrimes[i])) greens[nGreens++] = i;
	
	if (id != 0)
		MPI_Send(greens, nGreens, MPI_INT, 0, 20, MPI_COMM_WORLD);
	else
	{
		mySize = nGreens;

		for (i = 1; i < np; i++)
			{
				MPI_Probe(i, 20, MPI_COMM_WORLD, &status);
				MPI_Get_count(&status, MPI_INT, &size);//size = status.size;
				MPI_Recv(&greens[mySize], size, MPI_INT, i, 20, MPI_COMM_WORLD, &status);
				mySize += size;
			}
		sort(greens, mySize);	
	}
	
	if(id ==0)
	{
		elapsed_time += MPI_Wtime();
		printf("Elapsed time is %f \n", elapsed_time);
		for (i = 0; i < mySize; i++) fprintf(fOut,"%d\n", greens[i]);
		fclose(fOut);
	}

	MPI_Finalize();
	
	return 0;
}

