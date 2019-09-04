/*
 * Eratostenovo sito
 */
#include <mpi.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#define MIN(a,b) ((a) < (b) ? (a) : (b))

#define BLOCK_LOW(id, p, n) ((id) * (n) / (p))
#define BLOCK_HIGH(id, p, n) (BLOCK_LOW((id)+1, p, n) - 1)
#define BLOCK_SIZE(id, p, n) (BLOCK_HIGH(id, p, n) - BLOCK_LOW(id, p, n))
#define BLOCK_OWNER(index, p, n) (((p) * ((index) + 1) - 1) / (n))

const int A = 100000000;

int main(int argc, char *argv[])
{
	int count; /* Local prime count */
	double elapsed_time; /* Parallel execution time */
	int first; /* Index of first multiple */
	int global_count; /* Global prime count */
	int high_value; /* Highest value on this proc */
	int i;
	int id; /* Process ID number */
	int index; /* Index of current prime */
	int low_value; /* Lowest value on this proc */
	char *marked; /* Portion from 2, 3, ..., 'n' */
	int n; /* Sieving of 2, 3, ..., 'n' */
	int p; /* Number of process */
	int proc0_size; /* Size of proc 0's subarray  */
	int prime; /* Current prime */
	int size; /* size */

	MPI_Init(&argc, &argv);

	/* Start the timer */

	MPI_Barrier(MPI_COMM_WORLD);
	elapsed_time = - MPI_Wtime();

	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &p);

	n = A;

	/* svaki procesor saznaje pocetni i krajnji element sa kojim
	 * raspolaze, i velcinu tog niza */

	low_value = 2 + BLOCK_LOW(id, p, n-1);
	high_value = 2 + BLOCK_HIGH(id, p, n-1);
	size = BLOCK_SIZE(id, p, n-1);
	
	proc0_size = (n-1) / p;

	if((1 + proc0_size) < (int) sqrt((double)n))
	{
		if(!id) printf("To many processes\n");
		MPI_Finalize();
		exit(1);
	}
	
//	printf("id = %d \tlow = %d \t high = %d\n", id, low_value, high_value);
	
	/* Allocate this process's share of the array */
	marked = (char *) malloc (size);

	if(marked == NULL)
	{
		printf("Cannot allocate enough mamory\n");
		MPI_Finalize();
		exit(1);
	}

	for(i = 0; i < size; i++) marked[i] = 0;
	
	if(!id) index = 0;
	prime =  2;
	
	/*
	* Neka je d = sqrt(n), odnosno d*d=n. Kada bi postojao neki slozen broj r > d, tada bi on morao biti deljiv sa 
	* prirodnim brojem (b) manjim od d, sto znaci da bi vec bio markiran pri markiranju brojeva deljivih sa b.
	* Tako da smo markiranjem umnozaka svih prostih brojeva do d, izgubili mogucnost da se izmedju d i n nadje ne 
	* prost broj koji nije markiran.
	*/

	do
	{
	  if(prime * prime > low_value)
	    first = prime * prime - low_value;
	  else
	  {
	    if(!(low_value % prime)) 
	      first = 0;
	    else 
	      first = prime - (low_value % prime);
	  }
	  
	  for(i = first; i < size; i+= prime) marked[i] = 1;
	  if(!id)
	  {
		  while(marked[++index]);
		  prime = index + 2;
	  }
	  MPI_Bcast(&prime, 1, MPI_INT, 0, MPI_COMM_WORLD);
	} while (prime * prime <= n);
	count = 0;
	for(i = 0; i < size; i++)
		if(!marked[i]) count++;
	MPI_Reduce(&count, &global_count, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

	/* Stop timer */

	elapsed_time += MPI_Wtime();

	/* Print the results */

	if(!id)
	{
		printf("%d primes are less then or equal to %d\n", global_count, n);
		printf("Total elapsed time : %10.6f\n", elapsed_time);
	}
	
	MPI_Finalize();
	return 0;
}


