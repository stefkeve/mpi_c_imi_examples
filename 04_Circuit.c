/*
 * Circuit.c
 *
 *  Created on: Aug 20, 2010
 *      Author: Djordje
 */

#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

#define EXTRACT_BIT(n, i) ((n&(1<<i))?1:0)

int main(int argc, char *argv[])
{
  int i;
  int id; //Process rank
  int p; //Num of processes
  int solutions, solutionsI; /*Solutions found by this process*/
  int check_circuit (int, int);
  double elapsed_time;
  MPI_Status status;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
  
  if(id == 0)
    elapsed_time = - MPI_Wtime();

  solutions = 0;
  for(i = id; i < 65536; i += p)
    solutions += check_circuit(id, i);

  if(id == 0)
    for(i = 1; i < p; i++)
    {
      MPI_Recv(&solutionsI, 1, MPI_INT, i, 10, MPI_COMM_WORLD, &status);
      solutions += solutionsI;
    }
  else
    MPI_Send(&solutions, 1, MPI_INT, 0, 10, MPI_COMM_WORLD);
    
  if (id == 0)
    elapsed_time += MPI_Wtime();
  
  printf("Process %d is done\n", id);
  MPI_Finalize();
  

  if (id == 0) printf("There are %d different solutions. Elapsed time is %f\n", solutions, elapsed_time);
  fflush(stdout);
  return 0;
}

int check_circuit(int id, int z)
{
  int v[16];
  int i;

  for(i = 0; i <16; i++) v[i] = EXTRACT_BIT(z, i);

  if((v[0] || v[1]) && (!v[1] || !v[3]) && (v[2] || v[3])
      && (!v[3] || !v[4])  && (v[4] || !v[5])
      && (v[5] || !v[6])   && (v[5] || v[6])
      && (v[6] || !v[15])  && (v[7] || !v[8])
      && (!v[7] || !v[13]) && (v[8] || v[9])
      && (v[8] || !v[9])   && (!v[9] || !v[10])
      && (v[9] || v[11])   && (v[10] || v[11])
      && (v[12] || v[13])  && (v[13] || !v[14])
      && (v[14] || v[15]))
  {
  printf("%d) %d%d%d%d%d%d%d%d%d%d%d%d%d%d%d%d\n", id, v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10], v[11], v[12], v[13], v[14], v[15]);

    return 1;
  }
  
  return 0;
}
