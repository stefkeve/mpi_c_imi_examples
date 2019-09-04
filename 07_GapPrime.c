#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

const int A = 1000000;

int *allocate_1D(int n)
{
   int *a;
   a = (int *) malloc((unsigned)n*sizeof(int));
   return a;
}
      
/* f-ja koja vraca 1 ukoliko je n prost broj a 0 ako nije */
int prime(int n)
{
   int i;
   if((n == 0) || (n == 1)) return 0;
   if(n == 2) return 1;
   for(i = 2; i <= n/2; i++) if((n % i) == 0) return 0;
   return 1;
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
                                    
/* najveca razlika izmedju dva susedna broja u nizu */
int spaceMax(int *a, int n)
{
  int i, max;
                        
  max = a[1] - a[0];
  for(i = 1; i < n-1; i++)
      if(a[i+1] - a[i] > max) max = a[i+1] - a[i];
  return max;
}
                                           
int main(int argc, char *argv[])
{
  int i, j;
  int max; /* konacno resenje, najvece rastojanje izmedju dva prosta broja  */
  int lastCount, pomCount;
  int k; /* broj prostih brojeva po procesoru */
  int globalSolutionCount; /* ukupan broj prostih brojeva */
  int *partialPrimes; /* niz prostih brojeva na svakom procesoru */
  int *globalPrimes; /* konacan niz prostih brojeva */
  int id; /* Process rank */
  int p; /* Num of processes */
  double elapsed_time; /* proteklo vreme */
  MPI_Status status;
                                                
  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);
  MPI_Comm_size(MPI_COMM_WORLD, &p);
                                                                                
  if(id == 0) 
  {
    elapsed_time = - MPI_Wtime();
    globalPrimes = allocate_1D(A/2);
  }
  else partialPrimes = allocate_1D(A/(2*p));
                                                                                    
  k = 0;
                                                                                      
  for(i = 2 * id + 1; i < A; i += 2*p)
      if(prime(i))
      {
         if(id == 0) globalPrimes[k++] = i;
            else partialPrimes[k++] = i; 
      }
                                                                                                 
  MPI_Reduce(&k, &globalSolutionCount, 1 , MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
                                                                                                                 
  if(id != 0)
  {
     MPI_Send(&k, 1, MPI_INT, 0, 20, MPI_COMM_WORLD);
     MPI_Send(partialPrimes, k, MPI_INT, 0, 21, MPI_COMM_WORLD);
  }
     else
     {
         lastCount = k;
         for(j = 1; j < p; j++)
         {
             MPI_Recv(&pomCount, 1, MPI_INT, j, 20, MPI_COMM_WORLD, &status);
             MPI_Recv(&globalPrimes[lastCount], pomCount, MPI_INT, j, 21, MPI_COMM_WORLD, &status);
             lastCount += pomCount;
         }
     }
                                                                                                                                                 
  if(id == 0)
  {
    double t = - MPI_Wtime();
    sort(globalPrimes, globalSolutionCount);
    max = spaceMax(globalPrimes, globalSolutionCount);
    t += MPI_Wtime();
    elapsed_time += MPI_Wtime();
    printf("Sort time + find time is %lf\n", t);
    printf("%d, Elapsed time is %lf\n", max, elapsed_time);
  }																			   									

  if(id == 0) free(globalPrimes);
    else free(partialPrimes);
    
  MPI_Finalize();
    
  return 0;
}