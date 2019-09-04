#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

const int A = 0;
const int B = 1;
double f(double x)
{
   int i;
   double s = -2;
   double p = 1;
 
   for(i = 0; i < 1000; i++)
   {
      p *= x;
      s += sin(p);
   }
   return s;
}

int main( int argc, char *argv[])
{
   double limit = 1.E-11;
   double left, right, part;
   int  id, np, flag, i;
   int *flags;
   double globalL, globalR;
   double elapsed_time;

   MPI_Init(&argc,&argv);
   MPI_Comm_size(MPI_COMM_WORLD,&np);
   MPI_Comm_rank(MPI_COMM_WORLD,&id);

   if (id == 0)
   {
      flags = (int*)malloc((unsigned)np * sizeof(int));
      elapsed_time = -MPI_Wtime();
   }
   
   globalL = 0.0;
   globalR = 1.0;
   part = globalR - globalL;
   flag = 0;
   //printf("%13.12f\n",limit);
   
   while(part > limit)
   {
      part = part/np;
      left = globalL + id*part;
      right = left + part;
      //printf("id = %d and I work from %13.12f to %13.12f\n", id, left, right);
      if((f(left) * f(right)) < 0) flag = 1;
      //printf("id is %d, fl=%13.12f, fr=%13.12f\n", id, f(left), f(right));
      MPI_Gather(&flag, 1, MPI_INT, flags, 1, MPI_INT, 0, MPI_COMM_WORLD);
      
      if(id == 0)
      {
         i = 0;
         while(flags[i] == 0) i++; //printf("i=%d\n",i);
         globalL += i * part;
         globalR = globalL + part;
      }
      
      MPI_Bcast(&globalL, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      MPI_Bcast(&globalR, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
      
      flag = 0;
      //if(id == 0)printf("%13.12f\t%13.12f\n", globalL, globalR);
   }
   
   if (id == 0)
   {
      elapsed_time += MPI_Wtime();
      printf("smallest positive root of equation f is %10.9f, and the elapsed time is %13.12f, f(x) = %13.12f\n", globalL, elapsed_time, f(globalL));
   }
   MPI_Finalize();
   return 0;
}
