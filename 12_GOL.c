#include "mpi.h"   
#include <stdio.h> 
#include <stdlib.h>

#define NI 1000
#define NJ 1000
#define NSTEPS 5000

int* allocate_1D(int n)
{
  int *a;
  a = (int *) malloc((unsigned)n*sizeof(int));
  return a;
}

void printMat(int *a, int n, int m)
{
  int i,j;
  
  for(i = 0; i < n; i++)
  {
    for(j = 0; j < m; j++) printf("%3d", a[i*m+j]);
    printf("\n");
  }  
}

main(int argc, char* argv[])
{
  int i, j, n, nprocs, id, poproc, ii, steps;
  int alive, globSum; // broj zivih celija na svakom procesu i na ukupan broj zivih celija
  int ni, nj, nn; // dimenzije matrice sa host celijama, nn=ni*nj
  int nni, nnj, nnp; // dimenzije mpodmatrice na svakom procesu, nnp=nni*nnj
  int *old, *help, *myold, *mynew;
  MPI_Status status;
  float x;
  double elapsedTime;
              
  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD, &nprocs);
  MPI_Comm_rank(MPI_COMM_WORLD, &id);

  poproc = NI/nprocs;//broj vrsta po kolona po procesu....i ide po x osi a j po y osi

  if(id == 0)
  {
    n = NI*NJ;
    help = allocate_1D(n);//pomocna "matrica"
    
    srand ( time(NULL) );
    
    for(i = 0; i < n; i++)
    {
      x = (double)rand()/(RAND_MAX);
      if(x < 0.5) help[i] = 0;
      else help[i] = 1;
    }
    
    /* dimenzije nove matrice sa ghost prvom i poslednjom kolonom */
    ni = NI + 2;
    nj = NJ + 2;
    nn = ni*nj;
    old = allocate_1D(nn);
    
    for(i = 0; i < nn; i++) old[i] = 0;
    
    ii = 0;
    for(i = 1; i < ni-1; i++)
      for(j = 1; j < nj-1; j++) 
        old[i*nj + j] = help[ii++];    
    
    elapsedTime = -MPI_Wtime();
    
    for(j = 1; j < nj-1; j++) old[j] = old[(ni-2)*nj + j]; /* posledjnu ne ghost kolonu kopiramo u nultu (ghost) kolonu */
    for(j = 1; j < nj-1; j++) old[(ni-1)*nj + j] = old[nj + j]; /* prvu ne ghost kolonu kopiramo u poslednju(ghost) kolonu */
    
    for(i = 1; i < ni-1; i++) old[i*nj] = old[(i+1)*nj -2]; /* posledjnu ne ghost vrstu kopiramo u nultu (ghost) vrstu */
    for(i = 1; i < ni-1; i++) old[(i+1)*nj-1] = old[i*nj+1]; /* prvu ne ghost vrstu kopiramo u poslednju(ghost) vrstu */
    
    old[0] = old[(ni-1)*nj-2]; old[ni*nj-1] = old[nj+1]; // coskovi...glavna dijagonala
    old[nj-1] = old[(ni-2)*nj+1]; old[(ni-1)*nj] = old[2*nj-2]; // coskovi...sporedna dijagonala
    
    /* dimenzije matrice na svakom procesu */
    nni = poproc+2;
    nnj = nj;
    nnp = nni * nnj;
    myold = allocate_1D(nnp);
    mynew = allocate_1D(nnp);
    
    //spremanje kolona za sebe
    for(i = 0; i < nnp; i++) myold[i] = old[i];
    
    /*poslati informacije o broju elemenata koji ce da imaju. i deo niza
      sa kojim ce oni rukovati*/  
    for(i = 1; i < nprocs; i++)
    {
      MPI_Send(&nnp, 1, MPI_INT, i, 20, MPI_COMM_WORLD);
      MPI_Send(&nni, 1, MPI_INT, i, 20, MPI_COMM_WORLD);
      MPI_Send(&nnj, 1, MPI_INT, i, 20, MPI_COMM_WORLD);
      MPI_Send(&old[i*poproc*nnj], nnp, MPI_INT, i, 30, MPI_COMM_WORLD);
    }
  }
  else
  {
    MPI_Recv(&nnp, 1, MPI_INT, 0, 20, MPI_COMM_WORLD, &status);
    MPI_Recv(&nni, 1, MPI_INT, 0, 20, MPI_COMM_WORLD, &status);
    MPI_Recv(&nnj, 1, MPI_INT, 0, 20, MPI_COMM_WORLD, &status);
    
    myold = allocate_1D(nnp);
    mynew = allocate_1D(nnp);        
    MPI_Recv(myold, nnp, MPI_INT, 0, 30, MPI_COMM_WORLD, &status);
  }
  
  for(steps = 0; steps < NSTEPS; steps++)
  {
    /* u svakom koraku provera zivota za svaku celiju */
    for(i = nnj; i < (nnp-1); i++)
      if(!((i%nnj == 0)||(i%(nnj) == nnj-1))) // da ne pitamo za GHOST celije
      {
        alive = 0;
        alive = myold[i-nnj+1] + myold[i+1] + myold[i+nnj+1] + 
                myold[i-nnj]   +              myold[i+nnj]   +
                myold[i-nnj-1] + myold[i-1] + myold[i+nnj-1] ;
        switch(alive)
          {
            case 3:
              mynew[i] = 1; /* ako ima zive tri celije u okolini postaje ziv */
              break;
                                      
            case 2:
              mynew[i] = myold[i]; /* ako ima dve zive celije u okolini, ostaje kao i u prethodnom koraku */
              break;
                                                                   
            default:
              mynew[i] = 0; /* u suprotnom umire */
          }        
                                                                                                  
      }
    
    for(i = nnj; i < (nnp-1); i++) myold[i] = mynew[i]; //prepisivanje u staru matricu
    
    for(i = 1; i < nni-1; i++)
    {
      myold[i*nnj] = myold[(i+1)*nnj-2]; //     granicni
      myold[(i+1)*nnj-1] = myold[i*nnj+1]; //    uslovi 
    }
    
    if(id == 0)
    {
      //myold[0] = myold[(nni - 1)*nnj-2]; myold[nni*nnj-1] = myold[nnj+1];
      //myold[nnj-1] = myold[(nni-2)*nnj+1]; myold[(nni-1)*nnj] = myold[2*nnj-2];
      if(nprocs !=1 )
      {
        MPI_Send(&myold[nnj], nnj, MPI_INT, nprocs-1, 50, MPI_COMM_WORLD);
        MPI_Recv(&myold[0], nnj, MPI_INT, nprocs-1, 50, MPI_COMM_WORLD, &status);
      
        MPI_Send(&myold[(nni-2)*nnj], nnj, MPI_INT, id+1, 50, MPI_COMM_WORLD);
        MPI_Recv(&myold[(nni-1)*nnj], nnj, MPI_INT, id+1, 50, MPI_COMM_WORLD, &status);
      }
    }  
    else if(id == nprocs-1)
         {
           MPI_Recv(&myold[(nni-1)*nnj], nnj, MPI_INT, 0, 50, MPI_COMM_WORLD, &status);
           MPI_Send(&myold[(nni-2)*nnj], nnj, MPI_INT, 0, 50, MPI_COMM_WORLD);
         
           MPI_Recv(&myold[0], nnj, MPI_INT, id-1, 50, MPI_COMM_WORLD, &status);
           MPI_Send(&myold[nnj], nnj, MPI_INT, id-1, 50, MPI_COMM_WORLD);
         }       
	 else
	 {
		 MPI_Recv(&myold[0], nnj, MPI_INT, id-1, 50, MPI_COMM_WORLD, &status);
		 MPI_Send(&myold[nnj], nnj, MPI_INT, id-1, 50, MPI_COMM_WORLD);
		 
		 MPI_Send(&myold[(nni-2)*nnj], nnj, MPI_INT, id+1, 50, MPI_COMM_WORLD);
		 MPI_Recv(&myold[(nni-1)*nnj], nnj, MPI_INT, id+1, 50, MPI_COMM_WORLD, &status);
	 }
  }  
  
  alive = 0;
  
  /* prebrajamo zive celije */
  for(i = nnj; i < (nnp-1); i++)
    if(!((i%nnj == 0)||(i%(nnj) == nnj-1))) // da ne pitamo za GHOST celije
      alive += myold[i];
  
  MPI_Reduce(&alive, &globSum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
  elapsedTime += MPI_Wtime();                 
  if(id == 0) printf("Number of live cells = %d,  and elapsedTime = %f\n", globSum, elapsedTime);
  
  MPI_Finalize();
}
                     