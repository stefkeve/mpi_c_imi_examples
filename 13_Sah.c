#include <mpi.h>
#include <stdio.h>
#include <time.h>
#include <stdlib.h>

const double dt = 0.0001;

typedef struct Figura_s
{
	int id;
	int tip; // 0 - top, 1 - lovac, 2 - piun, -1 - pojedena
	int x;
	int y;
	int brojPojedenih;
	int brojPoteza;
} Figura;

typedef struct Polje_s
{
	int id;
	int tip;
} Polje;

Figura* allocate1DFigura(int n)
{
	Figura *figure;
	figure = (Figura*)malloc(n * sizeof(Figura));

	return figure;
}

Figura** allocate2DFigura(int m, int n)
{
	Figura **figure;
	figure = (Figura**)malloc(m * sizeof(Figura*));
	for (int i = 0; i < m; i++) figure[i] = (Figura*)malloc(n * sizeof(Figura));

	return figure;
}

Polje** allocate2DPolje(int m, int n)
{
	Polje **polje;
	polje = (Polje**)malloc(m * sizeof(Polje*));

	for(int i = 0; i < m; i++)
		polje[i] = (Polje*)malloc(n * sizeof(Polje));

	return polje;
}

void movePiun(int *x, int *y, int m)
{
	int r = rand() % 4;

	if (r == 0) 
	{
		if(x > 0) x--;
	}
	else if (r == 1)
	{
		if(y < m-1) y++;
	}
	else if (r == 2)
	{
		if(x < m-1) x++;
	}
	else if (r == 3)
	{
		if (y > 0) y--;
	}
}

void moveLovac(int *x, int *y, int m)
{
	int r = rand() % 4;

	if (r == 0) 
	{
		if(x > 0 && y > 0) 
		{
			x--;
			y--;
		}
	}
	else if (r == 1)
	{
		if(x > 0 && y < m-1) 
		{
			x--;
			y++;
		}
	}
	else if (r == 2)
	{
		if(x < m-1 && y < m-1) 
		{
			x++;
			y++;
		}
	}
	else if (r == 3)
	{
		if (x < m-1 && y > 0) 
		{
			x++;
			y--;
		}
	}
}

void moveTop(int *x, int *y, int m)
{
	int r = rand() % 2;

	if (r == 0) 
	{
		int p = rand() % m;
		x = p;
	}
	else if (r == 1)
	{
		int p = rand() % m;
		y = p;
	}
}

int* allocate1DInt(int n)
{
	int *a;
	a = (int *)malloc(n * sizeof(int));

	return a;
}

void addFigure(Figura *a, int *n, Figura *b, int m)
{
	int i;

	for(i = 0; i < m; i++)
		a[*n + i] = b[i];

	(*n) += m;
}

int deleteFigura(int id, Figura *figure, int *n)
{
	int i;

	for (i = 0; i < *n; i++)
	if (figure[i].id == id) 
	{
		figure[i] = figure[(*n)-1];
		(*n)--;
		return 1;	
	}
	return 0;
}

void UvecajBrojPojedenih(Figura *f, int n, int id)
{
	int i;

	for (i = 0; i < n; i++)
		if (f[i].id == id)
		{
			f[i].brojPojedenih++;
			printf("id = %d, pr poj = %d\n", f[i].id, f[i].brojPojedenih);
			break;
		}
}

int main(int argc, char *argv[])
{
	int id;
	int np;

	int i, j, x, y, r;

	int m, mm, nn;  //mm-broj vrsta po procesu, nn-broj kolona po procesu, nn=n;
	int brojFigura, ukupnoFigura;
	Figura **saljem, *primam;
	Polje **polja; //, **poljaNova;
	Figura *figure;
	int *kolikoSaljem, kolikoPrimam;
	MPI_Status status;

	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &id);
	MPI_Comm_size(MPI_COMM_WORLD, &np);

	// printf("caooooo %d\n",id);

	MPI_Datatype figuraType;

	MPI_Type_contiguous(5, MPI_INT, &figuraType);
	MPI_Type_commit(&figuraType);

	if (id == 0)
	{
		printf("Unesite koliko figura po procesu zelite: \n");
		scanf("%d", &brojFigura);
		printf("Unesite dimenzije kvadratne matrice m x m: \n");
		scanf("%d", &m);
	}

	MPI_Bcast(&brojFigura, 1, MPI_INT, 0, MPI_COMM_WORLD);
	MPI_Bcast(&m, 1, MPI_INT, 0, MPI_COMM_WORLD);

	mm = m;
	nn = m/np;

	polja = allocate2DPolje(mm, nn);
	figure = allocate1DFigura(np*brojFigura);
	saljem = allocate2DFigura(np, np*brojFigura);
	primam = allocate1DFigura(np*brojFigura);
	kolikoSaljem = allocate1DInt(np);

	for(i = 0; i < mm; i++)
		for (j = 0; j < nn; j++)
			polja[i][j].id = -1;
	
	srand(time(NULL)*id);

	for(i = 0; i < brojFigura; i++)
	{
		do 
		{
			x = rand() % mm;
			y = rand() % nn;
		} while (polja[x][y].id != -1);
		
		r = rand() % 3;

		figure[i].id = id*brojFigura + i;
		figure[i].x = x;
		figure[i].y = id*nn + y;
		figure[i].brojPojedenih = 0;
		figure[i].brojPoteza    = 0;
		figure[i].tip = r;

		polja[x][y].id = figure[i].id;
		polja[x][y].tip = figure[i].tip;
	}

	MPI_Allreduce(&brojFigura, &ukupnoFigura, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);

	printf("caooooo %d, ukupan broj figura je %d\n",id, ukupnoFigura);

	while (ukupnoFigura > 1)
	{
//		printf("Javljam se na raport, i imam %d figura\n", brojFigura);

/*		printf("Pre:\n");
		for (i = 0; i < mm; i++)
		{
			for (j = 0; j < nn; j++)
				if(polja[i][j].id != -1)printf("%d\t", polja[i][j].tip);
				else printf("%d\t", -1);
			printf("\n");
		}
*//*
		for (i = 0; i < brojFigura; i++) printf("%d\t",figure[i].tip);
		printf("\n");
*/
		for (i = 0; i < np; i++)
			kolikoSaljem[i] = 0;
		

		for(i = 0; i < brojFigura; i++)
		{
			polja[figure[i].x][figure[i].y%nn].id = -1;

			figure[i].brojPoteza++;

			switch(figure[i].tip)
			{
			case 0:
				moveTop(figure[i].x, figure[i].y, m);
				break;
			case 1:
				moveLovac(figure[i].x, figure[i].y, m);
				break;
			case 2:
				movePiun(figure[i].x, figure[i].y, m);
				break;
			}

			if ((figure[i].y < (id*nn)) || (figure[i].y >= ((id+1)*nn)))
			{
				saljem[y/nn][kolikoSaljem[y/nn]++] = figure[i];
				deleteFigura(figure[i].id, figure, &brojFigura);
				i--;
			}
			
		}

		if(np > 1)
		{
			for(i = 0; i < np; i++)
				if (i < id)//recv pa send
				{
					MPI_Probe(i, 10, MPI_COMM_WORLD, &status);
					MPI_Get_count(&status, figuraType, &kolikoPrimam);
					MPI_Recv(primam, kolikoPrimam, figuraType, i, 10, MPI_COMM_WORLD, &status);
					
					addFigure(figure, &brojFigura, primam, kolikoPrimam);

					MPI_Send(saljem[i], kolikoSaljem[i], figuraType, i, 10, MPI_COMM_WORLD);
				}
				else if (i > id)//send pa recv
				{
					MPI_Send(saljem[i], kolikoSaljem[i], figuraType, i, 10, MPI_COMM_WORLD);

					MPI_Probe(i, 10, MPI_COMM_WORLD, &status);
					MPI_Get_count(&status, figuraType, &kolikoPrimam);
					MPI_Recv(primam, kolikoPrimam, figuraType, i, 10, MPI_COMM_WORLD, &status);

					addFigure(figure, &brojFigura, primam, kolikoPrimam);
				}
		}

		for(i = 0; i < brojFigura; i++)
		{
			int xx = figure[i].x;
			int yy = figure[i].y%nn;
			if(polja[xx][yy].id == -1) 
			{
				polja[xx][yy].id = figure[i].id;
				polja[xx][yy].tip = figure[i].tip;
			}
			else 
			{
				switch (polja[xx][yy].tip)
				{
					// na polju se vec nalazi top
				case 0:
					switch (figure[i].tip)
					{					
					case 0:
						UvecajBrojPojedenih(figure, brojFigura, polja[xx][yy].id);
						deleteFigura(figure[i].id, figure, &brojFigura);
						i--;
						printf("Top jede topa\n");
						break;
					case 1:
						UvecajBrojPojedenih(figure, brojFigura, polja[xx][yy].id);
						deleteFigura(figure[i].id, figure, &brojFigura);
						i--;
						printf("Top jede lovca\n");
						break;
					case 2:
						deleteFigura(polja[xx][yy].id, figure, &brojFigura);
						UvecajBrojPojedenih(figure, brojFigura, figure[i].id);
						polja[xx][yy].id = figure[i].id;
						polja[xx][yy].tip = figure[i].tip;
						printf("Piun jede topa\n");
						i--;
						break;
					}				
					break;

					// na polju se vec nalazi lovac
				case 1:
					switch (figure[i].tip)
					{
					case 2:
						UvecajBrojPojedenih(figure, brojFigura, polja[xx][yy].id);
						deleteFigura(figure[i].id, figure, &brojFigura);
						i--;
						printf("Lovac jede piuna\n");
						break;
					case 1:
						UvecajBrojPojedenih(figure, brojFigura, polja[xx][yy].id);
						deleteFigura(figure[i].id, figure, &brojFigura);
						i--;
						printf("Lovac jede piuna\n");
						break;
					case 0:
						deleteFigura(polja[xx][yy].id, figure, &brojFigura);
						UvecajBrojPojedenih(figure, brojFigura, figure[i].id);
						polja[xx][yy].id = figure[i].id;
						polja[xx][yy].tip = figure[i].tip;
						printf("Top jede lovca\n");
						i--;
						break;
					}				
					break;

					// na polju se vec nalazi piun
				case 2:
					switch (figure[i].tip)
					{
					case 0:
						UvecajBrojPojedenih(figure, brojFigura, polja[xx][yy].id);
						deleteFigura(figure[i].id, figure, &brojFigura);
						i--;
						printf("Piun jede topa\n");
						break;
					case 2:
						UvecajBrojPojedenih(figure, brojFigura, polja[xx][yy].id);
						deleteFigura(figure[i].id, figure, &brojFigura);
						i--;
						printf("Piun jede piuna\n");
						break;
					case 1:
						deleteFigura(polja[xx][yy].id, figure, &brojFigura);
						UvecajBrojPojedenih(figure, brojFigura, figure[i].id);
						polja[xx][yy].id = figure[i].id;
						polja[xx][yy].tip = figure[i].tip;
						printf("Lovac jede piuna\n");
						i--;
						break;
					}				
					break;
				}			
			}
		}

		printf("Javljam se na raport, i imam %d figura\n", brojFigura);

/*		printf("Posle:\n");
		for (i = 0; i < mm; i++)
		{
			for (j = 0; j < nn; j++)
				if(polja[i][j].id != -1)printf("%d\t", polja[i][j].tip);
				else printf("%d\t", -1);
				printf("\n");
		}
*/

		MPI_Allreduce(&brojFigura, &ukupnoFigura, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
	}

	if (brojFigura > 0) 
	{
		switch (figure[0].tip)
		{
		case 0:
			printf("Ja sam proces %d. Pobednik je top! Pojeo sam %d figura!\n", id, figure[0].brojPojedenih);
			break;
		case 1:
			printf("Ja sam proces %d. Pobednik je lovac! Pojeo sam %d figura!\n", id, figure[0].brojPojedenih);
			break;
		case 2:
			printf("Ja sam proces %d. Pobednik je piun! Pojeo sam %d figura!\n", id, figure[0].brojPojedenih);
			break;
		}
	}

	printf("Zavrsio proces %d!\n", id);

	MPI_Finalize();

	return 0;
}
