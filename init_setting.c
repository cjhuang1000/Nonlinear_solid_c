#include "init_setting.h"
#include <math.h>
#include "stdlib.h"

int sumarray (int array[], int n);

void set_index(Index_S* ind, Grid_S* g, int **bs, char *fsineumanndirichlet){

	int N = g->N,	m = g->Nx,	n = g->Ny;
	int nint	=0,	nbnd	=0,		nfsi	=0, nneu	=0, ndir	=0;
	int n_xix	=0,	n_xiy	=0;
	int n_xfsi	=0,	n_yfsi	=0;
	int n_xneu	=0,	n_yneu	=0;
	int n_xdir	=0,	n_ydir	=0;
	int i,j,k,kk;

	unsigned int *cell_int,		*cell_bnd, *cell_fsi, *cell_neu, *cell_dir;
	int *xix_involved,	*xiy_involved;
	int *cell_involved_fsi_x,	*cell_involved_fsi_y;
	int *cell_involved_neu_x,	*cell_involved_neu_y;
	int *cell_involved_dir_x,	*cell_involved_dir_y;

	struct invol involved;

	cell_int = (unsigned int *) calloc(sizeof(unsigned int),N);
	cell_bnd = (unsigned int *) calloc(sizeof(unsigned int),N);
	cell_fsi = (unsigned int *) calloc(sizeof(unsigned int),N);
	cell_neu = (unsigned int *) calloc(sizeof(unsigned int),N);
	cell_dir = (unsigned int *) calloc(sizeof(unsigned int),N);

	xix_involved = (int *) calloc(sizeof(int),N);
	xiy_involved = (int *) calloc(sizeof(int),N);

	cell_involved_fsi_x = (int *) calloc(sizeof(int),N);
	cell_involved_fsi_y = (int *) calloc(sizeof(int),N);

	cell_involved_neu_x = (int *) calloc(sizeof(int),N);
	cell_involved_neu_y = (int *) calloc(sizeof(int),N);

	cell_involved_dir_x = (int *) calloc(sizeof(int),N);
	cell_involved_dir_y = (int *) calloc(sizeof(int),N);

	// part I:  Determine whether a cell is in, out or on the boundary
	for(j=1; j<n-1; j++)
		for(i=1; i<m-1;i++ )
		{

			// Not-external cells (at least 1 corner in)
			if ((bs[i][j] == -1) || (bs[i+1][j] == -1) ||
				(bs[i][j+1] == -1) || (bs[i+1][j+1] == -1))
			{
				kk = i+j*m;
				involved    = involvedIndices_grid( kk , m);
				for(k=0; k<6; k++)
				{
					xix_involved[involved.xix[k]] =1;
					xiy_involved[involved.xiy[k]] =1;
				}

				// interior
				if ((bs[i][j]+bs[i+1][j]+bs[i][j+1]+bs[i+1][j+1]) <= -3)
				{
					cell_int[nint] = kk;
					nint +=1;
				}
				// boundary cells
				else{
					cell_bnd[nbnd] = kk;
					nbnd +=1;

					switch (fsineumanndirichlet[kk])
					{
					case 'F':
						cell_fsi[nfsi] = kk;
						nfsi +=1;
						for(k=0; k<2; k++)
						{
							cell_involved_fsi_x[involved.stagx_cell[k]] =1;
							cell_involved_fsi_y[involved.stagy_cell[k]] =1;
						}
						break;
					case 'N':
						cell_neu[nneu] = kk;
						nneu +=1;
						for(k=0; k<2; k++)
						{
							cell_involved_neu_x[involved.stagx_cell[k]] =1;
							cell_involved_neu_y[involved.stagy_cell[k]] =1;
						}
						break;
					case 'D':
						cell_dir[ndir] = kk;
						ndir +=1;
						for(k=0; k<2; k++)
						{
							cell_involved_dir_x[involved.stagx_cell[k]] =1;
							cell_involved_dir_y[involved.stagy_cell[k]] =1;
						}
						break;

					}
				}

			}
		}

	// part II: Count the number of involved cells

	ind->cell_N_interior	= nint;
	ind->cell_N_boundary	= nbnd;
	ind->cell_N_Neumann	= nneu;
	ind->cell_N_Dirichlet	= ndir;
	ind->cell_N_FSI		= nfsi;

    ind->cell_interior = (unsigned int *) malloc(sizeof(unsigned int)*nint);
    for(i=0; i<nint; i++)
        ind->cell_interior[i]=cell_int[i];

    ind->cell_boundary = (unsigned int *) malloc(sizeof(unsigned int)*nbnd);
    for(i=0; i<nbnd; i++)
        ind->cell_boundary[i]=cell_bnd[i];

    ind->cell_fsi = (unsigned int *) malloc(sizeof(unsigned int)*nfsi);
    for(i=0; i<nfsi; i++)
        ind->cell_fsi[i]=cell_fsi[i];

    ind->cell_neumann= (unsigned int *) malloc(sizeof(unsigned int)*nneu);
    for(i=0; i<nneu; i++)
        ind->cell_neumann[i]=cell_neu[i];

    ind->cell_dirichlet = (unsigned int *) malloc(sizeof(unsigned int)*ndir);
    for(i=0; i<ndir; i++)
        ind->cell_dirichlet[i]=cell_dir[i];

	// part III: List the indices of the involved variables

	ind->xix_N = sumarray(xix_involved,N);
	ind->xix.g2G = (short int *) malloc(sizeof(short int)*ind->xix_N);
	ind->xix.G2g = (short int *) malloc(sizeof(short int)*N);

	ind->xiy_N = sumarray(xiy_involved,N);
	ind->xiy.g2G = (short int *) malloc(sizeof(short int)*ind->xiy_N);
	ind->xiy.G2g = (short int *) malloc(sizeof(short int)*N);

	ind->xix_Ncell_FSI  = sumarray(cell_involved_fsi_x,N);
	ind->xix.c2C_fsi = (short int *) malloc(sizeof(short int)*ind->xix_Ncell_FSI);
	ind->xix.C2c_fsi = (short int *) malloc(sizeof(short int)*N);

	ind->xiy_Ncell_FSI  = sumarray(cell_involved_fsi_y,N);
	ind->xiy.c2C_fsi = (short int *) malloc(sizeof(short int)*ind->xiy_Ncell_FSI);
	ind->xiy.C2c_fsi = (short int *) malloc(sizeof(short int)*N);

	ind->xix_Ncell_Neumann  = sumarray(cell_involved_neu_x,N);
	ind->xix.c2C_neumann = (short int *) malloc(sizeof(short int)*ind->xix_Ncell_Neumann);
	ind->xix.C2c_neumann = (short int *) malloc(sizeof(short int)*N);

	ind->xiy_Ncell_Neumann  = sumarray(cell_involved_neu_y,N);
	ind->xiy.c2C_neumann = (short int *) malloc(sizeof(short int)*ind->xiy_Ncell_Neumann);
	ind->xiy.C2c_neumann = (short int *) malloc(sizeof(short int)*N);

	ind->xix_Ncell_Dirichlet  = sumarray(cell_involved_dir_x,N);
	ind->xix.c2C_dirichlet = (short int *) malloc(sizeof(short int)*ind->xix_Ncell_Dirichlet);
	ind->xix.C2c_dirichlet = (short int *) malloc(sizeof(short int)*N);

	ind->xiy_Ncell_Dirichlet  = sumarray(cell_involved_dir_y,N);
	ind->xiy.c2C_dirichlet = (short int *) malloc(sizeof(short int)*ind->xiy_Ncell_Dirichlet);
	ind->xiy.C2c_dirichlet = (short int *) malloc(sizeof(short int)*N);

	for(i=0; i<N; i++)
	{
		if (xix_involved[i]==1)
		{
			ind->xix.g2G[n_xix] =(short int) i;
			ind->xix.G2g[i] =(short int) n_xix;
			n_xix +=1;
		}
		else
			ind->xix.G2g[i] =(short int) -1;

		if (xiy_involved[i]==1)
		{
			ind->xiy.g2G[n_xiy] =(short int) i;
			ind->xiy.G2g[i] =(short int) n_xiy;
			n_xiy +=1;
		}
		else
			ind->xiy.G2g[i] =(short int) -1;

		if (cell_involved_fsi_x[i]==1)
		{
			ind->xix.c2C_fsi[n_xfsi] = (short int)i;
			ind->xix.C2c_fsi[i] =(short int) n_xfsi;
			n_xfsi +=1;
		}
		else
			ind->xix.C2c_fsi[i] =(short int) -1;

		if (cell_involved_fsi_y[i]==1)
		{
			ind->xiy.c2C_fsi[n_yfsi] = (short int)i;
			ind->xiy.C2c_fsi[i] = (short int)n_yfsi;
			n_yfsi +=1;
		}
		else
			ind->xiy.C2c_fsi[i] =(short int) -1;

		if (cell_involved_neu_x[i]==1)
		{
			ind->xix.c2C_neumann[n_xneu] =(short int) i;
			ind->xix.C2c_neumann[i] = (short int)n_xneu;
			n_xneu +=1;
		}
		else
			ind->xix.C2c_neumann[i] =(short int) -1;

		if (cell_involved_neu_y[i]==1)
		{
			ind->xiy.c2C_neumann[n_yneu] =(short int) i;
			ind->xiy.C2c_neumann[i] = (short int)n_yneu;
			n_yneu +=1;
		}
		else
			ind->xiy.C2c_neumann[i] =(short int) -1;

		if (cell_involved_dir_x[i]==1)
		{
			ind->xix.c2C_dirichlet[n_xdir] =(short int) i;
			ind->xix.C2c_dirichlet[i] =(short int) n_xdir;
			n_xdir +=1;
		}
		else
			ind->xix.C2c_dirichlet[i] =(short int) -1;

		if (cell_involved_dir_y[i]==1)
		{
			ind->xiy.c2C_dirichlet[n_ydir] =(short int) i;
			ind->xiy.C2c_dirichlet[i] = (short int)n_ydir;
			n_ydir +=1;
		}
		else
			ind->xiy.C2c_dirichlet[i] = (short int)-1;

	}

}

struct invol involvedIndices_grid(int cell, int m){

	struct invol involved;
	int i,j;
	int kk,im1,ip1,jm1,jp1;
	int im1_jp1,ip1_jm1,ip1_jp1;

	i = cell%m;
	j = (int) floor((double)cell/m);

	kk	= i   +	j*m;
	im1	= i-1 + j*m;
	ip1	= i+1 + j*m;
	jm1	= i   + (j-1)*m;
	jp1	= i   + (j+1)*m;

	im1_jp1 = i-1   +	(j+1)*m;
	ip1_jm1 = i+1   +	(j-1)*m;
	ip1_jp1 = i+1   +	(j+1)*m;

    involved.i = i;
    involved.j = j;

	involved.xix[0] = jm1;
	involved.xix[1] = ip1_jm1;
	involved.xix[2] = kk;
	involved.xix[3] = ip1;
	involved.xix[4] = jp1;
	involved.xix[5] = ip1_jp1;

	involved.xiy[0] = im1;
	involved.xiy[1] = kk;
	involved.xiy[2] = ip1;
	involved.xiy[3] = im1_jp1;
	involved.xiy[4] = jp1;
	involved.xiy[5] = ip1_jp1;

	involved.stagx_cell[0] = jm1;
	involved.stagx_cell[1] = kk;

	involved.stagy_cell[0] = im1;
	involved.stagy_cell[1] = kk;

	return involved;
}

// summing up all the element in the array
int sumarray (int* a, int n){

	int i;
	int sum = 0;

	for(i=0; i<n; i++)
		sum += a[i];

	return sum;
}
