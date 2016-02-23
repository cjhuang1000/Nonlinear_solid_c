#include "init_setting.h"
#include <math.h>
#include "stdlib.h"

int sumarray (int array[], int target, int n);

void set_index(Field_S* s, AppCtx* ptr_u)
{

	int N = s->N,	m = s->Nx,	n = s->Ny;
	int nint	=0,	nbnd	=0,		nfsi	=0, nneu	=0, ndir	=0;
	int n_xix	=0,	n_xiy	=0;
	int n_xfsi	=0,	n_yfsi	=0;
	int n_xneu	=0,	n_yneu	=0;
	int n_xdir	=0,	n_ydir	=0;
	int ns_xixy	,	temp; // starting index of each processors
	int ns_xyfsi;
	int ns_xyneu;
	int ns_xydir;
	int i,j,k,kk;

	unsigned int 	*cell_int,	*cell_bnd, *cell_fsi, *cell_neu, *cell_dir;
	int 			*xix_involved,	*xiy_involved;
	int 			*cell_involved_fsi_x,	*cell_involved_fsi_y;
	int 			*cell_involved_neu_x,	*cell_involved_neu_y;
	int 			*cell_involved_dir_x,	*cell_involved_dir_y;
	int 			**bs = s->param.boundary_sign;

	struct 			invol involved;

	PetscInt 		temp2;
	PetscInt 		*pordering_x, *aordering_x, *pordering_y, *aordering_y;
	AO		 		ao_x, ao_y;
	Index_S			*ind = &(s->ind);

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


	s->dmx_s = 1e9;
	s->dmx_e = 0;
	s->dmy_s = 1e9;
	s->dmy_e = 0;


	for(j=1; j<n-1; j++)
		for(i=1; i<m-1;i++)
		{
			// Not-external cells (at least 1 corner in)
			if ((bs[i][j] == -1) || (bs[i+1][j] == -1) ||
				(bs[i][j+1] == -1) || (bs[i+1][j+1] == -1))
			{
				kk = i+j*m;
				involved    = involvedIndices_grid( kk , m);
				for(k=0; k<6; k++)
				{
					if (s->v2p[involved.xix[k]] == ptr_u->rank) // if the cell belongs to this processor
						xix_involved[involved.xix[k]] =1;
					else if (s->v2p[kk] == ptr_u->rank) // if the cell belongs to the next processor
						xix_involved[involved.xix[k]] =2;

					if (s->v2p[involved.xiy[k]] == ptr_u->rank)
						xiy_involved[involved.xiy[k]] =1;
					else if (s->v2p[kk] == ptr_u->rank)
						xiy_involved[involved.xiy[k]] =2;
				}

				// interior
				if ((bs[i][j]+bs[i+1][j]+bs[i][j+1]+bs[i+1][j+1]) <= -3)
				{
					if (s->v2p[kk] == ptr_u->rank)
					{
						cell_int[nint] = kk;
						nint +=1;
						if (i < s->dmx_s) s->dmx_s = i;
						if (i > s->dmx_e) s->dmx_e = i;
						if (j < s->dmy_s) s->dmy_s = j;
						if (j > s->dmy_e) s->dmy_e = j;
					}
				}
				// boundary cells
				else{

					if (s->v2p[kk] == ptr_u->rank)
					{
						cell_bnd[nbnd] = kk;
						nbnd +=1;
						if (i < s->dmx_s) s->dmx_s = i;
						if (i > s->dmx_e) s->dmx_e = i;
						if (j < s->dmy_s) s->dmy_s = j;
						if (j > s->dmy_e) s->dmy_e = j;
					}

					switch (s->con.fsineumanndirichlet[kk])
					{
					case 'F':

						if (s->v2p[kk] == ptr_u->rank)
						{
							cell_fsi[nfsi] = kk;
							nfsi +=1;
						}

						for(k=0; k<2; k++)
						{
							if (s->v2p[involved.stagx_cell[k]] == ptr_u->rank)
								cell_involved_fsi_x[involved.stagx_cell[k]] =1;
							else if (s->v2p[kk] == ptr_u->rank) // if the cell belongs to the next processor
								cell_involved_fsi_x[involved.stagx_cell[k]] =2;

							if (s->v2p[involved.stagy_cell[k]] == ptr_u->rank)
								cell_involved_fsi_y[involved.stagy_cell[k]] =1;
							else if (s->v2p[kk] == ptr_u->rank)
								cell_involved_fsi_y[involved.stagy_cell[k]] =2;
						}
						break;
					case 'N':

						if (s->v2p[kk] == ptr_u->rank)
						{
							cell_neu[nneu] = kk;
							nneu +=1;
						}

						for(k=0; k<2; k++)
						{
							if (s->v2p[involved.stagx_cell[k]] == ptr_u->rank)
								cell_involved_neu_x[involved.stagx_cell[k]] =1;
							else if (s->v2p[kk] == ptr_u->rank)
								cell_involved_neu_x[involved.stagx_cell[k]] =2;

							if (s->v2p[involved.stagy_cell[k]] == ptr_u->rank)
								cell_involved_neu_y[involved.stagy_cell[k]] =1;
							else if (s->v2p[kk] == ptr_u->rank)
								cell_involved_neu_y[involved.stagy_cell[k]] =2;
						}
						break;
					case 'D':

						if (s->v2p[kk] == ptr_u->rank)
						{
							cell_dir[ndir] = kk;
							ndir +=1;
						}

						for(k=0; k<2; k++)
						{
							if (s->v2p[involved.stagx_cell[k]] == ptr_u->rank)
								cell_involved_dir_x[involved.stagx_cell[k]] =1;
							else if (s->v2p[kk] == ptr_u->rank)
								cell_involved_dir_x[involved.stagx_cell[k]] =2;

							if (s->v2p[involved.stagy_cell[k]] == ptr_u->rank)
								cell_involved_dir_y[involved.stagy_cell[k]] =1;
							else if (s->v2p[kk] == ptr_u->rank)
								cell_involved_dir_y[involved.stagy_cell[k]] =2;
						}
						break;

					}
				}

			}
		}

	// part II: Count the number of involved cells

	ind->cell_N_interior	= nint;
	ind->cell_N_boundary	= nbnd;
	ind->cell_N_Neumann		= nneu;
	ind->cell_N_Dirichlet	= ndir;
	ind->cell_N_FSI			= nfsi;

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

    // to here
	// part II-I: reordering
    ind->xix_N = sumarray(xix_involved,1,N);
    ind->xiy_N = sumarray(xiy_involved,1,N);
    ind->xix_Ncell_FSI  = sumarray(cell_involved_fsi_x,1,N);
    ind->xiy_Ncell_FSI  = sumarray(cell_involved_fsi_y,1,N);
    ind->xix_Ncell_Neumann  = sumarray(cell_involved_neu_x,1,N);
    ind->xiy_Ncell_Neumann  = sumarray(cell_involved_neu_y,1,N);
    ind->xix_Ncell_Dirichlet  = sumarray(cell_involved_dir_x,1,N);
    ind->xiy_Ncell_Dirichlet  = sumarray(cell_involved_dir_y,1,N);

    ind->xix_ghoN = sumarray(xix_involved,2,N);
    ind->xiy_ghoN = sumarray(xiy_involved,2,N);

    temp = ind->xix_N + ind->xiy_N;
    MPI_Scan(&temp,&ns_xixy,1,MPIU_INT,MPI_SUM,PETSC_COMM_WORLD); // get the starting index of xi x,y on each processor
    ns_xixy -= temp;

    temp = ind->xix_Ncell_FSI + ind->xiy_Ncell_FSI;
    MPI_Scan(&temp,&ns_xyfsi,1,MPIU_INT,MPI_SUM,PETSC_COMM_WORLD);
    ns_xyfsi -= temp;

    temp = ind->xix_Ncell_Neumann + ind->xiy_Ncell_Neumann;
    MPI_Scan(&temp,&ns_xyneu,1,MPIU_INT,MPI_SUM,PETSC_COMM_WORLD);
    ns_xyneu -= temp;

    temp = ind->xix_Ncell_Dirichlet + ind->xiy_Ncell_Dirichlet;
    MPI_Scan(&temp,&ns_xydir,1,MPIU_INT,MPI_SUM,PETSC_COMM_WORLD);
    ns_xydir -= temp;

    // part III: List the indices of the involved variables

	ind->xix.l2g = (short int *) malloc(sizeof(short int)*(ind->xix_N + ind->xix_ghoN));
	ind->xix.l2G = (short int *) malloc(sizeof(short int)*(ind->xix_N + ind->xix_ghoN));
	ind->xix.G2g = (short int *) malloc(sizeof(short int)*N);

	ind->xiy.l2g= (short int *) malloc(sizeof(short int)*(ind->xiy_N + ind->xiy_ghoN));
	ind->xiy.l2G= (short int *) malloc(sizeof(short int)*(ind->xiy_N + ind->xiy_ghoN));
	ind->xiy.G2g = (short int *) malloc(sizeof(short int)*N);

	ind->xix.c2C_fsi = (short int *) malloc(sizeof(short int)*ind->xix_Ncell_FSI);
	ind->xix.C2c_fsi = (short int *) malloc(sizeof(short int)*N);

	ind->xiy.c2C_fsi = (short int *) malloc(sizeof(short int)*ind->xiy_Ncell_FSI);
	ind->xiy.C2c_fsi = (short int *) malloc(sizeof(short int)*N);

	ind->xix.c2C_neumann = (short int *) malloc(sizeof(short int)*ind->xix_Ncell_Neumann);
	ind->xix.C2c_neumann = (short int *) malloc(sizeof(short int)*N);

	ind->xiy.c2C_neumann = (short int *) malloc(sizeof(short int)*ind->xiy_Ncell_Neumann);
	ind->xiy.C2c_neumann = (short int *) malloc(sizeof(short int)*N);

	ind->xix.c2C_dirichlet = (short int *) malloc(sizeof(short int)*ind->xix_Ncell_Dirichlet);
	ind->xix.C2c_dirichlet = (short int *) malloc(sizeof(short int)*N);

	ind->xiy.c2C_dirichlet = (short int *) malloc(sizeof(short int)*ind->xiy_Ncell_Dirichlet);
	ind->xiy.C2c_dirichlet = (short int *) malloc(sizeof(short int)*N);


	for(i=0; i<N; i++)
	{
		if (xix_involved[i]==1)
		{
			ind->xix.l2G[n_xix] =(short int) i;
			ind->xix.l2g[n_xix] =(short int) n_xix + ns_xixy;
			ind->xix.G2g[i] =(short int) n_xix + ns_xixy;

			n_xix +=1;
		}
		else
			ind->xix.G2g[i] =(short int) -1;

		if (xiy_involved[i]==1)
		{
			ind->xiy.l2G[n_xiy] =(short int) i;
			ind->xiy.l2g[n_xiy] =(short int) n_xiy + ns_xixy + ind->xix_N;
			ind->xiy.G2g[i] =(short int) n_xiy + ns_xixy + ind->xix_N;
			n_xiy +=1;
		}
		else
			ind->xiy.G2g[i] =(short int) -1;

		if (cell_involved_fsi_x[i]==1)
		{
			ind->xix.c2C_fsi[n_xfsi] = (short int)i;
			ind->xix.C2c_fsi[i] =(short int) n_xfsi + ns_xyfsi;
			n_xfsi +=1;
		}
		else
			ind->xix.C2c_fsi[i] =(short int) -1;

		if (cell_involved_fsi_y[i]==1)
		{
			ind->xiy.c2C_fsi[n_yfsi] = (short int)i;
			ind->xiy.C2c_fsi[i] = (short int)n_yfsi + ns_xyfsi + ind->xix_Ncell_FSI;
			n_yfsi +=1;
		}
		else
			ind->xiy.C2c_fsi[i] =(short int) -1;

		if (cell_involved_neu_x[i]==1)
		{
			ind->xix.c2C_neumann[n_xneu] =(short int) i;
			ind->xix.C2c_neumann[i] = (short int)n_xneu + ns_xyneu;
			n_xneu +=1;
		}
		else
			ind->xix.C2c_neumann[i] =(short int) -1;

		if (cell_involved_neu_y[i]==1)
		{
			ind->xiy.c2C_neumann[n_yneu] =(short int) i;
			ind->xiy.C2c_neumann[i] = (short int)n_yneu + ns_xyneu + ind->xix_Ncell_Neumann;
			n_yneu +=1;
		}
		else
			ind->xiy.C2c_neumann[i] =(short int) -1;

		if (cell_involved_dir_x[i]==1)
		{
			ind->xix.c2C_dirichlet[n_xdir] =(short int) i;
			ind->xix.C2c_dirichlet[i] =(short int) n_xdir + ns_xydir;
			n_xdir +=1;
		}
		else
			ind->xix.C2c_dirichlet[i] =(short int) -1;

		if (cell_involved_dir_y[i]==1)
		{
			ind->xiy.c2C_dirichlet[n_ydir] =(short int) i;
			ind->xiy.C2c_dirichlet[i] = (short int)n_ydir + ns_xydir + ind->xix_Ncell_Dirichlet;
			n_ydir +=1;
		}
		else
			ind->xiy.C2c_dirichlet[i] = (short int)-1;

	}

	// part III-II: calculate global numbers

	MPI_Reduce(&ind->xix_N, &ind->xix_gloN, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ind->xix_gloN, 1, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Reduce(&ind->xiy_N, &ind->xi_gloN, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ind->xi_gloN, 1, MPI_INT, 0, MPI_COMM_WORLD);
	ind->xi_gloN += ind->xix_gloN;

	MPI_Reduce(&ind->xix_Ncell_Neumann, &ind->xix_gloNcell_Neumann, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ind->xix_gloNcell_Neumann, 1, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Reduce(&ind->xiy_Ncell_Neumann, &ind->xi_gloNcell_Neumann, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ind->xi_gloNcell_Neumann, 1, MPI_INT, 0, MPI_COMM_WORLD);
	ind->xi_gloNcell_Neumann += ind->xix_gloNcell_Neumann;

	MPI_Reduce(&ind->xix_Ncell_Dirichlet, &ind->xix_gloNcell_Dirichlet, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ind->xix_gloNcell_Dirichlet, 1, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Reduce(&ind->xiy_Ncell_Dirichlet, &ind->xi_gloNcell_Dirichlet, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ind->xi_gloNcell_Dirichlet, 1, MPI_INT, 0, MPI_COMM_WORLD);
	ind->xi_gloNcell_Dirichlet += ind->xix_gloNcell_Dirichlet;

	//printf("processor %d shows %d\n", ptr_u->rank, ind->xix_gloN);
	// for fsi (to be done)

	// part IV: Ghost cells
	// xix and xiy: append ghost cells at the end of l2G and l2g

	n_xix = 0; n_xiy = 0;

	PetscMalloc1(ind->xix_N, &pordering_x);
	PetscMalloc1(ind->xix_N, &aordering_x);
	PetscMalloc1(ind->xiy_N, &pordering_y);
	PetscMalloc1(ind->xiy_N, &aordering_y);

	for(i=0; i<ind->xix_N; i++)
	{
		pordering_x[i] = ind->xix.l2g[i];
		aordering_x[i] = ind->xix.l2G[i];
	}

	for(i=0; i<ind->xiy_N; i++)
	{
		pordering_y[i] = ind->xiy.l2g[i];
		aordering_y[i] = ind->xiy.l2G[i];
	}


	AOCreateMapping(MPI_COMM_WORLD,ind->xix_N,aordering_x,pordering_x,&ao_x);
	AOCreateMapping(MPI_COMM_WORLD,ind->xiy_N,aordering_y,pordering_y,&ao_y);

	//AOView(ao_x,PETSC_VIEWER_STDOUT_WORLD);
	//AOView(ao_y,PETSC_VIEWER_STDOUT_WORLD);

	// append ghost cell at the end of g2G and
	for(i=0; i<N; i++)
	{
		if (xix_involved[i]==2)
		{
			ind->xix.l2G[n_xix + ind->xix_N] =(short int) i;
			temp2 = i;
			AOApplicationToPetsc(ao_x,1,&temp2);
			ind->xix.G2g[i] =(short int) temp2;
			ind->xix.l2g[n_xix + ind->xix_N] =(short int) temp2;
			n_xix +=1;
		}
		if (xiy_involved[i]==2)
		{
			ind->xiy.l2G[n_xiy + ind->xiy_N] =(short int) i;
			temp2 = i;
			AOApplicationToPetsc(ao_y,1,&temp2);
			ind->xiy.G2g[i] =(short int) temp2;
			ind->xiy.l2g[n_xiy + ind->xiy_N] =(short int) temp2;
			n_xiy +=1;
		}
	}

	PetscFree(pordering_x);  PetscFree(aordering_x); PetscFree(pordering_y); PetscFree(aordering_y);
	AODestroy(&ao_x);
	AODestroy(&ao_y);



	// boundary cells - Neumann: do not need to append the ghost cells, only update the C2c

	PetscMalloc1(ind->xix_Ncell_Neumann, &pordering_x);
	PetscMalloc1(ind->xix_Ncell_Neumann, &aordering_x);
	PetscMalloc1(ind->xiy_Ncell_Neumann, &pordering_y);
	PetscMalloc1(ind->xiy_Ncell_Neumann, &aordering_y);

	for(i=0; i<ind->xix_Ncell_Neumann; i++)
	{
		pordering_x[i] = ind->xix.C2c_neumann[ind->xix.c2C_neumann[i]];
		aordering_x[i] = ind->xix.c2C_neumann[i];
	}

	for(i=0; i<ind->xiy_Ncell_Neumann; i++)
	{
		pordering_y[i] = ind->xiy.C2c_neumann[ind->xiy.c2C_neumann[i]];
		aordering_y[i] = ind->xiy.c2C_neumann[i];
	}

	AOCreateMapping(MPI_COMM_WORLD,ind->xix_Ncell_Neumann,aordering_x,pordering_x,&ao_x);
	AOCreateMapping(MPI_COMM_WORLD,ind->xiy_Ncell_Neumann,aordering_y,pordering_y,&ao_y);

	for(i=0; i<N; i++)
	{
		if (cell_involved_neu_x[i]==2)
		{
			temp2 = i;
			AOApplicationToPetsc(ao_x,1,&temp2);
			ind->xix.C2c_neumann[i] =(short int) temp2;
		}
		if (cell_involved_neu_y[i]==2)
		{
			temp2 = i;
			AOApplicationToPetsc(ao_y,1,&temp2);
			ind->xiy.C2c_neumann[i] =(short int) temp2;
		}
	}
	PetscFree(pordering_x);  PetscFree(aordering_x); PetscFree(pordering_y); PetscFree(aordering_y);
	AODestroy(&ao_x);
	AODestroy(&ao_y);

	// boundary cells - Dirichlet: do not need to append the ghost cells, only update the C2c

	PetscMalloc1(ind->xix_Ncell_Dirichlet, &pordering_x);
	PetscMalloc1(ind->xix_Ncell_Dirichlet, &aordering_x);
	PetscMalloc1(ind->xiy_Ncell_Dirichlet, &pordering_y);
	PetscMalloc1(ind->xiy_Ncell_Dirichlet, &aordering_y);

	for(i=0; i<ind->xix_Ncell_Dirichlet; i++)
	{
		pordering_x[i] = ind->xix.C2c_dirichlet[ind->xix.c2C_dirichlet[i]];
		aordering_x[i] = ind->xix.c2C_dirichlet[i];
	}

	for(i=0; i<ind->xiy_Ncell_Dirichlet; i++)
	{
		pordering_y[i] = ind->xiy.C2c_dirichlet[ind->xiy.c2C_dirichlet[i]];
		aordering_y[i] = ind->xiy.c2C_dirichlet[i];
	}

	AOCreateMapping(MPI_COMM_WORLD,ind->xix_Ncell_Dirichlet,aordering_x,pordering_x,&ao_x);
	AOCreateMapping(MPI_COMM_WORLD,ind->xiy_Ncell_Dirichlet,aordering_y,pordering_y,&ao_y);

	for(i=0; i<N; i++)
	{
		if (cell_involved_dir_x[i]==2)
		{
			temp2 = i;
			AOApplicationToPetsc(ao_x,1,&temp2);
			ind->xix.C2c_dirichlet[i] =(short int) temp2;
		}
		if (cell_involved_dir_y[i]==2)
		{
			temp2 = i;
			AOApplicationToPetsc(ao_y,1,&temp2);
			ind->xiy.C2c_dirichlet[i] =(short int) temp2;
		}
	}
	PetscFree(pordering_x);  PetscFree(aordering_x); PetscFree(pordering_y); PetscFree(aordering_y);
	AODestroy(&ao_x);
	AODestroy(&ao_y);

	free(cell_int);		free(cell_bnd);  free(cell_fsi); free(cell_neu);  free(cell_dir);
	free(xix_involved);	free(xiy_involved);
	free(cell_involved_fsi_x);	free(cell_involved_fsi_y);
	free(cell_involved_neu_x);	free(cell_involved_neu_y);
	free(cell_involved_dir_x);	free(cell_involved_dir_y);

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

// summing up all the ones in the array
int sumarray (int* a,int target, int n){

	int i;
	int sum = 0;

	for(i=0; i<n; i++)
		if (a[i] == target)
			sum ++;

	return sum;
}
