/*
 * interface_marker.c
 *
 * from matlab script:
 *
 * set_marker: set up the locations of the markers
 * ImmersedBoundary_info_cavity: update the sdf
 * force_calculation: compute the force distribution using the
 * information in the fluid domain
 *
 *  Created on: Oct 1, 2015
 *      Author: peggyhuang
 */

#include "interface_marker.h"
#include "dgefa3.c"

// offset_table.c
// ---------------------------
const offset2d_t offset2d[8] =
  {
    {-1,-1},{0,-1},{1,-1},
    {-1, 0}       ,{1, 0},
    {-1, 1},{0, 1},{1, 1}
  };

const offset3d_t offset3d[26] =
  {
    {-1,-1,-1},{0,-1,-1},{1,-1,-1},
    {-1, 0,-1},{0, 0,-1},{1, 0,-1},
    {-1, 1,-1},{0, 1,-1},{1, 1,-1},

    {-1,-1, 0},{0,-1, 0},{1,-1, 0},
    {-1, 0, 0},          {1, 0, 0},
    {-1, 1, 0},{0, 1, 0},{1, 1, 0},

    {-1,-1, 1},{0,-1, 1},{1,-1, 1},
    {-1, 0, 1},{0, 0, 1},{1, 0, 1},
    {-1, 1, 1},{0, 1, 1},{1, 1, 1},
  };

// ---------------------------


static int forcing_pts_list_init(forcing_point_list_t *list);
static int forcing_pts_list_addto(forcing_point_list_t *list,
                                  int i, int j, int k);

void marker_setup (Grid_S* ptr_g, Solid* ptr_s, Index_S* ptr_i,Lag_marker* ptr_mk, AppCtx* ptr_u)
{
	const int  	MAX_N = ptr_g->N;  // maximum number of markers
	double      sdf[9],x[9],y[9];
	double		n[2],t[2];
	int			index = 0, temp[4][4] = {{0,1,4,3},{1,2,5,4},{3,4,7,6},{4,5,8,7}};

	struct bound_elem{
		double sx,sy,ex,ey;
		int    tx,ty;
	}lg_marker[MAX_N];

	struct subgrid_info{
		double boundary_value[4];
		double x[4];
		double y[4];
	}symmetry[4];

	int     	bound_sign[4],side,side_next,ncorner1 =4,ncorner2;
	double  	corners[4][8];
	double  	new_corner[2];
	int 		i,j,cell,subcell,ii,jj,cell_i;

	PetscScalar *sx_array,*sy_array,*ex_array,*ey_array,*tx_array,*ty_array,*l0;
	Vec			tx_vec, ty_vec,l0_parallel;

	// ====== part I: Find intersection between interface and grid lines =======

	for(cell_i=0;cell_i<ptr_i->cell_N_boundary; cell_i++)
	{
		cell = ptr_i->cell_boundary[cell_i];

		ii = cell%ptr_g->Nx;
		jj = (int) floor((double)cell/ptr_g->Nx);

		// set the sdf of the doubly refined grid
		sdf[0] = ptr_s->boundary_value[ii][jj];
		sdf[2] = ptr_s->boundary_value[ii+1][jj];
		sdf[6] = ptr_s->boundary_value[ii][jj+1];
		sdf[8] = ptr_s->boundary_value[ii+1][jj+1];
		sdf[1] = (sdf[0]+sdf[2])/2;
		sdf[3] = (sdf[0]+sdf[6])/2;
		sdf[5] = (sdf[2]+sdf[8])/2;
		sdf[7] = (sdf[6]+sdf[8])/2;
		sdf[4] = (sdf[0]+sdf[2]+sdf[6]+sdf[8])/4;

		// normal vector (non-normalized) calculated by the SDF
		n[0]=sdf[5]-sdf[3]; n[1]=sdf[7]-sdf[1];

		// x- y- coordinates onf the doubly refined grid
		x[0] = ptr_g->x_grid[ii];
		x[2] = ptr_g->x_grid[ii+1];
		x[6] = ptr_g->x_grid[ii];
		x[8] = ptr_g->x_grid[ii+1];
		x[1] = (x[0]+x[2])/2;
		x[3] = (x[0]+x[6])/2;
		x[5] = (x[2]+x[8])/2;
		x[7] = (x[6]+x[8])/2;
		x[4] = (x[0]+x[2]+x[6]+x[8])/4;

		y[0] = ptr_g->y_grid[jj];
		y[2] = ptr_g->y_grid[jj];
		y[6] = ptr_g->y_grid[jj+1];
		y[8] = ptr_g->y_grid[jj+1];
		y[1] = (y[0]+y[2])/2;
		y[3] = (y[0]+y[6])/2;
		y[5] = (y[2]+y[8])/2;
		y[7] = (y[6]+y[8])/2;
		y[4] = (y[0]+y[2]+y[6]+y[8])/4;

		for(i=0; i<4; i++)
			for(j=0; j<4; j++)
			{
				symmetry[i].boundary_value[j]=sdf[temp[i][j]];
				symmetry[i].x[j]=x[temp[i][j]];
				symmetry[i].y[j]=y[temp[i][j]];
			}

		for(subcell=0;subcell<4; subcell++)
		{
			ncorner1 =4;

			for(i=0;i<4;i++)
			    bound_sign[i]=(int)SIGN(1,symmetry[subcell].boundary_value[i]);

			if( ((bound_sign[0]== -1) || (bound_sign[1]== -1) || (bound_sign[2]== -1) || (bound_sign[3]== -1))
			    && ( (bound_sign[0]+bound_sign[1]+bound_sign[2]+bound_sign[3]) > -3 ))
			{

			    for(i=0;i<4;i++)
			    {
			    	corners[0][i]=symmetry[subcell].x[i];
			    	corners[1][i]=symmetry[subcell].y[i];
			    	corners[2][i]=symmetry[subcell].boundary_value[i];
			    	corners[3][i]=i;
			    }

			    // Compute the vertices of subcells cut by the boundary
				for(side=3;side>-1;side--)
				{
				    side_next = (side+1)%4;
				    if( corners[2][side]*corners[2][side_next]<0)
				    {
				        ncorner1 +=1;

				        for(i=0;i<2;i++)
				            new_corner[i] = corners[i][side]+(corners[i][side_next]-corners[i][side])
				                        *corners[2][side]/(corners[2][side]-corners[2][side_next]);

				        for(i=0; i<4;i++)
				            for(j=ncorner1-1; j>side+1;j--)
				                corners[i][j] = corners[i][j-1];

				        corners[0][side+1] = new_corner[0];
				        corners[1][side+1] = new_corner[1];
				        corners[2][side+1] = 0;
				        corners[3][side+1] = side;
				    }
				}
				ncorner2 = ncorner1;

				// Delete the external corners
				for(side=ncorner1-1;side>-1;side--)
				{
				    if (corners[2][side] >0)
				    {
				        ncorner2 -= 1;
				        for(j=side; j<ncorner2;j++)
				            for(i=0; i<4; i++)
				                corners[i][j] = corners[i][j+1];
				    }
				}

				for(side=0; side<ncorner2;side++)
				{

					side_next = (side+1)%ncorner2;
				    if ((corners[2][side] == 0) && (corners[2][side_next] == 0))
				    {
				    	t[0]=corners[0][side_next]-corners[0][side]; // what am i writing about?
				    	t[1]=corners[1][side_next]-corners[1][side];

				    	if ((n[0]*t[1]-n[1]*t[0])>0 ) // n x t >0
				    	{
					   		lg_marker[index].sx = corners[0][side];
					   		lg_marker[index].ex = corners[0][side_next];
					   		lg_marker[index].sy = corners[1][side];
					   		lg_marker[index].ey = corners[1][side_next];
				    	}
				    	else
				    	{
					    	lg_marker[index].sx = corners[0][side_next];
					    	lg_marker[index].ex = corners[0][side];
					    	lg_marker[index].sy = corners[1][side_next];
					    	lg_marker[index].ey = corners[1][side];
				    	}

				    	switch (subcell)
				    	{
				    	case 0:
					    	lg_marker[index].tx = ii+(jj-1)*ptr_g->Nx;
					    	lg_marker[index].ty = ii-1+jj*ptr_g->Nx;
				    		break;
			    		case 1:
					    	lg_marker[index].tx = ii+(jj-1)*ptr_g->Nx;
					    	lg_marker[index].ty = ii+jj*ptr_g->Nx;
				    		break;
				    	case 2:
					    	lg_marker[index].tx = ii+jj*ptr_g->Nx;
					    	lg_marker[index].ty = ii-1+jj*ptr_g->Nx;
				    		break;
				    	case 3:
					    	lg_marker[index].tx = ii+jj*ptr_g->Nx;
					    	lg_marker[index].ty = ii+jj*ptr_g->Nx;
				    		break;
				    	} // end switch

				    	index++;

				    }
				 }
			} // end subcell intersects with boundary
		} // end subcell

	} // end boundary elements

	// ======== part II: build parallel marker Vecs ==============

	ptr_mk->nlocal_s = index;
	MPI_Reduce(&ptr_mk->nlocal_s, &ptr_mk->nglobal, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ptr_mk->nglobal, 1, MPI_INT, 0, MPI_COMM_WORLD);

    VecCreate(PETSC_COMM_WORLD,&ptr_mk->sp_x0);
    VecSetSizes(ptr_mk->sp_x0,ptr_mk->nlocal_s,ptr_mk->nglobal);
    VecSetFromOptions(ptr_mk->sp_x0);

    VecDuplicate(ptr_mk->sp_x0, &ptr_mk->sp_y0);
    VecDuplicate(ptr_mk->sp_x0, &ptr_mk->ep_x0);
    VecDuplicate(ptr_mk->sp_x0, &ptr_mk->ep_y0);
    VecDuplicate(ptr_mk->sp_x0, &tx_vec);
    VecDuplicate(ptr_mk->sp_x0, &ty_vec);
    VecDuplicate(ptr_mk->sp_x0, &l0_parallel);

    ptr_mk->tx 		= ivector(0, ptr_mk->nglobal-1);
    ptr_mk->ty 		= ivector(0, ptr_mk->nglobal-1);
    ptr_mk->lratio 	= dvector(0, ptr_mk->nglobal-1);

    // input the values in the Vec's
    VecGetArray(ptr_mk->sp_x0,&sx_array);
    VecGetArray(ptr_mk->sp_y0,&sy_array);
    VecGetArray(ptr_mk->ep_x0,&ex_array);
    VecGetArray(ptr_mk->ep_y0,&ey_array);
    VecGetArray(tx_vec		 ,&tx_array);
    VecGetArray(ty_vec		 ,&ty_array);
    VecGetArray(l0_parallel	 ,&l0);

    for(i=0; i<ptr_mk->nlocal_s; i++)
    {
    	sx_array[i]=lg_marker[i].sx;
    	sy_array[i]=lg_marker[i].sy;
    	ex_array[i]=lg_marker[i].ex;
    	ey_array[i]=lg_marker[i].ey;
    	tx_array[i]=ptr_i->xix.C2c_neumann[lg_marker[i].tx];
    	ty_array[i]=ptr_i->xiy.C2c_neumann[lg_marker[i].ty];
    	l0[i] = sqrt((sx_array[i]-ex_array[i])*(sx_array[i]-ex_array[i])
    			+ (sy_array[i]-ey_array[i])*(sy_array[i]-ey_array[i]));
    	printf("the %d panel is %d %f and %d %f \n",i, lg_marker[i].tx,tx_array[i],lg_marker[i].ty,ty_array[i]);
    }

    VecRestoreArray(ptr_mk->sp_x0,&sx_array);
    VecRestoreArray(ptr_mk->sp_y0,&sy_array);
    VecRestoreArray(ptr_mk->ep_x0,&ex_array);
    VecRestoreArray(ptr_mk->ep_y0,&ey_array);
    VecRestoreArray(tx_vec		 ,&tx_array);
    VecRestoreArray(ty_vec		 ,&ty_array);
    VecRestoreArray(l0_parallel	 ,&l0);

    // create sequential vec storing all the updated marker locations
	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nglobal,&ptr_mk->sp_x);
	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nglobal,&ptr_mk->sp_y);
	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nglobal,&ptr_mk->ep_x);
	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nglobal,&ptr_mk->ep_y);

	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nglobal,&ptr_mk->sp_u);
	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nglobal,&ptr_mk->sp_v);
	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nglobal,&ptr_mk->ep_u);
	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nglobal,&ptr_mk->ep_v);

	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nglobal,&ptr_mk->sp_ax);
	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nglobal,&ptr_mk->sp_ay);
	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nglobal,&ptr_mk->ep_ax);
	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nglobal,&ptr_mk->ep_ay);

    IS 			is;
    Vec			tx_vec_seq,ty_vec_seq;

	ISCreateStride	(PETSC_COMM_SELF,ptr_mk->nglobal,0,1,&is);
	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nglobal,&tx_vec_seq);
	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nglobal,&ty_vec_seq);
	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nglobal,&ptr_mk->l0);
	VecScatterCreate(ptr_mk->sp_x0,is,ptr_mk->sp_x,is,&ptr_mk->loc2all);

	VecScatterBegin	(ptr_mk->loc2all,l0_parallel,ptr_mk->l0,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd	(ptr_mk->loc2all,l0_parallel,ptr_mk->l0,INSERT_VALUES,SCATTER_FORWARD);

	VecScatterBegin	(ptr_mk->loc2all,tx_vec,tx_vec_seq,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd	(ptr_mk->loc2all,tx_vec,tx_vec_seq,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterBegin	(ptr_mk->loc2all,ty_vec,ty_vec_seq,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd	(ptr_mk->loc2all,ty_vec,ty_vec_seq,INSERT_VALUES,SCATTER_FORWARD);

	VecGetArray		(tx_vec_seq,&tx_array);
    VecGetArray		(ty_vec_seq,&ty_array);

    for(i=0; i<ptr_mk->nglobal; i++){
    	ptr_mk->tx[i] = tx_array[i];
    	ptr_mk->ty[i] = ty_array[i];
    }

    VecRestoreArray(tx_vec_seq,&tx_array);
    VecRestoreArray(ty_vec_seq,&ty_array);

    VecDestroy (&tx_vec_seq);
    VecDestroy (&ty_vec_seq);
    VecDestroy (&tx_vec);
    VecDestroy (&ty_vec);
    ISDestroy  (&is);

    // prepare the list_f
    ptr_mk->list_f 	= ivector(0, ptr_mk->nglobal-1);
    ptr_mk->nx 		= dvector(0, ptr_mk->nglobal-1);
    ptr_mk->ny 		= dvector(0, ptr_mk->nglobal-1);

    printf("Marker creation completed. \n");
    printf("The processor %d has %d markers. \n",ptr_u->rank,ptr_mk->nlocal_s);

}

// update the location of the current marker and build list
void build_current_marker(Lag_marker* ptr_mk, Field_S* ptr_s, Index_S* ptr_i, AppCtx* ptr_u){

	PetscReal 	*xix,*xiy,*xi; // local displacement array
	PetscReal	*sp_x,*sp_y,*ep_x,*ep_y,*l0;
	PetscReal   *mkp;
	PetscScalar	c = 1.0;
	Vec		   	dis_sp_x,dis_sp_y,dis_ep_x,dis_ep_y; // displacement and new marker locations
	Vec			xix_local, xiy_local;

	int			i_mk, comp, i, elem;
	int			stencil;
	struct interp2d *interp;

	VecCreateSeq(PETSC_COMM_SELF,ptr_i->xix_N +ptr_i->xix_ghoN,&xix_local);
    VecCreateSeq(PETSC_COMM_SELF,ptr_i->xiy_N +ptr_i->xiy_ghoN,&xiy_local);

	// Create and copy original locations of markers
	VecDuplicate(ptr_mk->sp_x0,&dis_sp_x);
	VecDuplicate(ptr_mk->sp_x0,&dis_sp_y);
	VecDuplicate(ptr_mk->sp_x0,&dis_ep_x);
	VecDuplicate(ptr_mk->sp_x0,&dis_ep_y);

    for(elem = 0; elem <3; elem ++){ // displacement, velocities, acceleration

    	switch(elem)
    	{
    	case 0:
    		VecScatterBegin(ptr_u->scatter_x,ptr_s->xi,xix_local,INSERT_VALUES,SCATTER_FORWARD);
    		VecScatterEnd  (ptr_u->scatter_x,ptr_s->xi,xix_local,INSERT_VALUES,SCATTER_FORWARD);

    		VecScatterBegin(ptr_u->scatter_y,ptr_s->xi,xiy_local,INSERT_VALUES,SCATTER_FORWARD);
    		VecScatterEnd  (ptr_u->scatter_y,ptr_s->xi,xiy_local,INSERT_VALUES,SCATTER_FORWARD);
    		break;

    	case 1:
    		VecScatterBegin(ptr_u->scatter_x,ptr_s->dxi,xix_local,INSERT_VALUES,SCATTER_FORWARD);
    	    VecScatterEnd  (ptr_u->scatter_x,ptr_s->dxi,xix_local,INSERT_VALUES,SCATTER_FORWARD);

    	    VecScatterBegin(ptr_u->scatter_y,ptr_s->dxi,xiy_local,INSERT_VALUES,SCATTER_FORWARD);
    	    VecScatterEnd  (ptr_u->scatter_y,ptr_s->dxi,xiy_local,INSERT_VALUES,SCATTER_FORWARD);
    	    break;

    	case 2:
    		VecScatterBegin(ptr_u->scatter_x,ptr_s->dxi,xix_local,INSERT_VALUES,SCATTER_FORWARD);
    	    VecScatterEnd  (ptr_u->scatter_x,ptr_s->dxi,xix_local,INSERT_VALUES,SCATTER_FORWARD);

    	    VecScatterBegin(ptr_u->scatter_y,ptr_s->dxi,xiy_local,INSERT_VALUES,SCATTER_FORWARD);
    	    VecScatterEnd  (ptr_u->scatter_y,ptr_s->dxi,xiy_local,INSERT_VALUES,SCATTER_FORWARD);
    	    break;
    	}

    	// Get local displacement xi = [x_local; y_local]
    	VecGetArray	(xix_local,&xix);
    	VecGetArray	(xiy_local,&xiy);

    	VecGetArray(dis_sp_x,&sp_x);
    	VecGetArray(dis_sp_y,&sp_y);
    	VecGetArray(dis_ep_x,&ep_x);
    	VecGetArray(dis_ep_y,&ep_y);

		for(comp=0; comp<4; comp++){

			switch (comp)
			{
			case 0:
				mkp   = sp_x;
				interp= ptr_mk->interp_s_x;
				xi = xix;
				break;
			case 1:
				mkp   = sp_y;
				interp= ptr_mk->interp_s_y;
				xi = xiy;
				break;
			case 2:
				mkp   = ep_x;
				interp= ptr_mk->interp_e_x;
				xi = xix;
				break;
			case 3:
				mkp   = ep_y;
				interp= ptr_mk->interp_e_y;
				xi = xiy;
				break;
			}

			for(i_mk=0; i_mk<ptr_mk->nlocal_s;i_mk++){ // loop through local markers

				mkp[i_mk] =0.0;

				for(i=0; i<interp[i_mk].n; i++){

					stencil 	=  interp[i_mk].sten[i];
					mkp[i_mk] 	+= interp[i_mk].coeff[i]*xi[stencil];

				}
			}
		}

		VecRestoreArray(dis_sp_x,&sp_x);
		VecRestoreArray(dis_sp_y,&sp_y);
		VecRestoreArray(dis_ep_x,&ep_x);
		VecRestoreArray(dis_ep_y,&ep_y);

		VecRestoreArray(xix_local,&xix);
		VecRestoreArray(xiy_local,&xiy);

		switch(elem)
		{
		case 0:

			// New location of markers
			VecAXPY(dis_sp_x,c,ptr_mk->sp_x0);
			VecAXPY(dis_sp_y,c,ptr_mk->sp_y0);
			VecAXPY(dis_ep_x,c,ptr_mk->ep_x0);
			VecAXPY(dis_ep_y,c,ptr_mk->ep_y0);

			// Distribute the updated location to all the processors
			VecScatterBegin	(ptr_mk->loc2all,dis_sp_x,ptr_mk->sp_x,INSERT_VALUES,SCATTER_FORWARD);
			VecScatterEnd	(ptr_mk->loc2all,dis_sp_x,ptr_mk->sp_x,INSERT_VALUES,SCATTER_FORWARD);

			VecScatterBegin	(ptr_mk->loc2all,dis_sp_y,ptr_mk->sp_y,INSERT_VALUES,SCATTER_FORWARD);
			VecScatterEnd	(ptr_mk->loc2all,dis_sp_y,ptr_mk->sp_y,INSERT_VALUES,SCATTER_FORWARD);

			VecScatterBegin	(ptr_mk->loc2all,dis_ep_x,ptr_mk->ep_x,INSERT_VALUES,SCATTER_FORWARD);
			VecScatterEnd	(ptr_mk->loc2all,dis_ep_x,ptr_mk->ep_x,INSERT_VALUES,SCATTER_FORWARD);

			VecScatterBegin	(ptr_mk->loc2all,dis_ep_y,ptr_mk->ep_y,INSERT_VALUES,SCATTER_FORWARD);
			VecScatterEnd	(ptr_mk->loc2all,dis_ep_y,ptr_mk->ep_y,INSERT_VALUES,SCATTER_FORWARD);

			break;

		case 1:

			VecScatterBegin	(ptr_mk->loc2all,dis_sp_x,ptr_mk->sp_u,INSERT_VALUES,SCATTER_FORWARD);
			VecScatterEnd	(ptr_mk->loc2all,dis_sp_x,ptr_mk->sp_u,INSERT_VALUES,SCATTER_FORWARD);

			VecScatterBegin	(ptr_mk->loc2all,dis_sp_y,ptr_mk->sp_v,INSERT_VALUES,SCATTER_FORWARD);
			VecScatterEnd	(ptr_mk->loc2all,dis_sp_y,ptr_mk->sp_v,INSERT_VALUES,SCATTER_FORWARD);

			VecScatterBegin	(ptr_mk->loc2all,dis_ep_x,ptr_mk->ep_u,INSERT_VALUES,SCATTER_FORWARD);
			VecScatterEnd	(ptr_mk->loc2all,dis_ep_x,ptr_mk->ep_u,INSERT_VALUES,SCATTER_FORWARD);

			VecScatterBegin	(ptr_mk->loc2all,dis_ep_y,ptr_mk->ep_v,INSERT_VALUES,SCATTER_FORWARD);
			VecScatterEnd	(ptr_mk->loc2all,dis_ep_y,ptr_mk->ep_v,INSERT_VALUES,SCATTER_FORWARD);
			break;

		case 2:

			VecScatterBegin	(ptr_mk->loc2all,dis_sp_x,ptr_mk->sp_ax,INSERT_VALUES,SCATTER_FORWARD);
			VecScatterEnd	(ptr_mk->loc2all,dis_sp_x,ptr_mk->sp_ax,INSERT_VALUES,SCATTER_FORWARD);

			VecScatterBegin	(ptr_mk->loc2all,dis_sp_y,ptr_mk->sp_ay,INSERT_VALUES,SCATTER_FORWARD);
			VecScatterEnd	(ptr_mk->loc2all,dis_sp_y,ptr_mk->sp_ay,INSERT_VALUES,SCATTER_FORWARD);

			VecScatterBegin	(ptr_mk->loc2all,dis_ep_x,ptr_mk->ep_ax,INSERT_VALUES,SCATTER_FORWARD);
			VecScatterEnd	(ptr_mk->loc2all,dis_ep_x,ptr_mk->ep_ax,INSERT_VALUES,SCATTER_FORWARD);

			VecScatterBegin	(ptr_mk->loc2all,dis_ep_y,ptr_mk->ep_ay,INSERT_VALUES,SCATTER_FORWARD);
			VecScatterEnd	(ptr_mk->loc2all,dis_ep_y,ptr_mk->ep_ay,INSERT_VALUES,SCATTER_FORWARD);
			break;
		}
    }

	//get the ratio of the current length and initial length
	VecGetArray(ptr_mk->l0,&l0);

	VecGetArray(ptr_mk->sp_x,&sp_x);
	VecGetArray(ptr_mk->sp_y,&sp_y);
	VecGetArray(ptr_mk->ep_x,&ep_x);
	VecGetArray(ptr_mk->ep_y,&ep_y);

	for(i_mk = 0; i_mk < ptr_mk->nglobal; i_mk++){

		ptr_mk->lratio[i_mk] = sqrt((sp_x[i_mk]-ep_x[i_mk])*(sp_x[i_mk]-ep_x[i_mk])
    			+ (sp_y[i_mk]-ep_y[i_mk])*(sp_y[i_mk]-ep_y[i_mk]))/l0[i_mk];

		//if(ptr_u->rank == 0)
		//	printf("%d %f %f %f %f\n",i_mk,sp_x[i_mk],sp_y[i_mk],ep_x[i_mk],ep_y[i_mk]);
	}

	VecRestoreArray(ptr_mk->sp_x,&sp_x);
	VecRestoreArray(ptr_mk->sp_y,&sp_y);
	VecRestoreArray(ptr_mk->l0,	 &l0);

	VecDestroy(&dis_sp_x);
	VecDestroy(&dis_sp_y);
	VecDestroy(&dis_ep_x);
	VecDestroy(&dis_ep_y);

}

void update_sdf(Grid_S* ptr_g, Lag_marker* ptr_mk, AppCtx* ptr_u, Field_F* ptr_f)
{

	PetscScalar **dist_u,  **dist_v,  **dist_p;
	PetscScalar **panel_u, **panel_v, **panel_p;
	PetscScalar **ratio_u, **ratio_v, **ratio_p;
	PetscScalar	*sx,*sy,*ex,*ey;

	double		xd[2],yd[2]; 	//domain boundary
	double 		x1,y1,x2,y2;
	double 		be_m,be_c,be_d;  		//y=mx+c and d=sqrt(m^2+1)
	double		c,d,xb,yb,dd;
	double 		sorted_x[2], sorted_y[2];
	double 		nc[2],n[2]; 	// normal vector
	double		gridx,gridy,dx = ptr_f->dx;

	int			search_x[2], search_y[2];
	int 		panel_i,bool1, bool2,s,i,j, count = 0;
	const int	search_range = 0;

	//if (ptr_u->rank == 0){
	//	VecView(ptr_mk->sp_x,PETSC_VIEWER_STDOUT_SELF	);
	//	VecView(ptr_mk->sp_y,PETSC_VIEWER_STDOUT_SELF	);
	//}

	// get the coordinates of the markers
	VecGetArray(ptr_mk->sp_x, &sx);
	VecGetArray(ptr_mk->sp_y, &sy);
	VecGetArray(ptr_mk->ep_x, &ex);
	VecGetArray(ptr_mk->ep_y, &ey);

	// get the limits of x- y- coordinates of the local processor
	xd[0]=ptr_f->idx_range.p_raw_ghost[0][0]*dx + ptr_f->x_grid[0];
	xd[1]=ptr_f->idx_range.p_raw_ghost[0][1]*dx + ptr_f->x_grid[0];

	yd[0]=ptr_f->idx_range.p_raw_ghost[1][0]*dx + ptr_f->y_grid[0];
	yd[1]=ptr_f->idx_range.p_raw_ghost[1][1]*dx + ptr_f->y_grid[0];

	// initialize the value in panel_u, panel_v, panel_p
	// if has value 0: means haven't updated sdf yet
	// if has negative value: not in the area covered by any panels
	// if has positive value: in the area covered by panels
	VecSet(ptr_f->panel_u,0.0);
	VecSet(ptr_f->panel_v,0.0);
	VecSet(ptr_f->panel_p,0.0);

	VecSet(ptr_f->ratio_u,0.0);
	VecSet(ptr_f->ratio_v,0.0);
	VecSet(ptr_f->ratio_p,0.0);

	// Get arrays
	// In C, the indexing is "backwards" from what expects: array[k][j][i] NOT array[i][j][k]!
	DMDAVecGetArray(ptr_f->da,	ptr_f->dist_u,	&dist_u);
	DMDAVecGetArray(ptr_f->da,	ptr_f->dist_v,	&dist_v);
	DMDAVecGetArray(ptr_f->da,	ptr_f->dist_p,	&dist_p);

	DMDAVecGetArray(ptr_f->da,	ptr_f->panel_u,	&panel_u);
	DMDAVecGetArray(ptr_f->da,	ptr_f->panel_v,	&panel_v);
	DMDAVecGetArray(ptr_f->da,	ptr_f->panel_p,	&panel_p);

	DMDAVecGetArray(ptr_f->da,	ptr_f->ratio_u,	&ratio_u);
	DMDAVecGetArray(ptr_f->da,	ptr_f->ratio_v,	&ratio_v);
	DMDAVecGetArray(ptr_f->da,	ptr_f->ratio_p,	&ratio_p);

	// loop through all the boundary element
	for (panel_i=0; panel_i<ptr_mk->nglobal; panel_i++)
	{
		x1 = sx[panel_i]; x2 = ex[panel_i];
		y1 = sy[panel_i]; y2 = ey[panel_i];

		bool1 = (x1<=xd[1]) && (x1>=xd[0]) && (y1<=yd[1]) && (y1>=yd[0]);
		bool2 = (x2<=xd[1]) && (x2>=xd[0]) && (y2<=yd[1]) && (y2>=yd[0]);

		if(bool1==0 && bool2==0) // if the boundary element is not in the range owned by current processor
			continue;

		ptr_mk->list_f[count] = panel_i;
		count++;

		// setup boundary element info
		be_m = (y2-y1)/(x1-x2);
		be_c = -y1-be_m*x1;
		be_d = sqrt(be_m*be_m+1);

		//find the min and max in x, y coord
		sorted_x[0]=DMIN(x1,x2)-ptr_f->x_grid[0];
		sorted_x[1]=DMAX(x1,x2)-ptr_f->x_grid[0];
		sorted_y[0]=DMIN(y1,y2)-ptr_f->y_grid[0];
		sorted_y[1]=DMAX(y1,y2)-ptr_f->y_grid[0];

		//////////////////////////
		nc[0] = y2-y1; nc[1]=-x2+x1;  // temporary normal vector
		n[0]  = nc[0]/sqrt(nc[0]*nc[0]+nc[1]*nc[1]);
		n[1]  = nc[1]/sqrt(nc[0]*nc[0]+nc[1]*nc[1]);

		if  (x1!=x2)// if it is not a vertical panel
		{
			nc[0]=be_m; nc[1]=1;
			s=SIGN(1,n[0]*nc[0]+n[1]*nc[1]);
			ptr_mk->nx[panel_i] = s*be_m/be_d;
			ptr_mk->ny[panel_i] = s/be_d;
		}
		else
		{
			nc[0]=1; nc[1]=0;
			s=SIGN(1,n[0]*nc[0]+n[1]*nc[1]);
			ptr_mk->nx[panel_i] = SIGN(1,n[0]);
			ptr_mk->ny[panel_i] = 0;
		}

		// searching nearby grids ============== p

		// set searching range
		search_x[0] = IMAX(floor(sorted_x[0]/dx-0.5)-search_range,	ptr_f->idx_range.p_raw[0][0]);
		search_x[0] = IMIN(search_x[0],								ptr_f->idx_range.p_raw[0][1]);

		search_x[1] = IMIN(ceil	(sorted_x[1]/dx-0.5)+search_range,	ptr_f->idx_range.p_raw[0][1]);
		search_x[1] = IMAX(search_x[1],								ptr_f->idx_range.p_raw[0][0]);

		search_y[0] = IMAX(floor(sorted_y[0]/dx-0.5)-search_range,	ptr_f->idx_range.p_raw[1][0]);
		search_y[0] = IMIN(search_y[0],								ptr_f->idx_range.p_raw[1][1]);

		search_y[1] = IMIN(ceil	(sorted_y[1]/dx-0.5)+search_range,	ptr_f->idx_range.p_raw[1][1]);
		search_y[1] = IMAX(search_y[1],								ptr_f->idx_range.p_raw[1][0]);


		for (i = search_x[0]; i<=search_x[1]; i++)
			for (j = search_y[0]; j<=search_y[1]; j++)
		    {
				gridx = ptr_f->xm_grid[i];   gridy = ptr_f->ym_grid[j];

				if (x1!=x2)
				{
		            d  = (be_m*gridx + gridy +be_c)/be_d;
		            xb = gridx-d/be_d*be_m;
		            yb = gridy-d/be_d;
		            c  = (xb-x1)/(x2-x1);
				}
				else
				{
		            d = gridx-x1;
		            xb = x1;  yb = gridy;
		            c = (yb-y1)/(y2-y1);
				}


				if (c >=0 && c<=1)
				{
					//if(ptr_u->rank==0) printf("i,j,s*d dist %d %d %f %f\n",i,j,fabs(s*d),fabs(dist_p[j][i]));
					if((panel_p[j][i]==0) || (fabs(s*d)<fabs(dist_p[j][i])))
					{
						dist_p[j][i] 	= s*d;
						panel_p[j][i] 	= panel_i+1;   // The panel numbering begins from 1 here!!!!!
						ratio_p[j][i] 	= c;
					}

				}
				else if (c<0)
				{
					dd = sqrt((gridx-x1)*(gridx-x1)+(gridy-y1)*(gridy-y1));
					if((panel_p[j][i]==0) || ( dd < fabs(dist_p[j][i])))
					{
						dist_p[j][i] 	= 	SIGN(1,s*d)*dd;
						panel_p[j][i] 	= -	(panel_i+1);
						ratio_p[j][i] 	= 	0;
					}
				}
				else
				{
					dd = sqrt((gridx-x2)*(gridx-x2)+(gridy-y2)*(gridy-y2));
					if((panel_p[j][i]==0) || (dd < fabs(dist_p[j][i])))
					{
						dist_p[j][i] 	= 	SIGN(1,s*d)*dd;
						panel_p[j][i] 	= -	(panel_i+1);
						ratio_p[j][i] 	= 	1;
					}
				}

				//if(ptr_u->rank==0) printf("i,j,panel,dist,ratio %d %d %f %f %f %d\n",i,j,panel_p[j][i],dist_p[j][i],ratio_p[j][i],s);
		    }

		// searching nearby grids ============== u

		// set searching range
		search_x[0] = IMAX(floor(sorted_x[0]/dx)-search_range,		ptr_f->idx_range.u_raw[0][0]);
		search_x[0] = IMIN(search_x[0],								ptr_f->idx_range.u_raw[0][1]);

		search_x[1] = IMIN(ceil	(sorted_x[1]/dx)+search_range,		ptr_f->idx_range.u_raw[0][1]);
		search_x[1] = IMAX(search_x[1],								ptr_f->idx_range.u_raw[0][0]);

		search_y[0] = IMAX(floor(sorted_y[0]/dx-0.5)-search_range,	ptr_f->idx_range.u_raw[1][0]);
		search_y[0] = IMIN(search_y[0],								ptr_f->idx_range.u_raw[1][1]);

		search_y[1] = IMIN(ceil	(sorted_y[1]/dx-0.5)+search_range,	ptr_f->idx_range.u_raw[1][1]);
		search_y[1] = IMAX(search_y[1],								ptr_f->idx_range.u_raw[1][0]);

		for (i = search_x[0]; i<=search_x[1]; i++)
			for (j = search_y[0]; j<=search_y[1]; j++)
		    {
				gridx = ptr_f->x_grid[i];   gridy = ptr_f->ym_grid[j];

				if (x1!=x2)
				{
		            d  = (be_m*gridx + gridy +be_c)/be_d;
		            xb = gridx-d/be_d*be_m;
		            yb = gridy-d/be_d;
		            c  = (xb-x1)/(x2-x1);
				}
				else
				{
		            d = gridx-x1;
		            xb = x1;  yb = gridy;
		            c = (yb-y1)/(y2-y1);
				}


				if (c >=0 && c<=1)
				{
					if((panel_u[j][i]==0) || (fabs(s*d)<fabs(dist_u[j][i])))
					{
						dist_u[j][i] 	= s*d;
						panel_u[j][i] 	= panel_i+1;   // The panel numbering begins from 1 here!!!!!
						ratio_u[j][i] 	= c;
					}

				}
				else if (c<0)
				{
					dd = sqrt((gridx-x1)*(gridx-x1)+(gridy-y1)*(gridy-y1));
					if((panel_u[j][i]==0) || (fabs(dd)<fabs(dist_u[j][i])))
					{
						dist_u[j][i] 	= 	SIGN(1,s*d)*dd;
						panel_u[j][i] 	= -	(panel_i+1);
						ratio_u[j][i] 	= 	0;
					}
				}
				else
				{
					dd = sqrt((gridx-x2)*(gridx-x2)+(gridy-y2)*(gridy-y2));
					if((panel_u[j][i]==0) || (fabs(dd)<fabs(dist_u[j][i])))
					{
						dist_u[j][i] 	= 	SIGN(1,s*d)*dd;
						panel_u[j][i] 	= -	(panel_i+1);
						ratio_u[j][i] 	= 	1;
					}
				}

		    }

		// searching nearby grids ============== p

		// set searching range
		search_x[0] = IMAX(floor(sorted_x[0]/dx-0.5)-search_range,	ptr_f->idx_range.v_raw[0][0]);
		search_x[0] = IMIN(search_x[0],								ptr_f->idx_range.v_raw[0][1]);

		search_x[1] = IMIN(ceil	(sorted_x[1]/dx-0.5)+search_range,	ptr_f->idx_range.v_raw[0][1]);
		search_x[1] = IMAX(search_x[1],								ptr_f->idx_range.v_raw[0][0]);

		search_y[0] = IMAX(floor(sorted_y[0]/dx)-search_range,		ptr_f->idx_range.v_raw[1][0]);
		search_y[0] = IMIN(search_y[0],								ptr_f->idx_range.v_raw[1][1]);

		search_y[1] = IMIN(ceil	(sorted_y[1]/dx)+search_range,		ptr_f->idx_range.v_raw[1][1]);
		search_y[1] = IMAX(search_y[1],								ptr_f->idx_range.v_raw[1][0]);

		for (i = search_x[0]; i<=search_x[1]; i++)
			for (j = search_y[0]; j<=search_y[1]; j++)
		    {
				gridx = ptr_f->xm_grid[i];   gridy = ptr_f->y_grid[j];

				if (x1!=x2)
				{
		            d  = (be_m*gridx + gridy +be_c)/be_d;
		            xb = gridx-d/be_d*be_m;
		            yb = gridy-d/be_d;
		            c  = (xb-x1)/(x2-x1);
				}
				else
				{
		            d = gridx-x1;
		            xb = x1;  yb = gridy;
		            c = (yb-y1)/(y2-y1);
				}


				if (c >=0 && c<=1) {
					if((panel_v[j][i]==0) || (fabs(s*d)<fabs(dist_v[j][i])))
					{
						dist_v[j][i] 	= s*d;
						panel_v[j][i] 	= panel_i+1;   // The panel numbering begins from 1 here!!!!!
						ratio_v[j][i] 	= c;
					}
				}
				else if (c<0){

					dd = sqrt((gridx-x1)*(gridx-x1)+(gridy-y1)*(gridy-y1));
					if((panel_v[j][i]==0) || (abs(dd)<abs(dist_v[j][i])))
					{
						dist_v[j][i] 	= 	SIGN(1,s*d)*dd;
						panel_v[j][i] 	= -	(panel_i+1);
						ratio_v[j][i] 	= 	0;
					}
				}
				else {

					dd = sqrt((gridx-x2)*(gridx-x2)+(gridy-y2)*(gridy-y2));
					if((panel_v[j][i]==0) || (fabs(dd)<fabs(dist_v[j][i])))
					{
						dist_v[j][i] 	= 	SIGN(1,s*d)*dd;
						panel_v[j][i] 	= -	(panel_i+1);
						ratio_v[j][i] 	= 	1;
					}
				}
		    }

	}

	ptr_mk->nlocal_f = count;
	// update the sdf
	DMDAVecRestoreArray(ptr_f->da,	ptr_f->dist_u,	&dist_u);
	DMDAVecRestoreArray(ptr_f->da,	ptr_f->dist_v,	&dist_v);
	DMDAVecRestoreArray(ptr_f->da,	ptr_f->dist_p,	&dist_p);

	DMDAVecRestoreArray(ptr_f->da,	ptr_f->panel_u,	&panel_u);
	DMDAVecRestoreArray(ptr_f->da,	ptr_f->panel_v,	&panel_v);
	DMDAVecRestoreArray(ptr_f->da,	ptr_f->panel_p,	&panel_p);

	DMDAVecRestoreArray(ptr_f->da,	ptr_f->ratio_u,	&ratio_u);
	DMDAVecRestoreArray(ptr_f->da,	ptr_f->ratio_v,	&ratio_v);
	DMDAVecRestoreArray(ptr_f->da,	ptr_f->ratio_p,	&ratio_p);

	VecRestoreArray(ptr_mk->sp_x, &sx);
	VecRestoreArray(ptr_mk->sp_y, &sy);
	VecRestoreArray(ptr_mk->ep_x, &ex);
	VecRestoreArray(ptr_mk->ep_y, &ey);

	PetscViewer viewer;
	// output dist function

    PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
    PetscViewerSetType(viewer,PETSCVIEWERASCII);
    PetscViewerBinarySkipInfo(viewer); // suppress .info file
    PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);
    PetscViewerFileSetName(viewer,"f_dist_p.dat");
    VecView(ptr_f->dist_p, viewer);
    PetscViewerDestroy(&viewer);

    PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
    PetscViewerSetType(viewer,PETSCVIEWERASCII);
    PetscViewerBinarySkipInfo(viewer); // suppress .info file
    PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);
    PetscViewerFileSetName(viewer,"f_ratio_p.dat");
    VecView(ptr_f->ratio_p, viewer);
    PetscViewerDestroy(&viewer);
}


// Create forcing point list
void ib_find_forcing_pts_2d(Field_F* ptr_f, Lag_marker* ptr_mk){

	DM       			da    = ptr_f->da;
	double 				**panel_u, **panel_v, **panel_p;
	double 				**ratio_u, **ratio_v, **ratio_p;
	const double		*sx,*sy,*ex,*ey;	// read only
	forcing_point_list_t  *list;
    Vec					local_dist;
	Vec      			dist_map;          /* forcing point map (no ghost pts) */
    Vec      			dist_map_ghosted;  /* forcing point map, with correct ghost pts */
    double 				**dist,  **ratio,  **panel;

    double				**dist_map_array;
    int 	 			id[2], jd[2], kd[2];
    int 	 			pn, comp,m,b,i,j;
    double   			rt, n[2];
    double   			*x,*y;

    /* Interpolation templates */
	int       n_offset  = sizeof(offset2d)/sizeof(offset2d_t);
	dist_t   *dist_list = (dist_t *)malloc(sizeof(dist_t)*n_offset);

	forcing_pts_list_init(&(ptr_f->flist_u));
	forcing_pts_list_init(&(ptr_f->flist_v));
	forcing_pts_list_init(&(ptr_f->flist_p));

	DMGetLocalVector(ptr_f->da, &dist_map);        // local distance function map
	DMGetLocalVector(ptr_f->da, &dist_map_ghosted);
	DMGetLocalVector(ptr_f->da, &local_dist);

	VecGetArrayRead(ptr_mk->sp_x, &sx);
	VecGetArrayRead(ptr_mk->sp_y, &sy);
	VecGetArrayRead(ptr_mk->ep_x, &ex);
	VecGetArrayRead(ptr_mk->ep_y, &ey);

	DMDAVecGetArray(ptr_f->da,	ptr_f->panel_u,	&panel_u);
	DMDAVecGetArray(ptr_f->da,	ptr_f->panel_v,	&panel_v);
	DMDAVecGetArray(ptr_f->da,	ptr_f->panel_p,	&panel_p);

	DMDAVecGetArray(ptr_f->da,	ptr_f->ratio_u,	&ratio_u);
	DMDAVecGetArray(ptr_f->da,	ptr_f->ratio_v,	&ratio_v);
	DMDAVecGetArray(ptr_f->da,	ptr_f->ratio_p,	&ratio_p);

	  for ( comp=1; comp<4; comp++) {

	        switch (comp)
	        {
	            case 1: // u component

	            	DMGlobalToLocalBegin(da,ptr_f->dist_u,INSERT_VALUES,local_dist);
	            	DMGlobalToLocalEnd	(da,ptr_f->dist_u,INSERT_VALUES,local_dist);
	                ratio= ratio_u; // added by Peggy
	                panel= panel_u; // added by Peggy
	                x    = ptr_f->x_grid;
	                y    = ptr_f->ym_grid;
	                list = &(ptr_f->flist_u);
	                index_range_get_rhs('u', ptr_f, id,  jd,  kd);
	                break;

	            case 2: // v component

	            	DMGlobalToLocalBegin(da,ptr_f->dist_v,INSERT_VALUES,local_dist);
	            	DMGlobalToLocalEnd	(da,ptr_f->dist_v,INSERT_VALUES,local_dist);
	                ratio= ratio_v; // added by Peggy
	                panel= panel_v; // added by Peggy
	                x    = ptr_f->xm_grid;
	                y    = ptr_f->y_grid;
	                list = &(ptr_f->flist_v);
	                index_range_get_rhs('v', ptr_f, id,  jd,  kd);
	                break;

	            case 3: // p component

	            	DMGlobalToLocalBegin(da,ptr_f->dist_p,INSERT_VALUES,local_dist);
	            	DMGlobalToLocalEnd	(da,ptr_f->dist_p,INSERT_VALUES,local_dist);
	                ratio= ratio_p; // added by Peggy
	                panel= panel_p; // added by Peggy
	                x    = ptr_f->xm_grid;
	                y    = ptr_f->ym_grid;
	                list = &(ptr_f->flist_p);
	                index_range_get_rhs('p', ptr_f, id,  jd,  kd);
	                break;

	        } // end of switch

	        // reset dist_map and gain array access
	        VecSet(dist_map, 0);

	        DMDAVecGetArray(da, local_dist, &dist);
	        DMDAVecGetArray(da, dist_map, 	&dist_map_array);

	        // identify forcing points
	        for ( j=jd[0]; j<=jd[1]; j++)
	            for ( i=id[0]; i<=id[1]; i++)
	            {
	            	if (dist[j][i] >= 0 &&
	                    (dist[j+1][i]<0 || dist[j-1][i]<0 || dist[j][i+1]<0 || dist[j][i-1]<0))
	                {
	                    dist_map_array[j][i] = 1.0;          // mark on the map
	                    forcing_pts_list_addto(list,i,j,0);  // add point
	                }
	            }

	        DMDAVecRestoreArray(da, dist_map, 	&dist_map_array);

	        // fill the ghost values in dist_map_ghosted.
	        DMLocalToLocalBegin(da, dist_map, INSERT_VALUES, dist_map_ghosted);
	        DMLocalToLocalEnd  (da, dist_map, INSERT_VALUES, dist_map_ghosted);
	        DMDAVecGetArray    (da, dist_map_ghosted, &dist_map_array);

	        for ( m=0; m<list->total; m++) {   // all forcing points

	        	int i = list->data[m].i[0];
	            int j = list->data[m].i[1];
	            int k = list->data[m].i[2];
	            int count = 0;

	            // loop through all offset points around one forcing point
	            for (b=0; b<n_offset; b++) {

	            	int ioff = offset2d[b].i;
	                int joff = offset2d[b].j;

	                //printf("rank, m, total %d %d %d\n",ptr_u->rank,i+ioff,j+joff);
	                if ( dist[j+joff][i+ioff] > 0 &&
	                    dist_map_array[j+joff][i+ioff] == 0 )
	                {
	                    dist_list[count].i = b;
	                    dist_list[count].d =   pow(x[i]-x[i+ioff],2)
	                    + pow(y[j]-y[j+joff],2);
	                    count ++;
	                } // end if
	            }// b loop

	            if (count < 2)
	                printf("!!! (%d,%d) has ONLY %d interplation template points\n",
	                       i,j,count);

	            // sort offset points using their distances to the forcing point
	            qsort(dist_list, count, sizeof(dist_t), (compfn)comp_distance);

	            // Pick 2 closest points, storing them in data[id].tp[].
	            for ( b=0; b<2; b++) {
	                list->data[m].tp[b] = dist_list[b].i;
	            }

	            // calculate the coordiante of the surface point whose normal
	            // vector passes the forcing point.

	            pn = fabs(panel[j][i])-1; //
	            rt = ratio[j][i];

	            list->data[m].surf[0]  = rt*(ex[pn]-sx[pn])+sx[pn];
	            list->data[m].surf[1]  = rt*(ey[pn]-sy[pn])+sy[pn];
	            list->data[m].panel 	= pn;
	            list->data[m].ratio 	= rt;

	            if (rt >=0 && rt<=1){
	            	list->data[m].n[0] = ptr_mk->nx[pn];
	            	list->data[m].n[1] = ptr_mk->ny[pn];
	            	list->data[m].n[2] = 0;
	            }
	            else {
	            	n[0] = x[i]-list->data[m].surf[0];
	            	n[1] = y[j]-list->data[m].surf[1];
	            	list->data[m].n[0] = n[0]/sqrt(n[0]*n[0]+n[1]*n[1]);
	            	list->data[m].n[1] = n[1]/sqrt(n[0]*n[0]+n[1]*n[1]);
	            	list->data[m].n[2] = 0;
	            }
	        } // id loop

	        DMDAVecRestoreArray(da, dist_map_ghosted, &dist_map_array);
	        DMDAVecRestoreArray(da, local_dist, 	  &dist);
	} // end of comp loop


	DMRestoreLocalVector(da, &local_dist);
	DMRestoreLocalVector(da, &dist_map);         // local distance function map
	DMRestoreLocalVector(da, &dist_map_ghosted);

	// -------------------------------------------------------

	DMDAVecRestoreArray(da,	ptr_f->panel_u,	&panel_u);
	DMDAVecRestoreArray(da,	ptr_f->panel_v,	&panel_v);
	DMDAVecRestoreArray(da,	ptr_f->panel_p,	&panel_p);

	DMDAVecRestoreArray(da,	ptr_f->ratio_u,	&ratio_u);
	DMDAVecRestoreArray(da,	ptr_f->ratio_v,	&ratio_v);
	DMDAVecRestoreArray(da,	ptr_f->ratio_p,	&ratio_p);

	VecRestoreArrayRead(ptr_mk->sp_x, &sx);
	VecRestoreArrayRead(ptr_mk->sp_y, &sy);
	VecRestoreArrayRead(ptr_mk->ep_x, &ex);
	VecRestoreArrayRead(ptr_mk->ep_y, &ey);
}


static int forcing_pts_list_addto(forcing_point_list_t *list,
                                  int i, int j, int k)
{
    /* add a forcing point to the list */
    list->total ++;

    /* allocate more space if necessary */
    if (list->max < list->total) {
        list->max += 10;
        list->data = (forcing_point_t *)
        realloc(list->data, sizeof(forcing_point_t)*(list->max));
    }

    /* add point to the list */
    list->data[list->total - 1].i[0] = i;
    list->data[list->total - 1].i[1] = j;
    list->data[list->total - 1].i[2] = k;

    return 0;
}

void apply_forcing_2d(Field_F *f, const double sor)
{
    //int rank, verb = f->verbose_mode;
    forcing_point_list_t  *list;
    Vec u_vec, v_vec;

    double **u_a, **v_a, **u_real, **v_real;
    double b[3], c[3], a[9], *x, *y;
    int 	id,i;

    //MPI_Comm_rank(f->comm, &rank);
    //if (rank==0 && verb) {
    //    printf("\n\n########## apply forcing to u ###########\n\n");
    //}

    DMGetLocalVector	(f->da, &u_vec);
    DMGlobalToLocalBegin(f->da, f->u, INSERT_VALUES, u_vec);
    DMGlobalToLocalEnd  (f->da, f->u, INSERT_VALUES, u_vec);
    DMDAVecGetArray		(f->da, u_vec, &u_a);

    DMGetLocalVector	(f->da, &v_vec);
    DMGlobalToLocalBegin(f->da, f->v, INSERT_VALUES, v_vec);
    DMGlobalToLocalEnd  (f->da, f->v, INSERT_VALUES, v_vec);
    DMDAVecGetArray		(f->da, v_vec, &v_a);

    DMDAVecGetArray		(f->da, f->u, &u_real);
    DMDAVecGetArray		(f->da, f->v, &v_real);

    /* forcing points for u */

    x = f->x_grid;
    y = f->ym_grid;
    list = &(f->flist_u);

    for (id=0; id < list->total; id++) {

        int i = list->data[id].i[0];
        int j = list->data[id].i[1];

        double xs = list->data[id].surf[0];    /* surface point */
        double ys = list->data[id].surf[1];

        int tp1 = list->data[id].tp[0];
        int tp2 = list->data[id].tp[1];
        int i1 = i + offset2d[tp1].i;
        int j1 = j + offset2d[tp1].j;
        int i2 = i + offset2d[tp2].i;
        int j2 = j + offset2d[tp2].j;

        double x1 = x[i1]; double y1 = y[j1];  /* fluid point 1 */
        double x2 = x[i2]; double y2 = y[j2];  /* fluid point 2 */

        /* coefficient matrix A, 3x3, column-major order */
        a[0] = 1; a[3] = xs;  a[6] = ys;
        a[1] = 1; a[4] = x1;  a[7] = y1;
        a[2] = 1; a[5] = x2;  a[8] = y2;

        PetscKernel_A_gets_inverse_A_3(a,0.0);

        b[0] = 0.0;           /* b.c. on the surface */
        b[1] = u_a[j1][i1];   /* velocity at (i1,j1) */
        b[2] = u_a[j2][i2];   /* velocity at (i1,j1) */

        A_times_b_3x3(a,b,c);

        /* apply forcing: u(interp) = c0 + c1*x + c2*y */

        double u_intrp = c[0] + c[1]*x[i] + c[2]*y[j];

        u_real[j][i] = (1.0-sor)*u_real[j][i] + sor*u_intrp;
    }

    /* forcing points for v */

    x = f->xm_grid;
    y = f->y_grid;
    list = &(f->flist_v);

    for (id=0; id < list->total; id++) {

        int i = list->data[id].i[0];
        int j = list->data[id].i[1];

        double xs = list->data[id].surf[0];    /* surface point */
        double ys = list->data[id].surf[1];

        int tp1 = list->data[id].tp[0];
        int tp2 = list->data[id].tp[1];
        int i1 = i + offset2d[tp1].i;
        int j1 = j + offset2d[tp1].j;
        int i2 = i + offset2d[tp2].i;
        int j2 = j + offset2d[tp2].j;

        double x1 = x[i1]; double y1 = y[j1];  /* fluid point 1 */
        double x2 = x[i2]; double y2 = y[j2];  /* fluid point 2 */

        /* coefficient matrix A, 3x3, column-major order */
        a[0] =  1; a[3] = xs; a[6] = ys;
        a[1] =  1; a[4] = x1; a[7] = y1;
        a[2] =  1; a[5] = x2; a[8] = y2;

        PetscKernel_A_gets_inverse_A_3(a,0.0);

        b[0] = 0.0;           /* b.c. on the surface */
        b[1] = v_a[j1][i1];   /* velocity at (i1,j1) */
        b[2] = v_a[j2][i2];   /* velocity at (i1,j1) */

        A_times_b_3x3(a,b,c);

        /* apply forcing: v(interp) = c0 + c1*x + c2*y */

        double v_intrp = c[0] + c[1]*x[i] + c[2]*y[j];

        v_real[j][i] =  (1.0-sor)*v_real[j][i] + sor*v_intrp;

    }

    DMDAVecRestoreArray	(f->da, u_vec, &u_a);
    DMDAVecRestoreArray	(f->da, v_vec, &v_a);
    DMDAVecRestoreArray	(f->da, f->u, &u_real);
    DMDAVecRestoreArray	(f->da, f->v, &v_real);
    DMRestoreLocalVector(f->da, &u_vec);
    DMRestoreLocalVector(f->da, &v_vec);

	PetscViewer viewer;
	// output dist function

    PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
    PetscViewerSetType(viewer,PETSCVIEWERASCII);
    PetscViewerBinarySkipInfo(viewer); // suppress .info file
    PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);
    PetscViewerFileSetName(viewer,"f_vel_u.dat");
    VecView(f->u, viewer);
    PetscViewerDestroy(&viewer);

}


void force_calculation(Index_S* ptr_i, Field_F* ptr_f, Field_S* ptr_s, Lag_marker* ptr_mk){

	double	*par_ux, *par_uy, *par_vx, *par_vy, *pressure;
	const double  *sx, *sy, *ex, *ey;
	int		id, panel;
	int 	i,j, tp1, tp2, i1, i2, j1, j2;
	double	interp_sten[3], a[9] ,c[3], ratio; //
	Vec		uvp_vec;
	double  **uvp_array;
	DM		da = ptr_f->da;
	forcing_point_list_t *list;

	par_ux 	 = dvector(0, ptr_f->flist_u.total-1);
	par_uy 	 = dvector(0, ptr_f->flist_u.total-1);
	par_vx 	 = dvector(0, ptr_f->flist_v.total-1);
	par_vy 	 = dvector(0, ptr_f->flist_v.total-1);
	pressure = dvector(0, ptr_f->flist_p.total-1);

	// ============ part I: getting velocity gradient and pressure =============

	DMGetLocalVector	(da, &uvp_vec);

	// u
    DMGlobalToLocalBegin(da, ptr_f->u, INSERT_VALUES, uvp_vec);
    DMGlobalToLocalEnd  (da, ptr_f->u, INSERT_VALUES, uvp_vec);
    DMDAVecGetArray		(da, uvp_vec, &uvp_array);

	VecGetArrayRead(ptr_mk->sp_u, &sx);
	VecGetArrayRead(ptr_mk->ep_u, &ex);

    list = &(ptr_f->flist_u);
	for( id = 0; id<ptr_f->flist_u.total; id++){

		i 	= list->data[id].i[0];
		j 	= list->data[id].i[1];
        tp1 = list->data[id].tp[0];
        tp2 = list->data[id].tp[1];
        i1 = i + offset2d[tp1].i;
        j1 = j + offset2d[tp1].j;
        i2 = i + offset2d[tp2].i;
        j2 = j + offset2d[tp2].j;
        panel = list->data[id].panel;
        ratio = list->data[id].ratio;

		interp_sten[0] = sx[panel]+ratio*(ex[panel]-sx[panel]); /// !!!!!!!!!!!!!!!!!!!!!!!! no velocity
		interp_sten[1] = uvp_array[j1][i1];
		interp_sten[2] = uvp_array[j2][i2];

        a[0] =  1; a[3] = list->data[id].surf[0]; a[6] = list->data[id].surf[1];
        a[1] =  1; a[4] = ptr_f->x_grid[i1]; 	  a[7] = ptr_f->ym_grid[j1];
        a[2] =  1; a[5] = ptr_f->x_grid[i2]; 	  a[8] = ptr_f->ym_grid[j2];

        PetscKernel_A_gets_inverse_A_3(a,0.0);
        A_times_b_3x3(a,interp_sten,c);

        par_ux[id]=c[1];  par_uy[id]=c[2];

	}
	DMDAVecRestoreArray		(da, uvp_vec, &uvp_array);

	VecRestoreArrayRead(ptr_mk->sp_u, &sx);
	VecRestoreArrayRead(ptr_mk->ep_u, &ex);

	// v
    DMGlobalToLocalBegin(da, ptr_f->v, INSERT_VALUES, uvp_vec);
    DMGlobalToLocalEnd  (da, ptr_f->v, INSERT_VALUES, uvp_vec);
    DMDAVecGetArray		(da, uvp_vec, &uvp_array);

	VecGetArrayRead(ptr_mk->sp_v, &sy);
	VecGetArrayRead(ptr_mk->ep_v, &ey);

    list = &(ptr_f->flist_v);
	for( id = 0; id<ptr_f->flist_v.total; id++){

		i 	= list->data[id].i[0];
		j 	= list->data[id].i[1];
        tp1 = list->data[id].tp[0];
        tp2 = list->data[id].tp[1];
        i1 = i + offset2d[tp1].i;
        j1 = j + offset2d[tp1].j;
        i2 = i + offset2d[tp2].i;
        j2 = j + offset2d[tp2].j;
        panel = list->data[id].panel;
        ratio = list->data[id].ratio;

		interp_sten[0] = sy[panel]+ratio*(ey[panel]-sy[panel]); /// !!!!!!!!!!!!!!!!!!!!!!!! no velocity !
		interp_sten[1] = uvp_array[j1][i1];
		interp_sten[2] = uvp_array[j2][i2];

        a[0] =  1; a[3] = list->data[id].surf[0]; a[6] = list->data[id].surf[1];
        a[1] =  1; a[4] = ptr_f->xm_grid[i1]; 	  a[7] = ptr_f->y_grid[j1];
        a[2] =  1; a[5] = ptr_f->xm_grid[i2]; 	  a[8] = ptr_f->y_grid[j2];

        PetscKernel_A_gets_inverse_A_3(a,0.0);
        A_times_b_3x3(a,interp_sten,c);

        par_vx[id]=c[1];  par_vy[id]=c[2];

	}
	DMDAVecRestoreArray		(da, uvp_vec, &uvp_array);

	VecRestoreArrayRead(ptr_mk->sp_v, &sy);
	VecRestoreArrayRead(ptr_mk->ep_v, &ey);

	// p

    DMGlobalToLocalBegin(da, ptr_f->p, INSERT_VALUES, uvp_vec);
    DMGlobalToLocalEnd  (da, ptr_f->p, INSERT_VALUES, uvp_vec);
    DMDAVecGetArray		(da, uvp_vec, &uvp_array);

	VecGetArrayRead(ptr_mk->sp_ax, &sx);
	VecGetArrayRead(ptr_mk->ep_ax, &ex);
	VecGetArrayRead(ptr_mk->sp_ay, &sy);
	VecGetArrayRead(ptr_mk->ep_ay, &ey);

    for( id = 0; id<ptr_f->flist_p.total; id++){

		i 	= list->data[id].i[0];
		j 	= list->data[id].i[1];
        tp1 = list->data[id].tp[0];
        tp2 = list->data[id].tp[1];
        i1 = i + offset2d[tp1].i;
        j1 = j + offset2d[tp1].j;
        i2 = i + offset2d[tp2].i;
        j2 = j + offset2d[tp2].j;
        panel = list->data[id].panel;
        ratio = list->data[id].ratio;

		interp_sten[0] = -list->data[id].n[0]*(sx[panel]+ratio*(ex[panel]-sx[panel]))
						 -list->data[id].n[1]*(sy[panel]+ratio*(ey[panel]-sy[panel])); /// !!!!!!!!!!!!!!!!!!!!!!!! no acceleration
		interp_sten[1] = uvp_array[j1][i1];
		interp_sten[2] = uvp_array[j2][i2];

        a[0] =  0; a[3] = list->data[id].n[0]; 	  a[6] = list->data[id].n[1];
        a[1] =  1; a[4] = ptr_f->xm_grid[i1]; 	  a[7] = ptr_f->ym_grid[j1];
        a[2] =  1; a[5] = ptr_f->xm_grid[i2]; 	  a[8] = ptr_f->ym_grid[j2];

        PetscKernel_A_gets_inverse_A_3(a,0.0);
    	A_times_b_3x3(a,interp_sten,c);

    	pressure[id]=c[0]+ptr_f->xm_grid[i]*c[1]+ptr_f->ym_grid[j]*c[2];
    }
    DMDAVecRestoreArray		(da, uvp_vec, &uvp_array);
	DMRestoreLocalVector	(da, &uvp_vec);

	VecRestoreArrayRead(ptr_mk->sp_ax, &sx);
	VecRestoreArrayRead(ptr_mk->ep_ax, &ex);
	VecRestoreArrayRead(ptr_mk->sp_ay, &sy);
	VecRestoreArrayRead(ptr_mk->ep_ay, &ey);

	// ================== part II: calculate traction on traction cell ===========

	int			ii;
	int			id_t, comp;
	short int	*C2c_neumann;
	double 		ref_pressure = 0, *traction_cell, *n;
	Vec			trac_n_vec,trac_temp_vec,trac_check_vec,trac_checktot_vec;
	PetscScalar cc;

	VecView(ptr_s->traction,PETSC_VIEWER_STDOUT_WORLD);

	VecDuplicate(ptr_s->traction, &trac_n_vec);
	VecDuplicate(ptr_s->traction, &trac_temp_vec);
	VecDuplicate(ptr_s->traction, &trac_check_vec);
	VecDuplicate(ptr_s->traction, &trac_checktot_vec);

	VecSet(trac_n_vec		,0.0);
	VecSet(trac_temp_vec	,0.0);
	VecSet(trac_check_vec	,0.0);
	VecSet(trac_checktot_vec,0.0);

	// du/dx du/dy part for tx
	for( ii = 0; ii<ptr_f->flist_u.total; ii++){

		id 		= ptr_f->flist_u.data[ii].panel;
		id_t 	= ptr_mk->tx[id];

		//printf("id,id_t %d %d \n",id,id_t);
		cc = 1;
		VecSetValue(trac_n_vec,		id_t,cc, ADD_VALUES);
		VecSetValue(trac_check_vec,	id_t,cc, INSERT_VALUES );

		cc = 2*par_ux[ii]*ptr_mk->nx[id] + par_uy[ii]*ptr_mk->ny[id];
		VecSetValue(trac_temp_vec,	id_t,cc, ADD_VALUES);
	}

	// du/dy part for ty
	for( ii = 0; ii<ptr_f->flist_u.total; ii++){

		id 		= ptr_f->flist_u.data[ii].panel;
		id_t 	= ptr_mk->ty[id];

		cc = 1;
		VecSetValue(trac_n_vec,		id_t,cc, ADD_VALUES);
		VecSetValue(trac_check_vec,	id_t,cc, INSERT_VALUES );

		cc = par_uy[ii]*ptr_mk->nx[id];
		VecSetValue(trac_temp_vec,	id_t,cc, ADD_VALUES);
	}

	VecAssemblyBegin(trac_n_vec);
	VecAssemblyEnd	(trac_n_vec);

	VecAssemblyBegin(trac_check_vec);
	VecAssemblyEnd	(trac_check_vec);

	VecAssemblyBegin(trac_temp_vec);
	VecAssemblyEnd	(trac_temp_vec);

	VecPointwiseDivide(trac_temp_vec,trac_temp_vec,trac_n_vec);
	VecCopy(trac_temp_vec, ptr_s->traction);

	cc = 1;
	VecAXPY(trac_checktot_vec,cc,trac_check_vec);

	VecSet(trac_n_vec		,0.0);
	VecSet(trac_temp_vec	,0.0);
	VecSet(trac_check_vec	,0.0);

	// dv/dx part for tx
	for( ii = 0; ii<ptr_f->flist_v.total; ii++){

		id 		= ptr_f->flist_v.data[ii].panel;
		id_t 	= ptr_mk->tx[id];

		cc = 1;
		VecSetValue(trac_n_vec,		id_t,cc, ADD_VALUES);
		VecSetValue(trac_check_vec,	id_t,cc, INSERT_VALUES );

		cc = par_vx[ii]*ptr_mk->ny[id];
		VecSetValue(trac_temp_vec,	id_t,cc, ADD_VALUES);
	}

	// dv/dx, dv/dy part for ty
	for( ii = 0; ii<ptr_f->flist_v.total; ii++){

		id 		= ptr_f->flist_v.data[ii].panel;
		id_t 	= ptr_mk->ty[id];

		cc = 1;
		VecSetValue(trac_n_vec,		id_t,cc, ADD_VALUES);
		VecSetValue(trac_check_vec,	id_t,cc, INSERT_VALUES );

		cc = 2*par_vy[ii]*ptr_mk->ny[id] + par_vx[ii]*ptr_mk->nx[id];
		VecSetValue(trac_temp_vec,	id_t,cc, ADD_VALUES);
	}

	VecAssemblyBegin(trac_n_vec);
	VecAssemblyEnd	(trac_n_vec);

	VecAssemblyBegin(trac_check_vec);
	VecAssemblyEnd	(trac_check_vec);

	VecAssemblyBegin(trac_temp_vec);
	VecAssemblyEnd	(trac_temp_vec);

	VecPointwiseDivide(trac_temp_vec,trac_temp_vec,trac_n_vec);
	cc = 1.0;
	VecAXPY(trac_temp_vec,cc,ptr_s->traction);
	VecAXPY(trac_checktot_vec,cc,trac_check_vec);

	cc = 1/ptr_f->reynolds_number;
	VecScale(ptr_s->traction, cc);

	// pressure part
	VecSet(trac_n_vec		,0.0);
	VecSet(trac_temp_vec	,0.0);
	VecSet(trac_check_vec	,0.0);

	// for tx
	for( ii = 0; ii<ptr_f->flist_p.total; ii++){

		id 		= ptr_f->flist_p.data[ii].panel;
		id_t 	= ptr_mk->tx[id];

		cc = 1;
		VecSetValue(trac_n_vec,		id_t,cc, ADD_VALUES);
		VecSetValue(trac_check_vec,	id_t,cc, INSERT_VALUES );

		cc = -(pressure[id])*ptr_mk->nx[id];
		VecSetValue(trac_temp_vec,	id_t,cc, ADD_VALUES);
	}

	for( ii = 0; ii<ptr_f->flist_p.total; ii++){

		id 		= ptr_f->flist_p.data[ii].panel;
		id_t 	= ptr_mk->ty[id];

		cc = 1;
		VecSetValue(trac_n_vec,		id_t,cc, ADD_VALUES);
		VecSetValue(trac_check_vec,	id_t,cc, INSERT_VALUES );

		cc = -(pressure[id])*ptr_mk->ny[id];
		VecSetValue(trac_temp_vec,	id_t,cc, ADD_VALUES);
	}

	VecAssemblyBegin(trac_n_vec);
	VecAssemblyEnd	(trac_n_vec);

	VecAssemblyBegin(trac_check_vec);
	VecAssemblyEnd	(trac_check_vec);

	VecAssemblyBegin(trac_temp_vec);
	VecAssemblyEnd	(trac_temp_vec);

	VecPointwiseDivide(trac_temp_vec,trac_temp_vec,trac_n_vec);
	cc = 1.0;
	VecAXPY(trac_temp_vec,cc,ptr_s->traction);
	VecAXPY(trac_checktot_vec,cc,trac_check_vec);

	// treatment of boundary cell with no velocity gradients or pressure


}

static int forcing_pts_list_init(forcing_point_list_t *list)
{
    list->total = 0; /* no element */
    list->max   = 1; /* this is arbitrarily chosen */
    list->data  = (forcing_point_t *)malloc(sizeof(forcing_point_t)*(list->max));
    return 0;
}

void A_times_b_3x3(double *a, double *b, double *c)
{
    /* compute c = A*b, where A is 3x3 (column-major order) and b is 3x1 */

    c[0] = a[0]*b[0] + a[3]*b[1] + a[6]*b[2];
    c[1] = a[1]*b[0] + a[4]*b[1] + a[7]*b[2];
    c[2] = a[2]*b[0] + a[5]*b[1] + a[8]*b[2];

}

void displacement_interpolation(Lag_marker* ptr_mk, Index_S* ptr_i, Grid_S* ptr_g ){ // 11/13/2015 by Peggy

	double		dx = ptr_g->dx, dx2 = dx/2;
	double		c1, c2, xx, yy;
	int			nx = ptr_g->Nx;
	int 		ii,jj,i_mk,comp,i,j;
	int			*G2l;
	double		*x, *y;
	struct interp2d *interp;
	PetscScalar	*sx_array,*sy_array,*ex_array,*ey_array;

	ptr_mk->interp_s_x = (struct interp2d*)calloc(ptr_mk->nlocal_s, sizeof(struct interp2d));
	ptr_mk->interp_s_y = (struct interp2d*)calloc(ptr_mk->nlocal_s, sizeof(struct interp2d));
	ptr_mk->interp_e_x = (struct interp2d*)calloc(ptr_mk->nlocal_s, sizeof(struct interp2d));
	ptr_mk->interp_e_y = (struct interp2d*)calloc(ptr_mk->nlocal_s, sizeof(struct interp2d));

    VecGetArray(ptr_mk->sp_x0,&sx_array);
    VecGetArray(ptr_mk->sp_y0,&sy_array);
    VecGetArray(ptr_mk->ep_x0,&ex_array);
    VecGetArray(ptr_mk->ep_y0,&ey_array);

    G2l = (int *) malloc(sizeof(int)*ptr_g->N);
    for(i=0; i<ptr_g->N; i++) G2l[i]= -1;
    for(i=0; i<ptr_i->xix_N +ptr_i->xix_ghoN; i++) G2l[ptr_i->xix.l2G[i]]=i;
    // ASSUMING uniform dx and dx=dy

	// for xi_x interpolation
    for (comp =0; comp<2; comp++){

    	switch(comp)
    	{
    	case 0: // starting point
    		x = sx_array;
    		y = sy_array;
    		interp = ptr_mk->interp_s_x;
    		break;
    	case 1: // ending point
    		x = ex_array;
    		y = ey_array;
    		interp = ptr_mk->interp_e_x;
    		break;
    	}

    	for (i_mk = 0; i_mk<ptr_mk->nlocal_s; i_mk++){

    		xx = x[i_mk]-ptr_g->x_grid[0];
    		yy = y[i_mk]-ptr_g->y_grid[0];

    		ii = floor(xx/dx);
    		jj = floor(yy/dx-0.5);

    		c1 = (xx	-ii*dx)/dx;	 	// local x coordinate defined by (x-x1)/(x2-x1)
    		c2 = (yy-dx2-jj*dx)/dx; 	// local y coordinate defined by (y-y1)/(y2-y1)

    		// if marker are outside of domain
    		// (to  be done)

    		interp[i_mk].sten[0] = G2l[ii 	+ nx*jj		];
    		interp[i_mk].sten[1] = G2l[ii+1 + nx*jj		];
    		interp[i_mk].sten[2] = G2l[ii 	+ nx*(jj+1) ];
    		interp[i_mk].sten[3] = G2l[ii+1 + nx*(jj+1) ];

    		interp[i_mk].coeff[0] = (1-c1)*(1-c2);
    		interp[i_mk].coeff[1] = c1	   *(1-c2);
    		interp[i_mk].coeff[2] = (1-c1)*c2	  ;
    		interp[i_mk].coeff[3] = c1	   *c2	  ;

    		interp[i_mk].n		  = 4;

    		for(i=3; i>=0; i--){

    			// get rid of the stencil that are not involved
    			if(interp[i_mk].sten[i] == -1 && interp[i_mk].coeff[i] < 1e-15){

    				interp[i_mk].n --;
    				for(j=i; j<interp[i_mk].n; j++) {
    					interp[i_mk].sten[j]  = interp[i_mk].sten[j+1];
    					interp[i_mk].coeff[j] = interp[i_mk].coeff[j+1];
    				}
    			}
    			else if (interp[i_mk].sten[i] == -1)
    				printf("wrong bilinear interpolation of marker at %f %f\n",x[i_mk],y[i_mk]);

    		}
    	} // end of mk loop
    } // end of comp loop

    G2l = (int *) malloc(sizeof(int)*ptr_g->N);
    for(i=0; i<ptr_g->N; i++) G2l[i]= -1;
    for(i=0; i<ptr_i->xiy_N +ptr_i->xiy_ghoN; i++) G2l[ptr_i->xiy.l2G[i]]=i;

	// for xi_y interpolation
    for (comp =0; comp<2; comp++){

    	switch(comp)
    	{
    	case 0: // starting point
    		x = sx_array;
    		y = sy_array;
    		interp = ptr_mk->interp_s_y;
    		break;
    	case 1: // ending point
    		x = ex_array;
    		y = ey_array;
    		interp = ptr_mk->interp_e_y;
    		break;
    	}

    	for (i_mk = 0; i_mk<ptr_mk->nlocal_s; i_mk++){

    		xx = x[i_mk]-ptr_g->x_grid[0];
    		yy = y[i_mk]-ptr_g->y_grid[0];

    		ii = floor(xx/dx-0.5);
    		jj = floor(yy/dx);

    		// if marker are outside of domain
    		// (to  be done)

    		c1 = (xx-dx2-ii*dx)/dx;		// local x coordinate defined by (x-x1)/(x2-x1)
    		c2 = (yy-	 jj*dx)/dx;  	// local y coordinate defined by (y-y1)/(y2-y1)

    		interp[i_mk].sten[0] = G2l[ii 	+ nx*jj		];
    		interp[i_mk].sten[1] = G2l[ii+1 + nx*jj		];
    		interp[i_mk].sten[2] = G2l[ii 	+ nx*(jj+1) ];
    		interp[i_mk].sten[3] = G2l[ii+1 + nx*(jj+1) ];

    		interp[i_mk].coeff[0] = (1-c1)*(1-c2);
    		interp[i_mk].coeff[1] = c1	   *(1-c2);
    		interp[i_mk].coeff[2] = (1-c1)*c2	  ;
    		interp[i_mk].coeff[3] = c1	   *c2	  ;

    		interp[i_mk].n		  = 4;

    		for(i=3; i>=0; i--){

    			// get rid of the stencil that are not involved
    			if(interp[i_mk].sten[i] == -1 && interp[i_mk].coeff[i] < 1e-15){

    				interp[i_mk].n --;
    				for(j=i; j<interp[i_mk].n; j++) {
    					interp[i_mk].sten[j]  = interp[i_mk].sten[j+1];
    					interp[i_mk].coeff[j] = interp[i_mk].coeff[j+1];
    				}
    			}
    			else if (interp[i_mk].sten[i] == -1)
    				printf("wrong bilinear interpolation of marker at %f %f\n",x[i_mk],y[i_mk]);

    		}
    	} // end of mk loop
    } // end of comp loop

    VecRestoreArray(ptr_mk->sp_x0,&sx_array);
    VecRestoreArray(ptr_mk->sp_y0,&sy_array);
    VecRestoreArray(ptr_mk->ep_x0,&ex_array);
    VecRestoreArray(ptr_mk->ep_y0,&ey_array);
}


