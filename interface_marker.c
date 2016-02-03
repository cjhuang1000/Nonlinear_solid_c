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

static int forcing_pts_list_addto(forcing_point_list_t *list,
                                  int i, int j, int k);

void marker_setup (Grid_S* ptr_g, Solid* ptr_s, Index_S* ptr_i,Lag_marker* ptr_mk, AppCtx* ptr_u)
{
	const int  	MAX_N = ptr_g->N;  // maximum number of markers
	double      sdf[9],x[9],y[9];
	double		n[2],t[2];
	int			index = 0, temp[4][4] = {{0,1,4,3},{1,2,5,4},{3,4,7,6},{4,5,8,7}};
	Vec			tx_vec, ty_vec, l0_vec;
	IS			is;
	VecScatter	scatter;

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
    VecDuplicate(ptr_mk->sp_x0, &l0_vec);

    ptr_mk->nx 	= dvector(0, ptr_mk->nglobal-1);
    ptr_mk->ny 	= dvector(0, ptr_mk->nglobal-1);

    // input the values in the Vec's
    VecGetArray(ptr_mk->sp_x0,&sx_array);
    VecGetArray(ptr_mk->sp_y0,&sy_array);
    VecGetArray(ptr_mk->ep_x0,&ex_array);
    VecGetArray(ptr_mk->ep_y0,&ey_array);
    VecGetArray(tx_vec		 ,&tx_array);
    VecGetArray(ty_vec   	 ,&ty_array);
    VecGetArray(l0_vec  	 ,&l0);

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
    }

    VecRestoreArray(ptr_mk->sp_x0,&sx_array);
    VecRestoreArray(ptr_mk->sp_y0,&sy_array);
    VecRestoreArray(ptr_mk->ep_x0,&ex_array);
    VecRestoreArray(ptr_mk->ep_y0,&ey_array);
    VecRestoreArray(tx_vec		 ,&tx_array);
    VecRestoreArray(ty_vec		 ,&ty_array);
    VecRestoreArray(l0_vec		 ,&l0);

    VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nglobal,&ptr_mk->tx);
    VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nglobal,&ptr_mk->ty);
    VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nglobal,&ptr_mk->l0);

	ISCreateStride	(PETSC_COMM_SELF,ptr_mk->nglobal,0,1,&is);
	VecScatterCreate(ptr_mk->sp_x0,is,ptr_mk->tx,is,&scatter);

	VecScatterBegin	(scatter,tx_vec ,ptr_mk->tx,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd	(scatter,tx_vec ,ptr_mk->tx,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterBegin	(scatter,ty_vec ,ptr_mk->ty,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd	(scatter,ty_vec ,ptr_mk->ty,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterBegin	(scatter,l0_vec ,ptr_mk->l0,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd	(scatter,l0_vec ,ptr_mk->l0,INSERT_VALUES,SCATTER_FORWARD);

    printf("Marker creation completed. \n");
    printf("The processor %d has %d markers. \n",ptr_u->rank,ptr_mk->nlocal_s);

    ISDestroy (&is);
    VecScatterDestroy(&scatter);
    VecDestroy(&tx_vec);
    VecDestroy(&ty_vec);
    VecDestroy(&l0_vec);
}


// Given the distributed array of the location of the panels,
// create the list of the indices of panel that are i

void list_setup(Lag_marker* ptr_mk, Field_F* ptr_f){

    IS 			is, is2;
    Vec			sx0, sy0, ex0, ey0;
    Vec			tx, ty;
    int			*panel_int, *panel_gho;
    VecScatter  all2loc;
    PetscScalar *sx,*sy,*ex,*ey;
    int 		panel, i0[2], j0[2], k0[2], i;
    double		x_int[2], y_int[2], x_gho[2], y_gho[2], mp[2];
    int 		bool1, bool2, bool3, bool4;

    // get the limits of x- y- coordinates of the local processor
    index_range_get_raw		 ('p', ptr_f, i0,  j0,  k0);
	x_int[0] = i0[0]*ptr_f->dx + ptr_f->x_grid[0];
	x_int[1] = (i0[1] + 1)*ptr_f->dx + ptr_f->x_grid[0];
	y_int[0] = j0[0]*ptr_f->dx + ptr_f->y_grid[0];
	y_int[1] = (j0[1] + 1)*ptr_f->dx + ptr_f->y_grid[0];

    index_range_get_raw_ghost('p', ptr_f, i0,  j0,  k0);
	x_gho[0] = i0[0]*ptr_f->dx + ptr_f->x_grid[0];
	x_gho[1] = (i0[1] + 1)*ptr_f->dx + ptr_f->x_grid[0];
	y_gho[0] = j0[0]*ptr_f->dx + ptr_f->y_grid[0];
	y_gho[1] = (j0[1] + 1)*ptr_f->dx + ptr_f->y_grid[0];

	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nglobal,&sx0);
	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nglobal,&sy0);
	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nglobal,&ex0);
	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nglobal,&ey0);

	ISCreateStride	(PETSC_COMM_SELF,ptr_mk->nglobal,0,1,&is);
	VecScatterCreate(ptr_mk->sp_x0,is,sx0,is,&all2loc);

	// Scatter all the panel location information to every processor
	VecScatterBegin	(all2loc,ptr_mk->sp_x0,sx0,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd	(all2loc,ptr_mk->sp_x0,sx0,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterBegin	(all2loc,ptr_mk->sp_y0,sy0,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd	(all2loc,ptr_mk->sp_y0,sy0,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterBegin	(all2loc,ptr_mk->ep_x0,ex0,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd	(all2loc,ptr_mk->ep_x0,ex0,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterBegin	(all2loc,ptr_mk->ep_y0,ey0,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd	(all2loc,ptr_mk->ep_y0,ey0,INSERT_VALUES,SCATTER_FORWARD);

	VecGetArray		(sx0,&sx);
    VecGetArray		(sy0,&sy);
	VecGetArray		(ex0,&ex);
    VecGetArray		(ey0,&ey);

    panel_int 	= ivector(0, ptr_mk->nglobal -1);
    panel_gho 	= ivector(0, ptr_mk->nglobal -1);
    ptr_mk->nlocal_f =0;
    ptr_mk->nlocal_gho_f =0;

    for(panel=0; panel<ptr_mk->nglobal; panel++)
    {
    	mp[0] = (sx[panel] + ex[panel])/2;
    	mp[1] = (sy[panel] + ey[panel])/2;

		bool1 = (sx[panel]<= x_gho[1]) && (sx[panel]> x_gho[0])
				&& (sy[panel]<= y_gho[1]) && (sy[panel]> y_gho[0]);
		bool2 = (ex[panel]<=x_gho[1]) && (ex[panel]> x_gho[0])
				&& (ey[panel]<=y_gho[1]) && (ey[panel]> y_gho[0]);

    	if ((mp[0]<= x_int[1]) && (mp[0]> x_int[0]) && (mp[1]<= y_int[1]) && (mp[1]> y_int[0])) // interior
    	{

    		panel_int[ptr_mk->nlocal_f] = panel;
    		ptr_mk->nlocal_f ++;

    	}
    	else if (bool1==1 || bool2==1) // ghost points
    	{
    		panel_gho[ptr_mk->nlocal_gho_f] = panel;
    		ptr_mk->nlocal_gho_f ++;
    	}
    }

	VecRestoreArray(sx0,&sx);
    VecRestoreArray(sy0,&sy);
	VecRestoreArray(ex0,&ex);
    VecRestoreArray(ey0,&ey);
    ISDestroy  	(&is);

    // prepare the list_f
    ptr_mk->list_f 	= ivector(0, ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f -1);

    for (i = 0; i < ptr_mk->nlocal_f; i++)
    	ptr_mk->list_f[i] = panel_int[i];

    for (i = 0; i < ptr_mk->nlocal_gho_f; i++)
    	ptr_mk->list_f[ptr_mk->nlocal_f + i] = panel_gho[i];

    // create VecScatter
    VecCreateSeq   (PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,&ptr_mk->sp_x);
    VecCreateSeq   (PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,&ptr_mk->sp_y);
    VecCreateSeq   (PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,&ptr_mk->ep_x);
    VecCreateSeq   (PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,&ptr_mk->ep_y);

    VecCreateSeq   (PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,&ptr_mk->sp_u);
    VecCreateSeq   (PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,&ptr_mk->sp_v);
    VecCreateSeq   (PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,&ptr_mk->ep_u);
    VecCreateSeq   (PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,&ptr_mk->ep_v);

    VecCreateSeq   (PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,&ptr_mk->sp_ax);
    VecCreateSeq   (PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,&ptr_mk->sp_ay);
    VecCreateSeq   (PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,&ptr_mk->ep_ax);
    VecCreateSeq   (PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,&ptr_mk->ep_ay);

    ISCreateGeneral(PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,ptr_mk->list_f,PETSC_COPY_VALUES,&is2);
	ISCreateStride (PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,0,1,&is);
	VecScatterCreate(ptr_mk->sp_x0,is2,ptr_mk->sp_x,is,&ptr_mk->glo2loc);

    printf("Marker list setup completed. \n");
    printf("There are %d interior marker and %d ghosted markers\n",ptr_mk->nlocal_f,ptr_mk->nlocal_gho_f);

    VecDestroy	(&sx0);
    VecDestroy	(&sy0);
    VecDestroy	(&ex0);
    VecDestroy	(&ey0);
    ISDestroy  	(&is);
    ISDestroy  	(&is2);
}

// update the location of the current marker and update the list_f
void build_current_marker(Lag_marker* ptr_mk, Field_S* ptr_s, Field_F* ptr_f, Index_S* ptr_i, AppCtx* ptr_u)
{

	PetscReal 	*xix,*xiy,*xi; // local displacement array
	PetscReal	*sp_x,*sp_y,*ep_x,*ep_y,*l0;
	PetscReal   *mkp;
	PetscScalar	c = 1.0;
	Vec		   	dis_sp_x,dis_sp_y,dis_ep_x,dis_ep_y; // displacement and new marker locations
	Vec 		sp_x_loc, sp_y_loc, ep_x_loc, ep_y_loc;  // markers new location which was previously in the domain
	Vec			xix_local, xiy_local;
	IS			is1, is2;
	MPI_Status 	status;

	int			i_mk, comp, i, elem;
	int			stencil, *mk_int, mkn_int = 0;
	int			*sl_e, *sl_w, *sl_n, *sl_s;  // send out as the ghost markers of the nearby processors
	int			sn_e = 0, sn_w = 0, sn_n = 0, sn_s = 0;
	int			*rl_e, *rl_w, *rl_n, *rl_s;  // receive the ghost markers from the nearby processors
	int			rn_e , rn_w , rn_n, rn_s;
	struct interp2d *interp;
	double		mp[2];
	int			bool1, bool2, bool3;
	double      *sbd = ptr_f->idx_range.shared_bound;
	double		*ibd = ptr_f->idx_range.int_bound;
	int			*spr = ptr_f->idx_range.shared_pro;
	int 		px = ptr_f->proc_coords[0], py = ptr_f->proc_coords[1];
	int			Nx = ptr_f->Nx, Ny = ptr_f->Ny;

	VecCreateSeq(PETSC_COMM_SELF,ptr_i->xix_N +ptr_i->xix_ghoN,&xix_local);
    VecCreateSeq(PETSC_COMM_SELF,ptr_i->xiy_N +ptr_i->xiy_ghoN,&xiy_local);

	// Create and copy original locations of markers
	VecDuplicate(ptr_mk->sp_x0,&dis_sp_x);
	VecDuplicate(ptr_mk->sp_x0,&dis_sp_y);
	VecDuplicate(ptr_mk->sp_x0,&dis_ep_x);
	VecDuplicate(ptr_mk->sp_x0,&dis_ep_y);

	// Calculate the current location of the
	VecScatterBegin(ptr_u->scatter_x,ptr_s->xi,xix_local,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd  (ptr_u->scatter_x,ptr_s->xi,xix_local,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterBegin(ptr_u->scatter_y,ptr_s->xi,xiy_local,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd  (ptr_u->scatter_y,ptr_s->xi,xiy_local,INSERT_VALUES,SCATTER_FORWARD);

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

	VecAXPY(dis_sp_x,c,ptr_mk->sp_x0);
	VecAXPY(dis_sp_y,c,ptr_mk->sp_y0);
	VecAXPY(dis_ep_x,c,ptr_mk->ep_x0);
	VecAXPY(dis_ep_y,c,ptr_mk->ep_y0);

	// Scatter the distributed updated marker location to the local vec

	VecCreateSeq(PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,&sp_x_loc);
	VecCreateSeq(PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,&sp_y_loc);
	VecCreateSeq(PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,&ep_x_loc);
	VecCreateSeq(PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,&ep_y_loc);

	VecScatterBegin(ptr_mk->glo2loc ,dis_sp_x,sp_x_loc,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd  (ptr_mk->glo2loc ,dis_sp_x,sp_x_loc,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterBegin(ptr_mk->glo2loc ,dis_sp_y,sp_y_loc,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd  (ptr_mk->glo2loc ,dis_sp_y,sp_y_loc,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterBegin(ptr_mk->glo2loc ,dis_ep_x,ep_x_loc,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd  (ptr_mk->glo2loc ,dis_ep_x,ep_x_loc,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterBegin(ptr_mk->glo2loc ,dis_ep_y,ep_y_loc,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd  (ptr_mk->glo2loc ,dis_ep_y,ep_y_loc,INSERT_VALUES,SCATTER_FORWARD);

	// Find out the updated interior markers and shared markers
	VecGetArray(sp_x_loc,&sp_x);
	VecGetArray(sp_y_loc,&sp_y);
	VecGetArray(ep_x_loc,&ep_x);
	VecGetArray(ep_y_loc,&ep_y);

	mk_int	= ivector(0, ptr_mk->nglobal-1);
	sl_e 	= ivector(0, ptr_mk->nglobal-1);
	sl_w 	= ivector(0, ptr_mk->nglobal-1);
	sl_n 	= ivector(0, ptr_mk->nglobal-1);
	sl_s 	= ivector(0, ptr_mk->nglobal-1);

	for(i=0; i < ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f; i++)
	{
	   	mp[0] = (sp_x[i] + ep_x[i])/2;
	    mp[1] = (sp_y[i] + ep_y[i])/2;

	    // interior markers
	    if (mp[0]>ibd[0] && mp[0]<=ibd[1] && mp[1]>ibd[2] && mp[1]<=ibd[3]){

	    	mk_int[mkn_int] = ptr_mk->list_f[i];
	    	mkn_int ++;

	    	// ghosted marker of left processor
	    	if ((sp_x[i]<= sbd[0]) || (ep_x[i] <= sbd[0])){
		    	sl_w[sn_w] = ptr_mk->list_f[i];
		    	sn_w ++;	}

	    	// right processor
	    	if ((sp_x[i] >= sbd[1]) || (ep_x[i] >= sbd[1])){
		    	sl_e[sn_e] = ptr_mk->list_f[i];
		    	sn_e ++;	}

	    	// lower processor
	    	if ((sp_y[i] <= sbd[2]) || (ep_y[i] <= sbd[2])){
		    	sl_s[sn_s] = ptr_mk->list_f[i];
		    	sn_s ++;	}

	    	// upper processor
	    	if ((sp_y[i] >= sbd[3]) || (ep_y[i] >= sbd[3])){
		    	sl_n[sn_n] = ptr_mk->list_f[i];
		    	sn_n ++;	}
	    }

	}
	VecRestoreArray(sp_x_loc,&sp_x);
	VecRestoreArray(sp_y_loc,&sp_y);
	VecRestoreArray(ep_x_loc,&ep_x);
	VecRestoreArray(ep_y_loc,&ep_y);

	printf("interior %d, shared %d %d %d %d\n", mkn_int, sn_w,sn_e,sn_s,sn_n);

	rn_e = 0;
	rn_w = 0;
	rn_s = 0;
	rn_n = 0;

	//printf("rank %d, %d %d %d %d\n", rank, spr[0],spr[1],spr[2],spr[3]);

	// sending and receiving
	// to left
	if (spr[0] != -1)  MPI_Send(&sn_w, 1, MPI_INT, spr[0],101, MPI_COMM_WORLD);
	if (spr[1] != -1)
	{
		MPI_Recv(&rn_e, 1, MPI_INT, spr[1],101, MPI_COMM_WORLD, &status);
		rl_e 	= ivector(0, rn_e-1);
	}
	if (spr[0] != -1) MPI_Send(sl_w, sn_w, MPI_INT,spr[0],101, MPI_COMM_WORLD);
	if (spr[1] != -1) MPI_Recv(rl_e, rn_e, MPI_INT, spr[1],101, MPI_COMM_WORLD, &status);

	// to right
	if (spr[1] != -1)  MPI_Send(&sn_e, 1, MPI_INT, spr[1],101, MPI_COMM_WORLD);
	if (spr[0] != -1)
	{
		MPI_Recv(&rn_w, 1, MPI_INT, spr[0],101, MPI_COMM_WORLD, &status);
		rl_w 	= ivector(0, rn_w-1);
	}
	if (spr[1] != -1)  MPI_Send(sl_e, sn_e, MPI_INT, spr[1],101, MPI_COMM_WORLD);
	if (spr[0] != -1)  MPI_Recv(rl_w, rn_w, MPI_INT, spr[0],101, MPI_COMM_WORLD, &status);

	// to lower
	if (spr[2] != -1)  MPI_Send(&sn_s, 1, MPI_INT, spr[2],101, MPI_COMM_WORLD);
	if (spr[3] != -1)
	{
		MPI_Recv(&rn_n, 1, MPI_INT, spr[3],101, MPI_COMM_WORLD, &status);
		rl_n 	= ivector(0, rn_n-1);
	}
	if (spr[2] != -1) MPI_Send(sl_s, sn_s, MPI_INT, spr[2],101, MPI_COMM_WORLD);
	if (spr[3] != -1) MPI_Recv(rl_n, rn_n, MPI_INT, spr[3],101, MPI_COMM_WORLD, &status);

	// to upper
	if (spr[3] != -1)  MPI_Send(&sn_n, 1, MPI_INT, spr[3],101, MPI_COMM_WORLD);
	if (spr[2] != -1)
	{
		MPI_Recv(&rn_s, 1, MPI_INT, spr[2],101, MPI_COMM_WORLD, &status);
		rl_s 	= ivector(0, rn_s-1);
	}
	if (spr[3] != -1)  MPI_Send(sl_n, sn_n, MPI_INT, spr[3],101, MPI_COMM_WORLD);
	if (spr[2] != -1)  MPI_Recv(rl_s, rn_s, MPI_INT, spr[2],101, MPI_COMM_WORLD, &status);


	//printf("rank %d marker %d %d %d %d\n",rank,rn_e,rn_w,rn_s,rn_n);
	// Updating list and new VecScatter

	free_ivector(ptr_mk->list_f, 0, ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f -1);
	ptr_mk->nlocal_f 	 = mkn_int;
	ptr_mk->nlocal_gho_f = rn_w + rn_e + rn_s + rn_n;
	ptr_mk->list_f 		 = ivector(0, ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f -1);

	i = 0; 		  for (i_mk = 0; i_mk < mkn_int; i_mk++ )  ptr_mk->list_f[i_mk+i] =  mk_int[i_mk];
	i += mkn_int; for (i_mk = 0; i_mk < rn_w; i_mk++ ) 	ptr_mk->list_f[i_mk+i] =  rl_w[i_mk];
	i += rn_w; 	  for (i_mk = 0; i_mk < rn_e; i_mk++ ) 	ptr_mk->list_f[i_mk+i] =  rl_e[i_mk];
	i += rn_e; 	  for (i_mk = 0; i_mk < rn_s; i_mk++ ) 	ptr_mk->list_f[i_mk+i] =  rl_s[i_mk];
	i += rn_s; 	  for (i_mk = 0; i_mk < rn_n; i_mk++ ) 	ptr_mk->list_f[i_mk+i] =  rl_n[i_mk];

	//for(i=0; i < ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f; i++)
		//printf("rank %d i %d %d\n",rank,i,ptr_mk->list_f[i]);

	VecDestroy(&ptr_mk->sp_x); VecDestroy(&ptr_mk->sp_y); VecDestroy(&ptr_mk->ep_x); VecDestroy(&ptr_mk->ep_y);
	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,&ptr_mk->sp_x);
	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,&ptr_mk->sp_y);
	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,&ptr_mk->ep_x);
	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,&ptr_mk->ep_y);

	ISCreateGeneral(PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,ptr_mk->list_f,PETSC_COPY_VALUES,&is2);
	ISCreateStride (PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,0,1,&is1);
	VecScatterCreate(ptr_mk->sp_x0,is2,ptr_mk->sp_x,is1,&ptr_mk->glo2loc);

	VecScatterBegin(ptr_mk->glo2loc ,dis_sp_x,ptr_mk->sp_x,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd  (ptr_mk->glo2loc ,dis_sp_x,ptr_mk->sp_x,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterBegin(ptr_mk->glo2loc ,dis_sp_y,ptr_mk->sp_y,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd  (ptr_mk->glo2loc ,dis_sp_y,ptr_mk->sp_y,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterBegin(ptr_mk->glo2loc ,dis_ep_x,ptr_mk->ep_x,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd  (ptr_mk->glo2loc ,dis_ep_x,ptr_mk->ep_x,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterBegin(ptr_mk->glo2loc ,dis_ep_y,ptr_mk->ep_y,INSERT_VALUES,SCATTER_FORWARD);
	VecScatterEnd  (ptr_mk->glo2loc ,dis_ep_y,ptr_mk->ep_y,INSERT_VALUES,SCATTER_FORWARD);

	// dealing with accelerations and velocities
	VecDestroy(&ptr_mk->sp_u); VecDestroy(&ptr_mk->sp_v);
	VecDestroy(&ptr_mk->ep_u); VecDestroy(&ptr_mk->ep_v);
	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,&ptr_mk->sp_u);
	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,&ptr_mk->sp_v);
	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,&ptr_mk->ep_u);
	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,&ptr_mk->ep_v);

	VecDestroy(&ptr_mk->sp_ax); VecDestroy(&ptr_mk->sp_ay);
	VecDestroy(&ptr_mk->ep_ax); VecDestroy(&ptr_mk->ep_ay);
	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,&ptr_mk->sp_ax);
	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,&ptr_mk->sp_ay);
	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,&ptr_mk->ep_ax);
	VecCreateSeq	(PETSC_COMM_SELF,ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f,&ptr_mk->ep_ay);

    for(elem = 1; elem <3; elem ++){ // displacement, velocities, acceleration

    	switch(elem)
    	{
    	case 1:
    		VecScatterBegin(ptr_u->scatter_x,ptr_s->dxi,xix_local,INSERT_VALUES,SCATTER_FORWARD);
    	    VecScatterEnd  (ptr_u->scatter_x,ptr_s->dxi,xix_local,INSERT_VALUES,SCATTER_FORWARD);

    	    VecScatterBegin(ptr_u->scatter_y,ptr_s->dxi,xiy_local,INSERT_VALUES,SCATTER_FORWARD);
    	    VecScatterEnd  (ptr_u->scatter_y,ptr_s->dxi,xiy_local,INSERT_VALUES,SCATTER_FORWARD);
    	    break;

    	case 2:
    		VecScatterBegin(ptr_u->scatter_x,ptr_s->ddxi,xix_local,INSERT_VALUES,SCATTER_FORWARD);
    	    VecScatterEnd  (ptr_u->scatter_x,ptr_s->ddxi,xix_local,INSERT_VALUES,SCATTER_FORWARD);

    	    VecScatterBegin(ptr_u->scatter_y,ptr_s->ddxi,xiy_local,INSERT_VALUES,SCATTER_FORWARD);
    	    VecScatterEnd  (ptr_u->scatter_y,ptr_s->ddxi,xiy_local,INSERT_VALUES,SCATTER_FORWARD);
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
		case 1:

			VecScatterBegin	(ptr_mk->glo2loc,dis_sp_x,ptr_mk->sp_u,INSERT_VALUES,SCATTER_FORWARD);
			VecScatterEnd	(ptr_mk->glo2loc,dis_sp_x,ptr_mk->sp_u,INSERT_VALUES,SCATTER_FORWARD);

			VecScatterBegin	(ptr_mk->glo2loc,dis_sp_y,ptr_mk->sp_v,INSERT_VALUES,SCATTER_FORWARD);
			VecScatterEnd	(ptr_mk->glo2loc,dis_sp_y,ptr_mk->sp_v,INSERT_VALUES,SCATTER_FORWARD);

			VecScatterBegin	(ptr_mk->glo2loc,dis_ep_x,ptr_mk->ep_u,INSERT_VALUES,SCATTER_FORWARD);
			VecScatterEnd	(ptr_mk->glo2loc,dis_ep_x,ptr_mk->ep_u,INSERT_VALUES,SCATTER_FORWARD);

			VecScatterBegin	(ptr_mk->glo2loc,dis_ep_y,ptr_mk->ep_v,INSERT_VALUES,SCATTER_FORWARD);
			VecScatterEnd	(ptr_mk->glo2loc,dis_ep_y,ptr_mk->ep_v,INSERT_VALUES,SCATTER_FORWARD);
			break;

		case 2:

			VecScatterBegin	(ptr_mk->glo2loc,dis_sp_x,ptr_mk->sp_ax,INSERT_VALUES,SCATTER_FORWARD);
			VecScatterEnd	(ptr_mk->glo2loc,dis_sp_x,ptr_mk->sp_ax,INSERT_VALUES,SCATTER_FORWARD);

			VecScatterBegin	(ptr_mk->glo2loc,dis_sp_y,ptr_mk->sp_ay,INSERT_VALUES,SCATTER_FORWARD);
			VecScatterEnd	(ptr_mk->glo2loc,dis_sp_y,ptr_mk->sp_ay,INSERT_VALUES,SCATTER_FORWARD);

			VecScatterBegin	(ptr_mk->glo2loc,dis_ep_x,ptr_mk->ep_ax,INSERT_VALUES,SCATTER_FORWARD);
			VecScatterEnd	(ptr_mk->glo2loc,dis_ep_x,ptr_mk->ep_ax,INSERT_VALUES,SCATTER_FORWARD);

			VecScatterBegin	(ptr_mk->glo2loc,dis_ep_y,ptr_mk->ep_ay,INSERT_VALUES,SCATTER_FORWARD);
			VecScatterEnd	(ptr_mk->glo2loc,dis_ep_y,ptr_mk->ep_ay,INSERT_VALUES,SCATTER_FORWARD);
			break;
		}
    }

	VecDestroy(&dis_sp_x); 	VecDestroy(&dis_sp_y);
	VecDestroy(&dis_ep_x); 	VecDestroy(&dis_ep_y);
	VecDestroy(&xix_local); VecDestroy(&xiy_local);
	ISDestroy (&is1);
	ISDestroy (&is2);
	free_ivector(mk_int, 0, ptr_mk->nglobal-1);
	free_ivector(sl_e, 0, ptr_mk->nglobal-1);
	free_ivector(sl_w, 0, ptr_mk->nglobal-1);
	free_ivector(sl_s, 0, ptr_mk->nglobal-1);
	free_ivector(sl_n, 0, ptr_mk->nglobal-1);

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
	int 		loc_i,panel_i,bool1, bool2,s,i,j;
	const int	search_range = 0;

	// get the coordinates of the markers
	VecGetArray(ptr_mk->sp_x, &sx);
	VecGetArray(ptr_mk->sp_y, &sy);
	VecGetArray(ptr_mk->ep_x, &ex);
	VecGetArray(ptr_mk->ep_y, &ey);

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
	for (loc_i=0; loc_i<ptr_mk->nlocal_f + ptr_mk->nlocal_gho_f; loc_i++)
	{
		x1 = sx[loc_i]; x2 = ex[loc_i];
		y1 = sy[loc_i]; y2 = ey[loc_i];

		//printf("rank %d i %d %f %f %f %f\n", ptr_u->rank,loc_i,x1,y1,x2,y2);
		//panel_i = ptr_mk->list_f[loc_i];

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

		if  (x1!=x2){// if it is not a vertical panel
			nc[0]=be_m; nc[1]=1;
			s=SIGN(1,n[0]*nc[0]+n[1]*nc[1]);
			ptr_mk->nx[loc_i] = s*be_m/be_d;
			ptr_mk->ny[loc_i] = s/be_d;
		}
		else {
			nc[0]=1; nc[1]=0;
			s=SIGN(1,n[0]*nc[0]+n[1]*nc[1]);
			ptr_mk->nx[loc_i] = SIGN(1,n[0]);
			ptr_mk->ny[loc_i] = 0;
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

		//printf("rank %d i %d search range %d %d %d %d\n", ptr_u->rank,loc_i,search_x[0],search_x[1],search_y[0],search_y[1] );
		for (i = search_x[0]; i<=search_x[1]; i++)
			for (j = search_y[0]; j<=search_y[1]; j++)
		    {
				gridx = ptr_f->xm_grid[i];   gridy = ptr_f->ym_grid[j];

				if (x1!=x2){
		            d  = (be_m*gridx + gridy +be_c)/be_d;
		            xb = gridx-d/be_d*be_m;
		            yb = gridy-d/be_d;
		            c  = (xb-x1)/(x2-x1);
				}
				else{
		            d = gridx-x1;
		            xb = x1;  yb = gridy;
		            c = (yb-y1)/(y2-y1);
				}

				if (c >=0 && c<=1){
					//if(ptr_u->rank==0) printf("i,j,s*d dist %d %d %f %f\n",i,j,fabs(s*d),fabs(dist_p[j][i]));
					if((panel_p[j][i]==0) || (fabs(s*d)<fabs(dist_p[j][i])))
					{
						dist_p[j][i] 	= s*d;
						panel_p[j][i] 	= loc_i+1;   // The local panel indexing begins from 1 here!!!!!
						ratio_p[j][i] 	= c;
					}

				}
				else if (c<0){
					dd = sqrt((gridx-x1)*(gridx-x1)+(gridy-y1)*(gridy-y1));
					if((panel_p[j][i]==0) || ( dd < fabs(dist_p[j][i])))
					{
						dist_p[j][i] 	= 	SIGN(1,s*d)*dd;
						panel_p[j][i] 	= -	(loc_i+1);
						ratio_p[j][i] 	= 	0;
					}
				}
				else{
					dd = sqrt((gridx-x2)*(gridx-x2)+(gridy-y2)*(gridy-y2));
					if((panel_p[j][i]==0) || (dd < fabs(dist_p[j][i]))){
						dist_p[j][i] 	= 	SIGN(1,s*d)*dd;
						panel_p[j][i] 	= -	(loc_i+1);
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
						panel_u[j][i] 	= loc_i+1;   // The panel numbering begins from 1 here!!!!!
						ratio_u[j][i] 	= c;
					}

				}
				else if (c<0)
				{
					dd = sqrt((gridx-x1)*(gridx-x1)+(gridy-y1)*(gridy-y1));
					if((panel_u[j][i]==0) || (fabs(dd)<fabs(dist_u[j][i])))
					{
						dist_u[j][i] 	= 	SIGN(1,s*d)*dd;
						panel_u[j][i] 	= -	(loc_i+1);
						ratio_u[j][i] 	= 	0;
					}
				}
				else
				{
					dd = sqrt((gridx-x2)*(gridx-x2)+(gridy-y2)*(gridy-y2));
					if((panel_u[j][i]==0) || (fabs(dd)<fabs(dist_u[j][i])))
					{
						dist_u[j][i] 	= 	SIGN(1,s*d)*dd;
						panel_u[j][i] 	= -	(loc_i+1);
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
						panel_v[j][i] 	= loc_i+1;   // The panel numbering begins from 1 here!!!!!
						ratio_v[j][i] 	= c;
					}
				}
				else if (c<0){

					dd = sqrt((gridx-x1)*(gridx-x1)+(gridy-y1)*(gridy-y1));
					if((panel_v[j][i]==0) || (abs(dd)<abs(dist_v[j][i])))
					{
						dist_v[j][i] 	= 	SIGN(1,s*d)*dd;
						panel_v[j][i] 	= -	(loc_i+1);
						ratio_v[j][i] 	= 	0;
					}
				}
				else {

					dd = sqrt((gridx-x2)*(gridx-x2)+(gridy-y2)*(gridy-y2));
					if((panel_v[j][i]==0) || (fabs(dd)<fabs(dist_v[j][i])))
					{
						dist_v[j][i] 	= 	SIGN(1,s*d)*dd;
						panel_v[j][i] 	= -	(loc_i+1);
						ratio_v[j][i] 	= 	1;
					}
				}
		    }

	}

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

	            	DMGlobalToLocalBegin(da, ptr_f->dist_u,INSERT_VALUES,local_dist);
	            	DMGlobalToLocalEnd	(da, ptr_f->dist_u,INSERT_VALUES,local_dist);
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
	                printf("!!! (%d,%d) has ONLY %d interpolation template points\n",
	                       i,j,count);

	            // sort offset points using their distances to the forcing point
	            qsort(dist_list, count, sizeof(dist_t), (compfn)comp_distance);

	            // Pick 2 closest points, storing them in data[id].tp[].
	            for ( b=0; b<2; b++) {
	                list->data[m].tp[b] = dist_list[b].i;
	            }

	            // calculate the coordinate of the surface point whose normal
	            // vector passes the forcing point.

	            pn = fabs(panel[j][i])-1; //
	            rt = ratio[j][i];

	            list->data[m].surf[0]  = rt*(ex[pn]-sx[pn])+sx[pn];
	            list->data[m].surf[1]  = rt*(ey[pn]-sy[pn])+sy[pn];
	            //list->data[m].panel 	= ptr_mk->list_f[pn];
	            //list->data[m].ratio 	= rt;

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

	        // construct global forcing point map from local map
	        switch (comp)
	        {
	        case 1: // u component

	        	DMLocalToGlobalBegin(da,dist_map_ghosted,INSERT_VALUES,ptr_f->dist_map_u);
	        	DMLocalToGlobalEnd	(da,dist_map_ghosted,INSERT_VALUES,ptr_f->dist_map_u);
	        	break;

	        case 2: // v component

	        	DMLocalToGlobalBegin(da,dist_map_ghosted,INSERT_VALUES,ptr_f->dist_map_v);
	        	DMLocalToGlobalEnd	(da,dist_map_ghosted,INSERT_VALUES,ptr_f->dist_map_v);
	        	break;

	        case 3: // p component

	        	DMLocalToGlobalBegin(da,dist_map_ghosted,INSERT_VALUES,ptr_f->dist_map_p);
	        	DMLocalToGlobalEnd	(da,dist_map_ghosted,INSERT_VALUES,ptr_f->dist_map_p);
	        	break;

	        } // end of switch
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

	const double  *sx, *sy, *ex,*ey,*panelu, *panelv, *panelax, *panelay, *l0;
	const double  *tx, *ty;
	int 		loc_panel, count, panel, id_t[2], rank;
	double		ref_pressure = 0, a[9], b[3], c[3];
	double		**u_ary, **v_ary, **p_ary;
	double		**fpm_u_ary, **fpm_v_ary, **fpm_p_ary;
	double		**dst_u_ary, **dst_v_ary, **dst_p_ary;
	double		mp[2], data[2][16], ux, uy, vx, vy, pp,lratio, trac[2],nx,ny;
	int 		pan_i,pan_j, i,j, ii1, ii2, jj1, jj2;
	dist_t2  	dist_list[16];
	Vec			panelu_vec, panelv_vec, panelax_vec, panelay_vec;
	Vec			u_vec, v_vec, p_vec;
	Vec			du_vec, dv_vec, dp_vec;
	Vec			fu_vec, fv_vec, fp_vec;

	PetscScalar cc1 = 1.0, cc2 = 0.5;

	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	VecDuplicate(ptr_mk->sp_u,&panelu_vec);
	VecDuplicate(ptr_mk->sp_u,&panelv_vec);
	VecDuplicate(ptr_mk->sp_u,&panelax_vec);
	VecDuplicate(ptr_mk->sp_u,&panelay_vec);

	VecWAXPY(panelu_vec, cc1, ptr_mk->sp_u, ptr_mk->ep_u);
	VecWAXPY(panelv_vec, cc1, ptr_mk->sp_v, ptr_mk->ep_v);
	VecWAXPY(panelax_vec, cc1, ptr_mk->sp_ax, ptr_mk->ep_ax);
	VecWAXPY(panelay_vec, cc1, ptr_mk->sp_ay, ptr_mk->ep_ay);

	VecScale(panelu_vec, cc2);
	VecScale(panelv_vec, cc2);
	VecScale(panelax_vec, cc2);
	VecScale(panelay_vec, cc2);

	VecGetArrayRead(panelu_vec, &panelu);
	VecGetArrayRead(panelv_vec, &panelv);
	VecGetArrayRead(panelax_vec, &panelax);
	VecGetArrayRead(panelay_vec, &panelay);

	VecGetArrayRead(ptr_mk->sp_x, &sx);
	VecGetArrayRead(ptr_mk->sp_y, &sy);
	VecGetArrayRead(ptr_mk->ep_x, &ex);
	VecGetArrayRead(ptr_mk->ep_y, &ey);

	VecGetArrayRead(ptr_mk->l0, &l0);
	VecGetArrayRead(ptr_mk->tx, &tx);
	VecGetArrayRead(ptr_mk->ty, &ty);

	DMGetLocalVector(ptr_f->da, &u_vec);
    DMGlobalToLocalBegin(ptr_f->da, ptr_f->u, INSERT_VALUES, u_vec);
    DMGlobalToLocalEnd  (ptr_f->da, ptr_f->u, INSERT_VALUES, u_vec);
    DMDAVecGetArrayRead (ptr_f->da, u_vec, &u_ary);

	DMGetLocalVector(ptr_f->da, &v_vec);
    DMGlobalToLocalBegin(ptr_f->da, ptr_f->v, INSERT_VALUES, v_vec);
    DMGlobalToLocalEnd  (ptr_f->da, ptr_f->v, INSERT_VALUES, v_vec);
    DMDAVecGetArrayRead (ptr_f->da, v_vec, &v_ary);

	DMGetLocalVector(ptr_f->da, &p_vec);
    DMGlobalToLocalBegin(ptr_f->da, ptr_f->p, INSERT_VALUES, p_vec);
    DMGlobalToLocalEnd  (ptr_f->da, ptr_f->p, INSERT_VALUES, p_vec);
    DMDAVecGetArrayRead (ptr_f->da, p_vec, &p_ary);

	DMGetLocalVector(ptr_f->da, &du_vec);
    DMGlobalToLocalBegin(ptr_f->da, ptr_f->dist_map_u, INSERT_VALUES, du_vec);
    DMGlobalToLocalEnd  (ptr_f->da, ptr_f->dist_map_u, INSERT_VALUES, du_vec);
    DMDAVecGetArrayRead (ptr_f->da, du_vec, &fpm_u_ary);

	DMGetLocalVector(ptr_f->da, &dv_vec);
    DMGlobalToLocalBegin(ptr_f->da, ptr_f->dist_map_v, INSERT_VALUES, dv_vec);
    DMGlobalToLocalEnd  (ptr_f->da, ptr_f->dist_map_v, INSERT_VALUES, dv_vec);
    DMDAVecGetArrayRead (ptr_f->da, dv_vec, &fpm_v_ary);

	DMGetLocalVector(ptr_f->da, &dp_vec);
    DMGlobalToLocalBegin(ptr_f->da, ptr_f->dist_map_p, INSERT_VALUES, dp_vec);
    DMGlobalToLocalEnd  (ptr_f->da, ptr_f->dist_map_p, INSERT_VALUES, dp_vec);
    DMDAVecGetArrayRead (ptr_f->da, dp_vec, &fpm_p_ary);

	DMGetLocalVector(ptr_f->da, &fu_vec);
    DMGlobalToLocalBegin(ptr_f->da, ptr_f->dist_u, INSERT_VALUES, fu_vec);
    DMGlobalToLocalEnd  (ptr_f->da, ptr_f->dist_u, INSERT_VALUES, fu_vec);
    DMDAVecGetArrayRead (ptr_f->da, fu_vec, &dst_u_ary);

	DMGetLocalVector(ptr_f->da, &fv_vec);
    DMGlobalToLocalBegin(ptr_f->da, ptr_f->dist_v, INSERT_VALUES, fv_vec);
    DMGlobalToLocalEnd  (ptr_f->da, ptr_f->dist_v, INSERT_VALUES, fv_vec);
    DMDAVecGetArrayRead (ptr_f->da, fv_vec, &dst_v_ary);

	DMGetLocalVector(ptr_f->da, &fp_vec);
    DMGlobalToLocalBegin(ptr_f->da, ptr_f->dist_p, INSERT_VALUES, fp_vec);
    DMGlobalToLocalEnd  (ptr_f->da, ptr_f->dist_p, INSERT_VALUES, fp_vec);
    DMDAVecGetArrayRead (ptr_f->da, fp_vec, &dst_p_ary);

	VecZeroEntries(ptr_s->traction);

	for( loc_panel = 0; loc_panel<ptr_mk->nlocal_f;loc_panel++)
	{
		mp[0] = (sx[loc_panel] + ex[loc_panel])/2; //midpoint of the panel
        mp[1] = (sy[loc_panel] + ey[loc_panel])/2;
        panel = ptr_mk->list_f[loc_panel];

        lratio = sqrt((sx[loc_panel]-ex[loc_panel])*(sx[loc_panel]-ex[loc_panel])
      			 +(sy[loc_panel]-ey[loc_panel])*(sy[loc_panel]-ey[loc_panel]))/l0[panel];
        nx = ptr_mk->nx[loc_panel];
        ny = ptr_mk->ny[loc_panel];
        /*
        // u gradient
        pan_i = floor((mp[0]-ptr_f->x_grid[0])/ptr_f->dx );
        pan_j = floor((mp[1]-ptr_f->y_grid[0])/ptr_f->dx -0.5);

        count = 0;
        for (i = pan_i-1; i<= pan_i+2; i++)
        	 for (j = pan_j-1; j<= pan_j+2; j++)
        	 {
        		 printf("%d %d %d %d\n",rank,panel,i,j);
        		 if (dst_u_ary[j][i]>0 && fpm_u_ary[j][i] == 0){
        	    	 dist_list[count].i = i;
        	    	 dist_list[count].j = j;
        	       	 dist_list[count].d = sqrt((ptr_f->x_grid[i]-mp[0])*(ptr_f->x_grid[i]-mp[0])
        	       			 +(ptr_f->ym_grid[j]-mp[1])*(ptr_f->ym_grid[j]-mp[1]));
        	         count ++;
        	     }
        		 printf("%d\n",count);
        	 }

        qsort(dist_list, count, sizeof(dist_t), (compfn)comp_distance);

        ii1 = dist_list[0].i;  jj1 = dist_list[0].j;
        ii2 = dist_list[1].i;  jj2 = dist_list[1].j;

        //choose two nearest fluid points
        a[0] = 1; a[3] = mp[0];  a[6] = mp[1];
        a[1] = 1; a[4] = ptr_f->x_grid[ii1];  a[7] = ptr_f->ym_grid[jj1];
        a[2] = 1; a[5] = ptr_f->x_grid[ii2];  a[8] = ptr_f->ym_grid[jj2];

        PetscKernel_A_gets_inverse_A_3(a,0.0);

        b[0] = panelu[loc_panel];
        b[1] = u_ary[jj1][ii1];
        b[2] = u_ary[jj2][ii2];

        A_times_b_3x3(a,b,c);
        ux = c[1]; uy = c[2];


        // v gradient
        pan_i = floor((mp[0]-ptr_f->x_grid[0])/ptr_f->dx -0.5);
        pan_j = floor((mp[1]-ptr_f->y_grid[0])/ptr_f->dx );

        printf("%d %d %d %d\n",rank,  panel, pan_i,pan_j);

        count = 0;
        for (i = pan_i-1; i<= pan_i+2; i++)
        	 for (j = pan_j-1; j<= pan_j+2; j++)
        	 {
        	     if (dst_v_ary[j][i]>0 && fpm_v_ary[j][i] == 0){
        	    	 dist_list[count].i = i;
        	    	 dist_list[count].j = j;
        	       	 dist_list[count].d = sqrt((ptr_f->xm_grid[i]-mp[0])*(ptr_f->xm_grid[i]-mp[0])
        	       			 +(ptr_f->y_grid[j]-mp[1])*(ptr_f->y_grid[j]-mp[1]));
        	         count ++;
        	     }
        	 }

        qsort(dist_list, count, sizeof(dist_t), (compfn)comp_distance);

        ii1 = dist_list[0].i;  jj1 = dist_list[0].j;
        ii2 = dist_list[1].i;  jj2 = dist_list[1].j;
        //choose two nearest fluid points
        a[0] = 1; a[3] = mp[0];  a[6] = mp[1];
        a[1] = 1; a[4] = ptr_f->xm_grid[ii1];  a[7] = ptr_f->y_grid[jj1];
        a[2] = 1; a[5] = ptr_f->xm_grid[ii2];  a[8] = ptr_f->y_grid[jj2];

        PetscKernel_A_gets_inverse_A_3(a,0.0);

        b[0] = panelv[loc_panel];
        b[1] = v_ary[jj1][ii1];
        b[2] = v_ary[jj2][ii2];

        A_times_b_3x3(a,b,c);
        vx = c[1]; vy = c[2];
        */
        if (rank == 0){
        // p gradient
        pan_i = floor((mp[0]-ptr_f->x_grid[0])/ptr_f->dx -0.5);
        pan_j = floor((mp[1]-ptr_f->y_grid[0])/ptr_f->dx -0.5);

        count = 0;
        for (i = pan_i-1; i<= pan_i+2; i++)
        	 for (j = pan_j-1; j<= pan_j+2; j++)
        	 {
        	     printf("rank %d %d %d %d\n", rank, panel,i,j );
        		 if (dst_p_ary[j][i]>0 && fpm_p_ary[j][i] == 0){
        	    	 dist_list[count].i = i;
        	    	 dist_list[count].j = j;
        	       	 dist_list[count].d = sqrt((ptr_f->xm_grid[i]-mp[0])*(ptr_f->xm_grid[i]-mp[0])
        	       			 +(ptr_f->ym_grid[j]-mp[1])*(ptr_f->ym_grid[j]-mp[1]));
        	         count ++;
        	     }
        	 }

        qsort(dist_list, count, sizeof(dist_t), (compfn)comp_distance);

        ii1 = dist_list[0].i;  jj1 = dist_list[0].j;
        ii2 = dist_list[1].i;  jj2 = dist_list[1].j;

        // choose two nearest fluid points
        a[0] = 0; a[3] = nx;  a[6] = ny;
        a[1] = 1; a[4] = ptr_f->xm_grid[ii1];  a[7] = ptr_f->ym_grid[jj1];
        a[2] = 1; a[5] = ptr_f->xm_grid[ii2];  a[8] = ptr_f->ym_grid[jj2];


        PetscKernel_A_gets_inverse_A_3(a,0.0);

        b[0] = -panelax[loc_panel]*nx -panelay[loc_panel]*ny ;
        b[1] = p_ary[jj1][ii1];
        b[2] = p_ary[jj2][ii2];

        A_times_b_3x3(a,b,c);
        pp = c[0] + mp[0]*c[1]+mp[1]*c[2] - ref_pressure;

        // calculating traction
        id_t[0] = (int)tx[panel];
        id_t[1] = (int)ty[panel];

        printf(" %d %d\n", id_t[0], id_t[1]);

        /*
        trac[0] = (-pp*nx + (nx*2*ux + ny*( vx + uy)/ptr_f->reynolds_number))*lratio;
        trac[1] = (-pp*ny + (ny*2*vy + nx*( vx + uy)/ptr_f->reynolds_number))*lratio;

        //VecSetValues(ptr_s->traction,2,id_t,trac, ADD_VALUES);
		*/
        }
	}

	VecAssemblyBegin(ptr_s->traction);
	VecAssemblyEnd	(ptr_s->traction);

	VecRestoreArrayRead(ptr_mk->l0, &l0);
	VecRestoreArrayRead(ptr_mk->tx, &tx);
	VecRestoreArrayRead(ptr_mk->ty, &ty);

	DMDAVecRestoreArrayRead(ptr_f->da,u_vec, &u_ary);
	DMDAVecRestoreArrayRead(ptr_f->da,v_vec, &v_ary);
	DMDAVecRestoreArrayRead(ptr_f->da,p_vec, &p_ary);

	DMDAVecRestoreArrayRead(ptr_f->da,du_vec, &fpm_u_ary);
	DMDAVecRestoreArrayRead(ptr_f->da,dv_vec, &fpm_v_ary);
	DMDAVecRestoreArrayRead(ptr_f->da,dp_vec, &fpm_p_ary);

	DMDAVecRestoreArrayRead(ptr_f->da,fu_vec, &dst_u_ary);
	DMDAVecRestoreArrayRead(ptr_f->da,fv_vec, &dst_v_ary);
	DMDAVecRestoreArrayRead(ptr_f->da,fp_vec, &dst_p_ary);

	VecRestoreArrayRead(panelu_vec, &panelu);
	VecRestoreArrayRead(panelv_vec, &panelv);
	VecRestoreArrayRead(panelax_vec, &panelax);
	VecRestoreArrayRead(panelay_vec, &panelay);

	VecRestoreArrayRead(ptr_mk->sp_ax, &sx);
	VecRestoreArrayRead(ptr_mk->ep_ax, &ex);
	VecRestoreArrayRead(ptr_mk->sp_ay, &sy);
	VecRestoreArrayRead(ptr_mk->ep_ay, &ey);

	VecDestroy(&panelu_vec); VecDestroy(&panelv_vec); VecDestroy(&panelax_vec); VecDestroy(&panelay_vec);

}

int forcing_pts_list_init(forcing_point_list_t *list)
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

void disp_interp_setup(Lag_marker* ptr_mk, Index_S* ptr_i, Grid_S* ptr_g ){ // 11/13/2015 by Peggy

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


