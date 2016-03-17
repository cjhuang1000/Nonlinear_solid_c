#include "fluid.h"

static int index_range_set_raw_ghost(char which_var, Field_F *f,
                                     const int *i, const int *j, const int *k);
static int index_range_set_raw(char which_var, Field_F *f,
                               const int *i, const int *j, const int *k);
static int index_range_set_rhs(char which_var, Field_F *f,
                               const int *i, const int *j, const int *k);


void fluid_setup(Field_F* f, Field_S* s)
{

	int  			s_width 	= 2; // stencil width
	int  			dof 		= 1;
	int  			nx 			= (s->Nx-1)*2+1, ny = (s->Ny-1)*2+1;  // fluid grid size is half of solid grid size
	int				i,rank;
	const int 		start		= 0, end=1;

	double 			dx			= s->dx/2;
	double			domain[6]	={-1.0,1.0,-1.0,1.0, 0.0, 0.0};
	DMDAStencilType s_type 		= DMDA_STENCIL_BOX;
	DMBoundaryType  b_type 		= DM_BOUNDARY_NONE;

	PetscReal		**dist_u,**dist_v,**dist_p;

	// Setup grid of fluid

	// Set the DMDA
	DMDACreate2d(MPI_COMM_WORLD, b_type, b_type, s_type, nx, ny,
	             PETSC_DECIDE, PETSC_DECIDE, dof, s_width,
	             PETSC_NULL, PETSC_NULL,
	             &f->da);


	// Set subdomain info, index range info
	int NX, NY, NZ, npx, npy, npz;
	DMDAGetInfo(f->da, PETSC_NULL,                  /* DA, no. of dims */
	            &NX, &NY, &NZ,                      /* global dimensions */
	            &npx, &npy, &npz,                   /* proc grid dimensions */
	            PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, PETSC_NULL, /* dof, width, boudanrytype */
	            PETSC_NULL);                        /* stencil type */

	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	f->proc_coords[0] = rank % npx;
	f->proc_coords[2] = rank / (npx * npy);
	f->proc_coords[1] = (rank - npx*npy*f->proc_coords[2]) / npx;
	f->dx = dx;

	const int px = f->proc_coords[0];
	const int py = f->proc_coords[1];
	const int pz = f->proc_coords[2];

	// ------------------ set up index ---------------

	int n1,n2,n3,i0,j0,k0,idl[2], jdl[2], kdl[2];
	DMDAGetCorners(f->da, &i0, &j0, &k0, &n1, &n2, &n3);

	idl[start] = i0; idl[end] = i0+n1-1;
	jdl[start] = j0; jdl[end] = j0+n2-1;
	kdl[start] = k0; kdl[end] = k0+n3-1;

	index_range_set_raw('u', f, idl, jdl, kdl);
	index_range_set_raw('v', f, idl, jdl, kdl);
	index_range_set_raw('p', f, idl, jdl, kdl);

	f->idx_range.int_bound[0] = (idl[start])*f->dx + domain[0];
	f->idx_range.int_bound[1] = (idl[end]+1)*f->dx + domain[0];
	f->idx_range.int_bound[2] = (jdl[start])*f->dx + domain[2];
	f->idx_range.int_bound[3] = (jdl[end]+1)*f->dx + domain[2];
	f->idx_range.int_bound[4] = (kdl[start])*f->dx + domain[4];
	f->idx_range.int_bound[5] = (kdl[end]+1)*f->dx + domain[4];

	// @@ Jan 2016, for panel redistribution
	idl[start] += 1;  idl[end] -= 1;
	jdl[start] += 1;  jdl[end] -= 1;
	kdl[start] += 1;  kdl[end] -= 1;
	//idl[start] += s_width;  idl[end] -= s_width;
	//jdl[start] += s_width;  jdl[end] -= s_width;
	//kdl[start] += s_width;  kdl[end] -= s_width;

	if (px == 0    ) idl[start] = 0;
	if (px == npx-1) idl[end]   = NX-1;

	if (py == 0    ) jdl[start] = 0;
	if (py == npy-1) jdl[end]   = NY-1;

	if (pz == 0    ) kdl[start] = 0;
	if (pz == npz-1) kdl[end]   = NZ-1;

	f->idx_range.shared_bound[0] = (idl[start])*f->dx + domain[0];
	f->idx_range.shared_bound[1] = (idl[end]+1)*f->dx + domain[0];
	f->idx_range.shared_bound[2] = (jdl[start])*f->dx + domain[2];
	f->idx_range.shared_bound[3] = (jdl[end]+1)*f->dx + domain[2];
	f->idx_range.shared_bound[4] = (kdl[start])*f->dx + domain[4];
	f->idx_range.shared_bound[5] = (kdl[end]+1)*f->dx + domain[4];

	// Neighboring processor's rank. If there is no neighboring -> -1.

	// left
	if (px == 0) 	 f->idx_range.shared_pro[0] = -1;
	else 		 	 f->idx_range.shared_pro[0] = pz*npx*npy + py*npx + px -1;

	// right
	if (px == npx-1) f->idx_range.shared_pro[1] = -1;
	else 		 	 f->idx_range.shared_pro[1] = pz*npx*npy + py*npx + px +1;

	// bottom
	if (py == 0) 	 f->idx_range.shared_pro[2] = -1;
	else 		 	 f->idx_range.shared_pro[2] = pz*npx*npy + (py-1)*npx + px;

	// top
	if (py == npy-1) f->idx_range.shared_pro[3] = -1;
	else 		 	 f->idx_range.shared_pro[3] = pz*npx*npy + (py+1)*npx + px;

	// front and back (for 3D)
	if (pz == 0) 	 f->idx_range.shared_pro[4] = -1;
	else 		 	 f->idx_range.shared_pro[4] = (pz-1)*npx*npy + py*npx + px;

	if (pz == npz-1) f->idx_range.shared_pro[5] = -1;
	else 		 	 f->idx_range.shared_pro[5] = (pz+1)*npx*npy + py*npx + px;
	// @@

	int m1,m2,m3;

	DMDAGetGhostCorners(f->da,&i0,&j0,&k0,&m1, &m2, &m3);

	// raw ghosted ranges
	idl[start] = i0; idl[end] = i0+m1-1;
	jdl[start] = j0; jdl[end] = j0+m2-1;
	kdl[start] = k0; kdl[end] = k0+m3-1;

	index_range_set_raw_ghost('u', f, idl, jdl, kdl);
	index_range_set_raw_ghost('v', f, idl, jdl, kdl);
	index_range_set_raw_ghost('p', f, idl, jdl, kdl);

	// local rhs range
	DMDAGetCorners(f->da,&i0,&j0,&k0,&n1,&n2,&n3);

	idl[start] = i0; idl[end] = i0+n1-1;
	jdl[start] = j0; jdl[end] = j0+n2-1;
	kdl[start] = k0; kdl[end] = k0+n3-1;

	if (px == 0    ) idl[start] = 2;
	if (px == npx-1) idl[end]   = NX-2;

	if (py == 0    ) jdl[start] = 1;
	if (py == npy-1) jdl[end]   = NY-2;

	if (pz == 0    ) kdl[start] = 1;
	if (pz == npz-1) kdl[end]   = NZ-2;

	index_range_set_rhs('u',f,idl,jdl,kdl);

	idl[start] = i0; idl[end] = i0+n1-1;
	jdl[start] = j0; jdl[end] = j0+n2-1;
	kdl[start] = k0; kdl[end] = k0+n3-1;

	if (px == 0    ) idl[start] = 1;
	if (px == npx-1) idl[end]   = NX-2;

	if (py == 0    ) jdl[start] = 2;
	if (py == npy-1) jdl[end]   = NY-2;

    if (pz == 0    ) kdl[start] = 1;
    if (pz == npz-1) kdl[end]   = NZ-2;

	index_range_set_rhs('v',f,idl,jdl,kdl);

	idl[start] = i0; idl[end] = i0+n1-1;
	jdl[start] = j0; jdl[end] = j0+n2-1;
	kdl[start] = k0; kdl[end] = k0+n3-1;

    if (px==0    ) idl[start] = 1;
    if (px==npx-1) idl[end]   = NX-2;

    if (py==0    ) jdl[start] = 1;
    if (py==npy-1) jdl[end]   = NY-2;

    if (pz==0    ) kdl[start] = 2;
    if (pz==npz-1) kdl[end]   = NZ-2;

    index_range_set_rhs('w',f,idl,jdl,kdl);

    if (px == 0    ) idl[start] = 1;
    if (px == npx-1) idl[end]   = NX-2;

    if (py == 0    ) jdl[start] = 1;
    if (py == npy-1) jdl[end]   = NY-2;

    if (pz == 0    ) kdl[start] = 1;
    if (pz == npz-1) kdl[end]   = NZ-2;

    index_range_set_rhs('p',f,idl,jdl,kdl);

    // --------------- Create global velocity vectors

    DMCreateGlobalVector(f->da, &f->u    );
	DMCreateGlobalVector(f->da, &f->v    );
	DMCreateGlobalVector(f->da, &f->p    );

	VecSet(f->u ,0.05);
	VecSet(f->v ,0.05);
	VecSet(f->p ,0.0);

	// --------------- Create global sdf vector layout
	DMCreateGlobalVector(f->da, &f->dist_u     );
	DMCreateGlobalVector(f->da, &f->dist_v     );
	DMCreateGlobalVector(f->da, &f->dist_p     );

	DMCreateGlobalVector(f->da, &f->panel_u    );
	DMCreateGlobalVector(f->da, &f->panel_v    );
	DMCreateGlobalVector(f->da, &f->panel_p    );

	DMCreateGlobalVector(f->da, &f->ratio_u    );
	DMCreateGlobalVector(f->da, &f->ratio_v    );
	DMCreateGlobalVector(f->da, &f->ratio_p    );

	DMCreateGlobalVector(f->da, &f->dist_map_u     );
	DMCreateGlobalVector(f->da, &f->dist_map_v     );
	DMCreateGlobalVector(f->da, &f->dist_map_p     );

	// set grid c
	f->x_grid  = (double *) malloc(sizeof(double) * nx);
	f->y_grid  = (double *) malloc(sizeof(double) * ny);

	f->xm_grid = (double *) malloc(sizeof(double) * nx);
	f->ym_grid = (double *) malloc(sizeof(double) * ny);

	f->dx 	   = dx;
	// uniform grid size
	for(i=0; i<nx; i++)
	{
		f->x_grid[i]  = i*dx + domain[0];
		f->xm_grid[i] = (i+0.5)*dx + domain[0];
	}
	for(i=0; i<ny; i++)
	{
		f->y_grid[i]  = i*dx + domain[2];
		f->ym_grid[i] = (i+0.5)*dx + domain[2];
	}

	f->ncv[0]=nx;
	f->ncv[1]=ny;

	// initialize the sdf

	DMDAVecGetArray(f->da,	f->dist_u,	&dist_u);
	DMDAVecGetArray(f->da,	f->dist_v,	&dist_v);
	DMDAVecGetArray(f->da,	f->dist_p,	&dist_p);

	// @@
	set_sdf	(s, f->idx_range.u_raw[1][0], f->idx_range.u_raw[1][1]+1, f->idx_range.u_raw[0][0],
	       	  f->idx_range.u_raw[0][1]+1, f->ym_grid, f->x_grid, dist_u, 'f');
	set_sdf	(s, f->idx_range.v_raw[1][0], f->idx_range.v_raw[1][1]+1, f->idx_range.v_raw[0][0],
	       	  f->idx_range.v_raw[0][1]+1, f->y_grid, f->xm_grid, dist_v,'f');
	set_sdf	(s, f->idx_range.p_raw[1][0], f->idx_range.p_raw[1][1]+1, f->idx_range.p_raw[0][0],
	       	  f->idx_range.p_raw[0][1]+1, f->ym_grid, f->xm_grid, dist_p,'f');

	DMDAVecRestoreArray(f->da,	f->dist_u,	&dist_u);
	DMDAVecRestoreArray(f->da,	f->dist_v,	&dist_v);
	DMDAVecRestoreArray(f->da,	f->dist_p,	&dist_p);

	forcing_pts_list_init(&(f->flist_u));
	forcing_pts_list_init(&(f->flist_v));
	forcing_pts_list_init(&(f->flist_p));

}


// ------------------------------------------------------------------

static int index_range_set_raw(char which_var, Field_F *f,
                               const int *i, const int *j, const int *k)
{
  const int x=0, y=1, z=2;
  int rank,s;

  /* see also index_range_t in field.h */
  switch (which_var)
    {
    case 'u':
      for (s=0; s<=1; s++) {
        f->idx_range.u_raw[x][s] = i[s];
        f->idx_range.u_raw[y][s] = j[s];
        f->idx_range.u_raw[z][s] = k[s];
      }
      break;

    case 'v':
      for (s=0; s<=1; s++) {
        f->idx_range.v_raw[x][s] = i[s];
        f->idx_range.v_raw[y][s] = j[s];
        f->idx_range.v_raw[z][s] = k[s];
      }
      break;

    case 'w':
      for (s=0; s<=1; s++) {
        f->idx_range.w_raw[x][s] = i[s];
        f->idx_range.w_raw[y][s] = j[s];
        f->idx_range.w_raw[z][s] = k[s];
      }
      break;

    case 'p':
      for (s=0; s<=1; s++) {
        f->idx_range.p_raw[x][s] = i[s];
        f->idx_range.p_raw[y][s] = j[s];
        f->idx_range.p_raw[z][s] = k[s];
      }
      break;

    default:
      MPI_Comm_rank(f->comm, &rank);
      if (rank==0)
        printf("[%s] unknown variable type!\n",__func__);
      abort();
      break;
    } /* end of switch */
  return 0;
}

static int index_range_set_raw_ghost(char which_var, Field_F *f,
                                     const int *i, const int *j, const int *k)
{
  const int x=0, y=1, z=2;
  int rank,s;

  /* see also index_range_t in field.h */
  switch (which_var)
    {
    case 'u':
      for (s=0; s<=1; s++) {
        f->idx_range.u_raw_ghost[x][s] = i[s];
        f->idx_range.u_raw_ghost[y][s] = j[s];
        f->idx_range.u_raw_ghost[z][s] = k[s];
      }
      break;

    case 'v':
      for (s=0; s<=1; s++) {
        f->idx_range.v_raw_ghost[x][s] = i[s];
        f->idx_range.v_raw_ghost[y][s] = j[s];
        f->idx_range.v_raw_ghost[z][s] = k[s];
      }
      break;

    case 'w':
      for (s=0; s<=1; s++) {
        f->idx_range.w_raw_ghost[x][s] = i[s];
        f->idx_range.w_raw_ghost[y][s] = j[s];
        f->idx_range.w_raw_ghost[z][s] = k[s];
      }
      break;

    case 'p':
      for (s=0; s<=1; s++) {
        f->idx_range.p_raw_ghost[x][s] = i[s];
        f->idx_range.p_raw_ghost[y][s] = j[s];
        f->idx_range.p_raw_ghost[z][s] = k[s];
      }
      break;

    default:
      MPI_Comm_rank(f->comm, &rank);
      if (rank==0)
        printf("[%s] unknown variable type!\n",__func__);
      abort();
      break;
    } /* end of switch */
  return 0;
}

//
void index_range_get_raw(char which_var, Field_F *f,
                        int *i, int *j, int *k)
{
  const int x=0, y=1, z=2;
  int rank,s;

  /* see also index_range_t in field.h */

  switch (which_var)
    {
    case 'u':
      for ( s=0; s<=1; s++) {
        i[s] = f->idx_range.u_raw[x][s];
        j[s] = f->idx_range.u_raw[y][s];
        k[s] = f->idx_range.u_raw[z][s];
      }
      break;

    case 'v':
      for ( s=0; s<=1; s++) {
        i[s] = f->idx_range.v_raw[x][s];
        j[s] = f->idx_range.v_raw[y][s];
        k[s] = f->idx_range.v_raw[z][s];
      }
      break;

    case 'w':
      for ( s=0; s<=1; s++) {
        i[s] = f->idx_range.w_raw[x][s];
        j[s] = f->idx_range.w_raw[y][s];
        k[s] = f->idx_range.w_raw[z][s];
      }
      break;

    case 'p':
      for ( s=0; s<=1; s++) {
        i[s] = f->idx_range.p_raw[x][s];
        j[s] = f->idx_range.p_raw[y][s];
        k[s] = f->idx_range.p_raw[z][s];
      }
      break;

    } /* end of switch */
}

void index_range_get_raw_ghost(char which_var, Field_F *f,
                        int *i, int *j, int *k)
{
  const int x=0, y=1, z=2;
  int rank,s;

  /* see also index_range_t in field.h */

  switch (which_var)
    {
    case 'u':
      for ( s=0; s<=1; s++) {
        i[s] = f->idx_range.u_raw_ghost[x][s];
        j[s] = f->idx_range.u_raw_ghost[y][s];
        k[s] = f->idx_range.u_raw_ghost[z][s];
      }
      break;

    case 'v':
      for ( s=0; s<=1; s++) {
        i[s] = f->idx_range.v_raw_ghost[x][s];
        j[s] = f->idx_range.v_raw_ghost[y][s];
        k[s] = f->idx_range.v_raw_ghost[z][s];
      }
      break;

    case 'w':
      for ( s=0; s<=1; s++) {
        i[s] = f->idx_range.w_raw_ghost[x][s];
        j[s] = f->idx_range.w_raw_ghost[y][s];
        k[s] = f->idx_range.w_raw_ghost[z][s];
      }
      break;

    case 'p':
      for ( s=0; s<=1; s++) {
        i[s] = f->idx_range.p_raw_ghost[x][s];
        j[s] = f->idx_range.p_raw_ghost[y][s];
        k[s] = f->idx_range.p_raw_ghost[z][s];
      }
      break;

    } /* end of switch */
}


static int index_range_set_rhs(char which_var, Field_F *f,
                               const int *i, const int *j, const int *k)
{
  const int x=0, y=1, z=2;
  int rank,s;

  /* see also index_range_t in field.h */

  switch (which_var)
    {
    case 'u':
      for ( s=0; s<=1; s++) {
        f->idx_range.urhs[x][s] = i[s];
        f->idx_range.urhs[y][s] = j[s];
        f->idx_range.urhs[z][s] = k[s];
      }
      break;

    case 'v':
      for ( s=0; s<=1; s++) {
        f->idx_range.vrhs[x][s] = i[s];
        f->idx_range.vrhs[y][s] = j[s];
        f->idx_range.vrhs[z][s] = k[s];
      }
      break;

    case 'w':
      for ( s=0; s<=1; s++) {
        f->idx_range.wrhs[x][s] = i[s];
        f->idx_range.wrhs[y][s] = j[s];
        f->idx_range.wrhs[z][s] = k[s];
      }
      break;

    case 'p':
      for ( s=0; s<=1; s++) {
        f->idx_range.prhs[x][s] = i[s];
        f->idx_range.prhs[y][s] = j[s];
        f->idx_range.prhs[z][s] = k[s];
      }
      break;

    default:

        printf("[%s] unknown variable type!\n",__func__);

      break;
    } /* end of switch */
  return 0;
}


void index_range_get_rhs(char which_var, Field_F *f,
                        int *i, int *j, int *k)
{
  const int x=0, y=1, z=2;
  int rank,s;

  /* see also index_range_t in field.h */

  switch (which_var)
    {
    case 'u':
      for (s=0; s<=1; s++) {
        i[s] = f->idx_range.urhs[x][s];
        j[s] = f->idx_range.urhs[y][s];
        k[s] = f->idx_range.urhs[z][s];
      }
      break;

    case 'v':
      for (s=0; s<=1; s++) {
        i[s] = f->idx_range.vrhs[x][s];
        j[s] = f->idx_range.vrhs[y][s];
        k[s] = f->idx_range.vrhs[z][s];
      }
      break;

    case 'w':
      for (s=0; s<=1; s++) {
        i[s] = f->idx_range.wrhs[x][s];
        j[s] = f->idx_range.wrhs[y][s];
        k[s] = f->idx_range.wrhs[z][s];
      }
      break;

    case 'p':
      for (s=0; s<=1; s++) {
        i[s] = f->idx_range.prhs[x][s];
        j[s] = f->idx_range.prhs[y][s];
        k[s] = f->idx_range.prhs[z][s];
      }
      break;

    default:

        printf("[%s] unknown variable type!\n",__func__);

      break;
    } /* end of switch */

}
