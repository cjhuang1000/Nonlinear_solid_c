/*
 * fluid.h = field.h by SCH
 *
 *  Created on: Oct 25, 2015
 *      Author: peggyhuang
 */
//DONE
#ifndef FLUID_H_
#define FLUID_H_ 1

#include "petsc.h"
#include "petscdm.h"
#include "petscsys.h"
#include "petscvec.h"
#include "petscmat.h"

#include "mpi.h"

#include "Field_s.h"
#include "user_param.h"

typedef struct{

  int is_first_proc_x1;
  int is_first_proc_x2;
  int is_first_proc_x3;

  int is_last_proc_x1;
  int is_last_proc_x2;
  int is_last_proc_x3;

} End_proc;


/* index range for each process */

typedef struct {
  int start;
  int end;
} Int_bounds;

typedef struct {
  Int_bounds u_rhs, u_plot;
  Int_bounds v_rhs, v_plot;
  Int_bounds w_rhs, w_plot;
  Int_bounds cv;
} Idx_range;


/* index range */

typedef struct {

  /* index range for RHS's */
  int urhs[3][2], urhs_gbl[3][2],
  // @@ Jan 2016 add 'shared', means the boundary and the processor between
  // shared panels (which are the ghosted panel for neighboring processor)
      vrhs[3][2], vrhs_gbl[3][2],
      wrhs[3][2], wrhs_gbl[3][2],
      prhs[3][2], prhs_gbl[3][2];

  /* index range for unknowns */
  int u[3][2], u_gbl[3][2], u_raw[3][2],
      v[3][2], v_gbl[3][2], v_raw[3][2],
      w[3][2], w_gbl[3][2], w_raw[3][2],
      p[3][2], p_gbl[3][2], p_raw[3][2];

  int    shared_pro[6];
  double shared_bound[6], int_bound[6];
  // @@ Jan 2016 add 'shared', means the boundary and the processor between
  // shared panels (which are the ghosted panel for neighboring processor)
  // and the pure interior panels.
  // four entries mean the value of the left, right, bottom, top, front and
  // back sides of the subdomain

  /* index range for with ghost points */
  int u_ghost[3][2], u_raw_ghost[3][2];
  int v_ghost[3][2], v_raw_ghost[3][2];
  int w_ghost[3][2], w_raw_ghost[3][2];
  int p_ghost[3][2], p_raw_ghost[3][2];

  /* index range for grids */
  int x[2], y[2], z[2], xm[2], ym[2], zm[2];

} index_range_t;

/* benchmark timings */

typedef struct {

  double convect;
  double viscous;
  double pgrad;
  double copy_rhs;
  double divg;
  double poisson;
  double u_update;

} timing_t;




/* data structure to hold forcing points information:
   i[3]    : (i,j,k) of forcing point.
   tp[3]   : indices of interpolation template points, in terms of the
             indices of the offset2d/offset3d arrays.
   n[3]    : normal vector passing the forcing point.
   surf[3] : coordinates of the surface point.

   // Added by Peggy Oct.2015
   panel   : corresponding global index of boundary element
   ratio   :

*/
typedef struct {
  int                 i[3];
  int                tp[3];
  double              n[3];
  double           surf[3];
  //int				 panel; // added by Peggy Nov,1,2015
  //double			 ratio; // added by Peggy Nov,1,2015
  //double			 value; // boundary point property for traction calculation.
  	  	  	  	  	  	  	// It is velocities for velocity forcing points, and acceleration for pressure forcing points.
} forcing_point_t;

/* data  : pointer to data in the forcing point list.
   total : current total number elements in the list.
   max   : maximum number of elements in the list (for memory allocation).
*/
typedef struct {
  forcing_point_t  *data;
  int               total;
  int               max;
} forcing_point_list_t;


/* boundary condition type */

enum bc_type {PERIODIC=0, DIRICHLET, NEUMANN, TRACTION, CONVECTIVE};

typedef struct
{
  double x_value[2];
  double y_value[2];
  double z_value[2];

  char x_type_name[2][32];
  char y_type_name[2][32];
  char z_type_name[2][32];

  int x_type[2];
  int y_type[2];
  int z_type[2];

} bc_t;



/* Field data structure : this is the central part of the whole
   computation */

typedef struct {

  MPI_Comm comm;

  DM da;           /* DA for velocities */
//  DA da_s;         /* DA for CV-center scalars */

  DM da2d5p; // for poisson solvers
  DM da3d7p;

  int Nx, Ny, Nz;

  /* number of control volumes */
  int ncv[3];

  /* process-grid coordinates */
  int proc_coords[3];

  int x1_is_periodic, x2_is_periodic, x3_is_periodic;
  int verbose_mode;

  double cfl_number, reynolds_number, time_step_size;
  double cfl_max;
  double current_time;

  int time_step_index_start;
  int num_time_steps, export_every_n_steps;
  int have_immersed_boundary;
  int read_ic_from_file;

  double dx;
  double *x_grid, *y_grid, *z_grid;
  double *xm_grid, *ym_grid, *zm_grid;

  Vec u, u_rhs, u_rhs2;
  Vec v, v_rhs, v_rhs2;
  Vec w, w_rhs, w_rhs2;
  Vec phi, p, p_rhs;

  index_range_t idx_range;

  End_proc u_id_range; // v_id_range, w_id_range;

  //Idx_range index_range;

  Mat A_poisson;

  timing_t timing;

  /* mesh file names */
  char xmesh_file[80], ymesh_file[80], zmesh_file[80];
  char ic_type_name[80];

  /*--- immersed boundary --- */

  /* distance functions */
  Vec dist_u, dist_v, dist_w, dist_p;
  Vec dist_map_u, dist_map_v, dist_map_w, dist_map_p; // map of forcing points (if 1, forcing point. If 0, solid or fluid point)
  Vec panel_u, panel_v, panel_w, panel_p; // added by Peggy Nov,1,2015
  Vec ratio_u, ratio_v, ratio_w, ratio_p; // added by Peggy Nov,1,2015
  /* forcing point lists */
  forcing_point_list_t flist_u, flist_v, flist_w, flist_p;

  /* boundary conditions */

  bc_t u_bc, v_bc, w_bc;

  /* HYPRE stuff */
  //HYPRE_StructMatrix  A_hypre;
  //HYPRE_StructVector  x_hypre,b_hypre;
  //HYPRE_StructSolver  solver_hypre, precond_hypre;
  //HYPRE_StructGrid    grid_hypre;

} Field_F;

// -----------------------
// This section is added just for the testing of the
// solid solver and interface code in c


void fluid_setup(Field_F* f, Field_S* s);
void index_range_get_raw(char which_var, Field_F *f,
                        int *i, int *j, int *k);
void index_range_get_rhs(char which_var, Field_F *f,
                        int *i, int *j, int *k);
void index_range_get_raw_ghost(char which_var, Field_F *f,
                        int *i, int *j, int *k);
//------------------------

#endif /* FLUID_H_ */
