#ifndef FIELD_S_H
#define FIELD_S_H 1

#include <petsc.h>
#include <petscsys.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscis.h>
#include <petscdmda.h>

typedef struct {

	PetscMPIInt	rank,size;

}AppCtx;

typedef struct
{

	short int *G2g,			*l2G;
	short int *G2g_before,	*l2g;
	short int *C2c_fsi,		*c2C_fsi;
	short int *C2c_neumann,	*c2C_neumann;
	short int *C2c_dirichlet,*c2C_dirichlet;
    short int *keep;

    // Gx,Gy: global index on structured grid
    // g: global index of the involved grid with new order
    //    x1->y1->x2->y2->... where the number denotes the
    //    processor where the node is on.
    // l: local index
} Index_xixy;

typedef struct {

  /* list of the internal/boundary cells indices*/
  unsigned int *cell_interior;
  unsigned int *cell_boundary;
  unsigned int *cell_fsi;
  unsigned int *cell_neumann;
  unsigned int *cell_dirichlet;

  IS is_xi, is_dir,is_neu;
  /* list of indices of involved variables*/
  Index_xixy xix;

  /* Position of the involved variables in the list*/
  Index_xixy xiy;

  /* Number of involved cells/variables (local)*/
  int cell_N_interior,cell_N_FSI,cell_N_boundary,cell_N_Neumann,cell_N_Dirichlet;
  int xix_N         ,xiy_N;
  int xix_N_before  ,xiy_N_before;
  int xix_Ncell_FSI ,xiy_Ncell_FSI;
  int xix_Ncell_Neumann     ,xiy_Ncell_Neumann;
  int xix_Ncell_Dirichlet   ,xiy_Ncell_Dirichlet;

  /* Number of involved cells/variables (global)*/
  int xi_gloN, xi_gloN_before;
  int xix_gloN, xix_gloN_before;
  int xix_ghoN, xiy_ghoN;
  int xi_gloNcell_Neumann     ,xi_gloNcell_Dirichlet;
  int xix_gloNcell_Neumann    ,xix_gloNcell_Dirichlet;

} Index_S;


/* boundary type at the solid boundary */
typedef struct {

	// body force
	double *body_funct_x,  *body_funct_y;

	// FSI Surface Forces on the staggerred grid
	double *fsi_func_x, *fsi_func_y;

	// Neumann Suface Forces Definition
    double *neumann_funct_x, *neumann_funct_y;

	// Dirichlet Suface Constraints Definition
	double *dirichlet_dxi_x,*dirichlet_dxi_y;

	// Projection of the surface specification on the main grid
	char *fsineumanndirichlet;

	Vec xi_body, xi_Neumann, xi_Dirichlet;
} Constraint_S;

typedef struct {

	double  length;
	double	rho;
	double	mu;
	double	lambda;

	char	type[80];
	char	initial[80];
	char	boundcond;

	int		**boundary_sign;
	double	**boundary_value;

	double  threshold;
	double  damping[2]; // damping alpha and beta

}Param_S;


typedef struct {

	double	initial;
	double	dt;
	int		nstep;
	int		nstep_output;
	double	delta;
	double	theta;

}TimeMarching_S;

typedef struct {

    // ----  displacement,velocity, acceleration [xi.x1; xi.y1; xi_x2;xi_y2;...] ---
	Vec xi,	xi_old, xi_k;
    Vec dxi, dxi_old, dxi_k;
    Vec ddxi, ddxi_old, ddxi_k; 			// at current step, previous step, current iteration
    Vec inc_dxi;							// increment analysis
	Vec traction; 							// traction on the solid boundary

	// The index of the solid staggered grid is defined in the following manner
	//  +------------+
	//  |            |
	//  |            |
	// ---   i,j     |
	//  |            |
	//  |            |
	//  +------|-----+
	//

	Index_S ind;

	// ----- Governing matrices: mass, stiffness matrices, load vectors
	Mat MS; // constant part of the LHS
		// MS.xx  |
		// -------+-------
		//        | MS.yy
	Mat DS;
		//  DS.xc |
		// -------+-------
		//        | DS.yc
	Mat KLS,KNS;
		// KLS.xx | KLS.xy
		// -------+--------
		// KLS.xy'| KLS.yy
	Mat NS;
		// NS.xc  |
		// -------+--------
		//        | NS.yc
    Vec FS;

    Constraint_S 	con;
    Param_S			param;
    TimeMarching_S	time;

	PetscInt    *v2p;  // the index of the processor each cell belongs to
	VecScatter  scatter_x, scatter_y;
	int 		dmx_s, dmy_s, dmx_e, dmy_e; // subdomain boundary

	// global domain dimension
	double		char_lengh;  // characteristic length of the problem
	double		domain_x[2],domain_y[2]; // domain interval
	double		domain_size[2];
	double		res,dx;   // resolution and grid spacing
	int			Nx,Ny,N;     // number of grid in x- and y- direction, total grid number
	double		*x_grid, *y_grid;

	// for solving Ax = RHS
	Mat			A,subA[4];
	Vec			RHS,subRHS[2],x,subx[2], tempvec, tempvec2;
	IS 			isg[2];
	KSP			ksp;
	PC			pc;

}Field_S;


#endif  /* FIELD_S_H */
