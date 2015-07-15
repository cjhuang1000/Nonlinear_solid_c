#ifndef BOUN_FUNC_H
#define BOUN_FUNC_H 1

#include "Field_s.h"
#include "init_setting.h"

#include <petsc.h>

typedef struct {

	double KS_hxxhxx[6][6];
	double KS_hyyhyy[6][6];
	double KS_hxyhxy[6][6];
	double KS_hyxhyx[6][6];
	double KS_hxxhyx[6][6];
	double KS_hxxhyy[6][6];
	double KS_hxyhyx[6][6];
	double KS_hxyhyy[6][6];
	double KS_hxxhxy[6][6];
	double KS_hyyhyx[6][6];
	double FS_hxx[6];
	double FS_hxy[6];
	double FS_hyy[6];
	double FS_hyx[6];

} boundfunc_cutcell;

typedef struct {

	/* % Internal Cell Elemental Matrices*/
	// integration within the internal grid
	double MS_xx[6][6];
	double MS_yy[6][6];
	double KS_hxxhxx[6][6];
	double KS_hyyhyy[6][6];
	double KS_hxyhxy[6][6];
	double KS_hyxhyx[6][6];
	double KS_hxxhyx[6][6];
	double KS_hxxhyy[6][6];
	double KS_hxyhyx[6][6];
	double KS_hxyhyy[6][6];
	double KS_hxxhxy[6][6];
	double KS_hyyhyx[6][6];
	double FS_hxx[6];
	double FS_hxy[6];
	double FS_hyy[6];
	double FS_hyx[6];
	double NS_xc[6][2];
	double NS_yc[6][2];
	double DS_xc[6][2];
	double DS_yc[6][2];

} boundfunc_cell;

typedef struct {

	/*Symetry consideration*/
	short int cellx;
	short int celly;
	short int signx;
	short int signy;
	short int x[6];
	short int y[6];
	short int subcell[4];

} bound_func_sym;

typedef struct {

	/* Elementary functions coefficients */
	/* % Coefficients which characterize the bilinear elementary functions N^xi
	 involved on the subcell omega1 (N_x1, N_x2, N_x3, N_x4, N_y7, N_y8,
	 N_y10, N_y11), such that, for (x,y) in omega1:
                       N(x,y) = eps*(x+a)*(y+b)
	The array below contains the [eps; a; b] coefficients for each elementary
	 functions N*/
	double coeff[3][8];
    bound_func_sym sym[4];

}boundary_function;

struct shapefunction{

    int               **info;
    boundfunc_cutcell *cutcell;

};
/*Internal Cell Elemental Matrices*/

void set_boundfunc(double dx);
void compute_matricesNonlinearStructure(Matrices_S* ptr_ms, Index_S* ptr_i, Grid_S* ptr_g, Solid* ptr_s, char* fnd);

#endif  /* BOUN_FUNC_H */
