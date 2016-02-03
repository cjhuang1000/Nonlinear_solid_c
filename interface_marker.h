#ifndef INTERFACE_MARKER_H
#define INTERFACE_MARKER_H 1

#include "Field_s.h"
#include "fluid.h"
#include "nrutil.h"

#include <math.h>

#include <petsc.h>
#include <petscvec.h>
#include <petscis.h>
#include <petscdmda.h>
#include <petscdm.h>

// part of ib.h
// ---------------------

typedef struct {
  int              i;
  int              j;
} offset2d_t;

typedef struct {
  int              i;
  int              j;
  int              k;
} offset3d_t;

// ----------------------


// part of ib.c
// --------------------------
typedef struct {
    int i;
    double d;
} dist_t;

// pointer to the comparison function
typedef int (*compfn)(const void*, const void*);

static int comp_distance(dist_t *a, dist_t *b)
{
    if (a->d > b->d)
        return 1;
    else if (a->d < b->d)
        return -1;
    else
        return 0;
}
// --------------------------

struct interp2d{

	int		n; 			// number of stencils
	int		sten[4];
	double	coeff[4];
};

typedef struct {
    int i;
    int j;
    double d;
} dist_t2;

typedef struct {

	int 		nglobal, nlocal_s;              // local number of panels in the solid subdomain
	int			nlocal_f, nlocal_gho_f;			// local number of interior and ghosted panels in the fluid subdomain
	VecScatter 	glo2loc;						// scatter the distributed panels to the processor which has the panel										// in the fluid subdomain

	Vec 			sp_x0, sp_y0, ep_x0, ep_y0;
	Vec 			sp_x , sp_y , ep_x , ep_y ;
	Vec 			sp_u , sp_v , ep_u , ep_v ;
	Vec 			sp_ax, sp_ay, ep_ax, ep_ay;

	// compact marker point representation (be don ein the future)
	// distributed vec
	//Vec 			mkp_x0; 						// initial location of markers

	// local sequential vec
	//Vec			mkp_x;							// current velocity of all markers (sequential)
	//Vec			mkp_dx;							// current acceleration of all markers (sequential)
	//Vec			mkp_ddx; 						// current location of all markers

	// contains all the information
	Vec			l0;							// initial length
	Vec			tx, ty;						// indices of traction cell of all markers

	// have nonzero values if the marker is owned by local processor
	int			*list_f;
	double		*nx, *ny;	  					// nx, ny only have value on the panels in the current fluid processor

	struct interp2d *interp_s_x,*interp_s_y,*interp_e_x,*interp_e_y;

} Lag_marker;

// before the iteration
int forcing_pts_list_init(forcing_point_list_t *list);

void marker_setup (Grid_S* ptr_g, Solid* ptr_s, Index_S* ptr_i, Lag_marker* ptr_mk, AppCtx* ptr_u); // v
void list_setup(Lag_marker* ptr_mk, Field_F* ptr_f);
void disp_interp_setup (Lag_marker* ptr_mk, Index_S* ptr_i, Grid_S* ptr_g ); // v

void build_current_marker (Lag_marker* ptr_mk, Field_S* ptr_s, Field_F* ptr_f, Index_S* ptr_i, AppCtx* ptr_u);
void update_sdf (Grid_S* ptr_g, Lag_marker* ptr_mk, AppCtx* ptr_u, Field_F* ptr_f);
void force_calculation (Index_S* ptr_i, Field_F* ptr_f, Field_S* ptr_s, Lag_marker* ptr_mk);

void ib_find_forcing_pts_2d (Field_F* ptr_f, Lag_marker* ptr_mk);
void apply_forcing_2d (Field_F *f, const double sor);
void A_times_b_3x3 (double *a, double *b, double *c); // v


#endif  /* INTERFACE_MARKER_H */
