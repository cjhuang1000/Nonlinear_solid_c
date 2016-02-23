/*
 * interface.c
 *
 *  Created on: Feb 18, 2016
 *      Author: peggyhuang
 */


#include "interface.h"

// initialize, setting up
// marker initialization should be after the system reduction in
// compute_matricesNonlinearStructure(...)
void InterfaceInitialize(Lag_marker* mk, Field_S* s, Field_F* f, AppCtx* ptr_u){

    marker_setup(mk, s, ptr_u); // should be after reduce_system
    list_setup  (mk, f);
    disp_interp_setup(mk, s);
}


// Iteration, update the displacement
void InterfaceUpdate(Lag_marker* mk, Field_S* s, Field_F* f){

    build_current_marker(mk, s, f);
	update_sdf (mk,f);

}

void InterfaceTraction(Lag_marker* mk, Field_S* s, Field_F* f)
{
	ib_find_forcing_pts_2d (f, mk);
	apply_forcing_2d(f, 1.0);
	force_calculation(f, s, mk);
}

// output the field to day file
void InterfaceOutput(Lag_marker* mk){

	printf("output of interface");

}

void InterfaceFinalize(Lag_marker* mk){

	printf("finalize interface");

}
