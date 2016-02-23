/*

Nonlinear solid solver implemented in C and Petsc

There are three group of class (objects):
Solid(s), Interface(mk), Fluid(f)

Peggy Huang
1. Finished on Feb 03
2. Modularized on Feb 18

*/

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <petscsys.h>
#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscis.h>


#include "Field_s.h"
#include "user_param.h"
#include "init_setting.h"
#include "boun_func.h"
#include "init_cond.h"
#include "fluid.h"
#include "solidsolver.h"
#include "interface_marker.h"

static char help[] = "Nonlinear solid solver.\n\n";

int main(int argc,char **argv)
{
	Field_S 	solid;
	Lag_marker 	marker;
	Field_F		flow;
	AppCtx		user;

	int			step;
	double		timestep;

    // =============== Setting up PETSC =================
	MPI_Init(&argc, &argv);
    PetscInitialize(&argc,&argv,(char*)0,help);
    MPI_Comm_rank(MPI_COMM_WORLD,&(user.rank));
    MPI_Comm_size(MPI_COMM_WORLD,&(user.size));

    // =============== Initialization ====

    //fluid_setup(&flow, &solid);
    SolidInitialize(&solid, &user);
    //InterfaceInitialize(&marker,&solid, &flow, &user);

    // iteration
    //for(step = 1; step <= solid.time.nstep; step++)
    //{
    	timestep = step*solid.time.dt;

    	//InterfaceUpdate(&marker, &solid, &flow);	// update marker
    	// fluid solver

    	//InterfaceTraction(&marker, &solid, &flow); 	// force calculation
    	SolidSolver(&solid);
    	// check

    //}

    // finalize
    SolidFinalize(&solid);

	PetscFinalize();
	MPI_Finalize();
    return 0;
}

