/*
Nonlinear solid solver implemented in C and Petsc

Peggy Huang
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
#include "interface_marker.h"
#include "fluid.h"

static char help[] = "Nonlinear solid solver.\n\n";

// Copy the value in field_s f1 to f2
void fieldCopy (Field_S* ptr_f1, Field_S* ptr_f2)
{
	VecCopy(ptr_f1->xi,ptr_f2->xi);
	VecCopy(ptr_f1->dxi,ptr_f2->dxi);
	VecCopy(ptr_f1->ddxi,ptr_f2->ddxi);
	VecCopy(ptr_f1->inc_dxi,ptr_f2->inc_dxi);
}

void setScattering(AppCtx* ptr_u, Index_S* ptr_i,Vec f)
{

	Vec				local_x,local_y;
    IS 				is_xix,is_xiy;
    PetscInt		*l2g_x,*l2g_y;
    PetscInt		nx = ptr_i->xix_N + ptr_i->xix_ghoN, ny = ptr_i->xiy_N + ptr_i->xiy_ghoN;
    int				i;

    PetscMalloc1(nx,&l2g_x);
    PetscMalloc1(ny,&l2g_y);

    for(i = 0; i<nx; i++)
    	l2g_x[i] = ptr_i->xix.l2g[i];

    for(i = 0; i<ny; i++)
        l2g_y[i] = ptr_i->xiy.l2g[i];

	VecCreateSeq(PETSC_COMM_SELF,nx,&local_x);
    VecCreateSeq(PETSC_COMM_SELF,ny,&local_y);

    ISCreateGeneral(PETSC_COMM_SELF,nx,l2g_x,PETSC_OWN_POINTER,&is_xix);
    ISCreateGeneral(PETSC_COMM_SELF,ny,l2g_y,PETSC_OWN_POINTER,&is_xiy);

    VecScatterCreate(f,is_xix,local_x,NULL,&ptr_u->scatter_x);
    VecScatterCreate(f,is_xiy,local_y,NULL,&ptr_u->scatter_y);

    VecDestroy(&local_x); VecDestroy(&local_y);
    ISDestroy(&is_xix);  ISDestroy(&is_xiy);
}

int main(int argc,char **argv)
{

    int     step;
    int     i,j,n;
    
    double  timestep;
    double  residue0, residue1;

    Grid_S	    	grid;
    Solid	    	solid;
    TimeMarching    timem;
    Index_S	    	ind;
    Constraint_S    const_s;
    Field_S         field_s,field_s_old,field_s_k;
    Field_F			fluid;
    Matrices_S      mat_s;
    AppCtx          user;
    Lag_marker      marker0,marker;  	// marker0 represent

    Mat				A,subA[4];
    Vec				RHS,subRHS[2], tempvec, tempvec2,x,subx[2];
    IS 				isg[2];
    KSP				ksp;
    PC				pc;
    PetscErrorCode  ierr;
    PetscScalar  	c,c1;

    PetscInt        mm,nn;
    PetscScalar 	shift;

    // =============== Setting up PETSC =================

    ierr = PetscInitialize(&argc,&argv,(char*)0,help);CHKERRQ(ierr);
    MPI_Comm_rank(MPI_COMM_WORLD,&user.rank);
    MPI_Comm_size(MPI_COMM_WORLD,&user.size);

    // =============== Parameter setting ================

    // set simulation parameters
    user_param(&grid,&solid,&timem,&const_s,&user);

    // set fluid parameters
    fluid_setup(&fluid, &grid, &solid);

    // set indices
    set_index(&ind,&grid,&user,solid.boundary_sign,const_s.fsineumanndirichlet);

    // initialize boundary function
    set_boundfunc(grid.dx);

    // =============== Constructing governing matrices and vectors ================

    // construct the governing matrices
    compute_matricesNonlinearStructure(&mat_s, &ind, &grid, &solid, const_s.fsineumanndirichlet);

    n = ind.xix_N + ind.xiy_N;
    VecCreate(PETSC_COMM_WORLD,&mat_s.FS);
    VecSetSizes(mat_s.FS,n,ind.xi_gloN);
    VecSetFromOptions(mat_s.FS);

    //================== Setting constraints ======================

    VecDuplicate(mat_s.FS, &const_s.xi_body);
    VecDuplicate(mat_s.FS, &const_s.xi_Neumann);

    VecCreate(PETSC_COMM_WORLD,&const_s.xi_Dirichlet);
    VecSetSizes(const_s.xi_Dirichlet,ind.xix_Ncell_Dirichlet +ind.xiy_Ncell_Dirichlet,ind.xi_gloNcell_Dirichlet);
    VecSetFromOptions(const_s.xi_Dirichlet);

    VecCreate(PETSC_COMM_WORLD,&tempvec);
    VecSetSizes(tempvec,ind.xix_Ncell_Neumann +ind.xiy_Ncell_Neumann,ind.xi_gloNcell_Neumann);
    VecSetFromOptions(tempvec);

    for(i=0; i<ind.xix_N; i++)
    	VecSetValue(const_s.xi_body,ind.xix.l2g[i],const_s.body_funct_x[ind.xix.l2G[i]], INSERT_VALUES);

    for(i=0; i<ind.xiy_N; i++)
        VecSetValue(const_s.xi_body,ind.xiy.l2g[i],const_s.body_funct_y[ind.xiy.l2G[i]], INSERT_VALUES);

    for(i=0; i<ind.xix_Ncell_Dirichlet; i++)
    	VecSetValue(const_s.xi_Dirichlet,ind.xix.C2c_dirichlet[ind.xix.c2C_dirichlet[i]],const_s.dirichlet_dxi_x[ind.xix.c2C_dirichlet[i]], INSERT_VALUES);

    for(i=0; i<ind.xiy_Ncell_Dirichlet; i++)
        VecSetValue(const_s.xi_Dirichlet,ind.xiy.C2c_dirichlet[ind.xiy.c2C_dirichlet[i]],const_s.dirichlet_dxi_y[ind.xiy.c2C_dirichlet[i]], INSERT_VALUES);

    for(i=0; i<ind.xix_Ncell_Neumann; i++)
        VecSetValue(tempvec,ind.xix.C2c_neumann[ind.xix.c2C_neumann[i]],const_s.neumann_funct_x[ind.xix.c2C_neumann[i]], INSERT_VALUES);

    for(i=0; i<ind.xiy_Ncell_Neumann; i++)
        VecSetValue(tempvec,ind.xiy.C2c_neumann[ind.xiy.c2C_neumann[i]],const_s.neumann_funct_y[ind.xiy.c2C_neumann[i]], INSERT_VALUES);

    VecAssemblyBegin(const_s.xi_body);
    VecAssemblyEnd(const_s.xi_body);
    VecAssemblyBegin(const_s.xi_Dirichlet);
    VecAssemblyEnd(const_s.xi_Dirichlet);
    VecAssemblyBegin(tempvec);
    VecAssemblyEnd(tempvec);

    MatMult(mat_s.NS,tempvec,const_s.xi_Neumann);
    VecDestroy(&tempvec);

    // =============== Initial setting of field_s ================

    set_initial(&ind, &solid, &field_s);

    // initialize the field_s_old
    VecDuplicate(field_s.xi,&field_s_old.xi);
    VecDuplicate(field_s.xi,&field_s_old.dxi);
    VecDuplicate(field_s.xi,&field_s_old.ddxi);
    VecDuplicate(field_s.xi,&field_s_old.inc_dxi);

    // initialize the field_s_k
    VecDuplicate(field_s.xi,&field_s_k.xi);
    VecDuplicate(field_s.xi,&field_s_k.dxi);
    VecDuplicate(field_s.xi,&field_s_k.ddxi);
    VecDuplicate(field_s.xi,&field_s_k.inc_dxi);

    fieldCopy(&field_s,&field_s_old);

    // set up interface markers

    marker_setup(&grid, &solid, &ind, &marker, &user); // should be after reduce_system
    list_setup  (&marker, &fluid);
    disp_interp_setup(&marker, &ind, &grid );

    // ================ Setting Scattering ==========================

    setScattering(&user, &ind, field_s.xi);

    build_current_marker(&marker, &field_s, &fluid, &ind, &user);

	update_sdf (&grid, &marker, &user, &fluid);

	ib_find_forcing_pts_2d (&fluid, &marker);

	apply_forcing_2d(&fluid, 1.0);

	printf("works!");
	force_calculation(&ind, &fluid, &field_s, &marker);

	/*
    // ================ Time marching ===============================

    // Setup vec and mat
    VecDuplicate(field_s.xi, &subRHS[0]);
    VecDuplicate(const_s.xi_Dirichlet, &subRHS[1]);
    VecDuplicate(field_s.xi, &subx[0]);
    VecDuplicate(const_s.xi_Dirichlet, &subx[1]);
    VecDuplicate(field_s.xi,&tempvec);

    VecCreateNest(PETSC_COMM_WORLD,2,NULL,subx,&x);

    n = ind.xix_Ncell_Dirichlet +ind.xiy_Ncell_Dirichlet;
	MatCreate(PETSC_COMM_WORLD,&subA[3]);
	MatSetSizes(subA[3],n,n,ind.xi_gloNcell_Dirichlet,ind.xi_gloNcell_Dirichlet);
	MatSetFromOptions(subA[3]);
	MatSetUp(subA[3]);
	MatZeroEntries(subA[3]);
	MatAssemblyBegin(subA[3],MAT_FINAL_ASSEMBLY);
	MatAssemblyEnd(subA[3],MAT_FINAL_ASSEMBLY);

	MatDuplicate(mat_s.DS,MAT_COPY_VALUES,&subA[1]); // DS
	MatTranspose(mat_s.DS,MAT_INITIAL_MATRIX,&subA[2]);			// DS^T
	MatDuplicate(mat_s.MS,MAT_DO_NOT_COPY_VALUES,&subA[0]);

	// ----------------------------

	for(step = 1; step <= timem.nstep; step++)
    {
    	// first iteration

		timestep = step*timem.dt;
	    compute_matricesNonlinearStructure_update(&mat_s, &ind, &grid, &solid, &user, &field_s_old);

	    // Governing matrix A
	    MatCopy(mat_s.MS,subA[0],DIFFERENT_NONZERO_PATTERN);

	    c = solid.damping[0] + 1/timem.dt/timem.delta;
	    MatScale(subA[0],c);
	    c = solid.damping[1] + timem.dt * timem.theta / timem.delta;
	    MatAXPY(subA[0],c,mat_s.KLS,DIFFERENT_NONZERO_PATTERN);
	    c = timem.dt * timem.theta / timem.delta;
	    MatAXPY(subA[0],c,mat_s.KNS,DIFFERENT_NONZERO_PATTERN);

	    MatCreateNest(PETSC_COMM_WORLD,2,NULL,2,NULL,subA,&A); // !!!!
	    MatNestGetISs(A,isg,NULL);

	    // RHS
	    c = 1.0;
	    VecWAXPY(subRHS[0],c,const_s.xi_Neumann,const_s.xi_body);
	    c = -solid.damping[0]; c1 = 1/timem.delta -1;
	    VecCopy(field_s_old.ddxi,tempvec);
	    VecAXPBY(tempvec,c,c1,field_s_old.dxi);
	    MatMultAdd(mat_s.MS,tempvec,subRHS[0],subRHS[0]);

	    c = -timem.dt - solid.damping[1]; c1 = -timem.dt*timem.dt* (1-2*timem.theta/timem.delta);
	    VecCopy(field_s_old.ddxi,tempvec);
	    VecAXPBY(tempvec,c,c1,field_s_old.dxi);
	    MatMultAdd(mat_s.KLS,tempvec,subRHS[0],subRHS[0]);

	    c = -timem.dt; c1 = -timem.dt*timem.dt* (1-2*timem.theta/timem.delta);
	    VecCopy(field_s_old.ddxi,tempvec);
	    VecAXPBY(tempvec,c,c1,field_s_old.dxi);
	   	MatMultAdd(mat_s.KNS,tempvec,subRHS[0],subRHS[0]);
	   	c = -1.0;
	    ierr = VecAXPY(subRHS[0],c,mat_s.FS); CHKERRQ(ierr);

	   	VecCopy(field_s_old.dxi,tempvec);
	   	VecScale(tempvec,c);
	   	MatMult(subA[2],tempvec,subRHS[1]);

	   	VecCreateNest(PETSC_COMM_WORLD,2,NULL,subRHS,&RHS);

		// seting the KSP solver
	    KSPCreate(PETSC_COMM_WORLD,&ksp);
	    KSPSetType(ksp, KSPGMRES);
		KSPSetOperators(ksp,A,A);
		//KSPSetTolerances(ksp,PETSC_DEFAULT,1.e-20,PETSC_DEFAULT,PETSC_DEFAULT);
		KSPSetFromOptions(ksp);
		KSPGetPC(ksp,&pc);
		PCFieldSplitSetIS(pc,"0",isg[0]);
		PCFieldSplitSetIS(pc,"1",isg[1]);

		// solve
		ierr = KSPSolve(ksp,RHS,x); CHKERRQ(ierr);

		// check ###############
		//Vec checkvec;
		//VecDuplicate(RHS,&checkvec);
		//c = -1.0;
		//MatMult(A,x,checkvec);
		//VecAXPY(checkvec,c,RHS);
		//VecNestGetSubVec(checkvec,0,&tempvec2);
		//VecView(tempvec2,PETSC_VIEWER_STDOUT_WORLD);
		// #####################

		VecNestGetSubVec(x,0,&tempvec2);
		VecCopy(tempvec2,field_s.inc_dxi);

		//update
		c = (timem.delta - 1.0)*timem.dt;  c1 = 1/timem.delta/timem.dt;
		VecWAXPY(field_s.ddxi,c,field_s_old.ddxi,field_s.inc_dxi);
		VecScale(field_s.ddxi,c1);

		c = timem.dt;
		VecWAXPY(field_s.xi,c,field_s_old.dxi,field_s_old.xi);

		c = timem.theta*timem.dt*timem.dt;  c1 = (1-2*timem.theta)*timem.dt*timem.dt/2;
		VecAXPY(field_s.xi,c,field_s.ddxi);
		VecAXPY(field_s.xi,c1,field_s_old.ddxi);

		c = 1.0;
		VecWAXPY(field_s.dxi,c,field_s_old.dxi,field_s.inc_dxi);
		// residue
		VecNorm(field_s.xi,NORM_2,&residue0);
		residue1 = residue0;

		fieldCopy(&field_s,&field_s_k);
		//VecView(field_s.xi,PETSC_VIEWER_STDOUT_WORLD);

		// while loop
		while ( (residue1/residue0) > solid.threshold)
		{
			//printf("looping\n");
			compute_matricesNonlinearStructure_update(&mat_s, &ind, &grid, &solid, &user, &field_s_k);

			// Governing matrix A
			MatCopy(mat_s.MS,subA[0],DIFFERENT_NONZERO_PATTERN); CHKERRQ(ierr);
			c = solid.damping[0] + 1/timem.dt/timem.delta;
			MatScale(subA[0],c);
			c = solid.damping[1] + timem.dt * timem.theta / timem.delta;
			MatAXPY(subA[0],c,mat_s.KLS,DIFFERENT_NONZERO_PATTERN);
			c = timem.dt * timem.theta / timem.delta;
			MatAXPY(subA[0],c,mat_s.KNS,DIFFERENT_NONZERO_PATTERN);

			MatCreateNest(PETSC_COMM_WORLD,2,NULL,2,NULL,subA,&A);

			// RHS
			c = 1.0;
			VecWAXPY(subRHS[0],c,const_s.xi_Neumann,const_s.xi_body);

			c = -1.0; c1 = 1.0/timem.dt/timem.delta;
			VecWAXPY(tempvec,c,field_s_k.dxi,field_s_old.dxi);
			VecScale(tempvec,c1);

			c = -solid.damping[0]; c1 = 1.0/timem.delta -1.0;
			VecAXPY(tempvec,c,field_s_k.dxi);
			VecAXPY(tempvec,c1,field_s_old.ddxi);
			MatMultAdd(mat_s.MS,tempvec,subRHS[0],subRHS[0]); // ok to here

			c = -1.0; c1 = -timem.dt*timem.theta/timem.delta - solid.damping[1];
			VecWAXPY(tempvec,c,field_s_old.xi,field_s_k.xi);
			VecAXPY(tempvec,c1,field_s_k.dxi);
			c = -timem.dt*(1.0-timem.theta/timem.delta); c1 = -timem.dt*timem.dt* (1-2*timem.theta/timem.delta);
			VecAXPY(tempvec,c1,field_s_old.dxi);
			VecAXPY(tempvec,c,field_s_old.ddxi);
			MatMultAdd(mat_s.KLS,tempvec,subRHS[0],subRHS[0]);

			c = -1.0; c1 = -timem.dt*timem.theta/timem.delta;
			VecWAXPY(tempvec,c,field_s_old.xi,field_s_k.xi);
			VecAXPY(tempvec,c1,field_s_k.dxi);
			c = -timem.dt*(1.0-timem.theta/timem.delta); c1 = -timem.dt*timem.dt* (1-2*timem.theta/timem.delta);
			VecAXPY(tempvec,c1,field_s_old.dxi);
			VecAXPY(tempvec,c,field_s_old.ddxi);
			MatMultAdd(mat_s.KNS,tempvec,subRHS[0],subRHS[0]);

			c = -1.0;
			VecAXPY(subRHS[0],c,mat_s.FS);
			//VecView(subRHS[0],PETSC_VIEWER_STDOUT_WORLD);

			VecCopy(field_s_old.dxi,tempvec);
			VecScale(tempvec,c);
			MatMult(subA[2],tempvec,subRHS[1]);

		   	VecCreateNest(PETSC_COMM_WORLD,2,NULL,subRHS,&RHS);

			// seting the KSP solver
			KSPDestroy(&ksp);
			KSPCreate(PETSC_COMM_WORLD,&ksp);
			KSPSetType(ksp, KSPFGMRES);
			KSPSetOperators(ksp,A,A);
			KSPSetFromOptions(ksp);
			KSPGetPC(ksp,&pc);
			PCFieldSplitSetIS(pc,"0",isg[0]);CHKERRQ(ierr);
			PCFieldSplitSetIS(pc,"1",isg[1]);CHKERRQ(ierr);

			// solve
			KSPSolve(ksp,RHS,x);
			VecNestGetSubVec(x,0,&tempvec2);
			VecCopy(tempvec2,field_s.inc_dxi);

			//update
			c = 1.0; c1  = -1.0;
			VecWAXPY(field_s.dxi,c,field_s_k.dxi,field_s.inc_dxi);
			VecWAXPY(tempvec,c1,field_s_old.dxi,field_s.dxi);

			c = (timem.delta - 1.0)*timem.dt;  c1 = 1/timem.delta/timem.dt;
			VecWAXPY(field_s.ddxi,c,field_s_old.ddxi,tempvec);
			VecScale(field_s.ddxi,c1);

			c = timem.dt;
			VecWAXPY(field_s.xi,c,field_s_old.dxi,field_s_old.xi);
			c = timem.theta*timem.dt*timem.dt;  c1 = (1-2*timem.theta)*timem.dt*timem.dt/2.0;
			VecAXPY(field_s.xi,c,field_s.ddxi);
			VecAXPY(field_s.xi,c1,field_s_old.ddxi);

			// residue
			c = -1.0;
			VecWAXPY(tempvec,c,field_s_k.xi,field_s.xi);
			VecNorm(tempvec,NORM_2,&residue1);

			fieldCopy(&field_s,&field_s_k);

			if (user.rank == 0)
				printf("residue1 ratio: %e %e %e\n", residue1,residue0,residue1/residue0);
			//residue1 = residue0*0.00001;

	}
	// output
	//if(step%timem.nstep_output == 0)
	//MatView(mat_s.KLS,PETSC_VIEWER_STDOUT_WORLD); // --------------
	//PetscObjectSetName((PetscObject)field_s.xi,"xi");
	//VecView(field_s.xi,PETSC_VIEWER_MATLAB_WORLD);
	VecView(field_s.xi,PETSC_VIEWER_STDOUT_WORLD); // --------------
	fieldCopy(&field_s,&field_s_old);
    }
	*/

	PetscFinalize();
    return 0;
}


