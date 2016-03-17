/*
 * solidsolver.c
 *
 *  Created on: Feb 17, 2016
 *      Author: peggyhuang
 */

#include "solidsolver.h"

void fieldCopy_current_to_old (Field_S* s);
void fieldCopy_current_to_k   (Field_S* s);
void setScattering			  (Field_S* s);

// initialize, setting up
void SolidInitialize(Field_S* s, AppCtx* ptr_u){

	int			n, i;

	user_param(s, ptr_u); 	// set simulation parameters
    set_index(s, ptr_u); 	// set indices
    set_boundfunc(s->dx); 	// initialize boundary function



    // =============== Constructing governing matrices and vectors ================

    compute_matricesNonlinearStructure(s);	// construct the governing matrices

    n = s->ind.xix_N + s->ind.xiy_N;
    VecCreate(PETSC_COMM_WORLD,&s->FS);
    VecSetSizes(s->FS,n,s->ind.xi_gloN);
    VecSetFromOptions(s->FS);

    //================== Setting constraints ======================

    printf("[Setting constraints] started\n");

     VecDuplicate(s->FS, &s->con.xi_body);
     VecDuplicate(s->FS, &s->con.xi_Neumann);

     VecCreate(PETSC_COMM_WORLD,&s->con.xi_Dirichlet);
     VecSetSizes(s->con.xi_Dirichlet,s->ind.xix_Ncell_Dirichlet +s->ind.xiy_Ncell_Dirichlet,s->ind.xi_gloNcell_Dirichlet);
     VecSetFromOptions(s->con.xi_Dirichlet);

     VecCreate(PETSC_COMM_WORLD,&(s->tempvec));
     VecSetSizes(s->tempvec,s->ind.xix_Ncell_Neumann +s->ind.xiy_Ncell_Neumann,s->ind.xi_gloNcell_Neumann);
     VecSetFromOptions(s->tempvec);

     for(i=0; i<s->ind.xix_N; i++)
     	VecSetValue(s->con.xi_body,s->ind.xix.l2g[i],s->con.body_funct_x[s->ind.xix.l2G[i]], INSERT_VALUES);

     for(i=0; i<s->ind.xiy_N; i++)
        VecSetValue(s->con.xi_body,s->ind.xiy.l2g[i],s->con.body_funct_y[s->ind.xiy.l2G[i]], INSERT_VALUES);

     for(i=0; i<s->ind.xix_Ncell_Dirichlet; i++)
     	VecSetValue(s->con.xi_Dirichlet,s->ind.xix.C2c_dirichlet[s->ind.xix.c2C_dirichlet[i]],s->con.dirichlet_dxi_x[s->ind.xix.c2C_dirichlet[i]], INSERT_VALUES);

     for(i=0; i<s->ind.xiy_Ncell_Dirichlet; i++)
        VecSetValue(s->con.xi_Dirichlet,s->ind.xiy.C2c_dirichlet[s->ind.xiy.c2C_dirichlet[i]],s->con.dirichlet_dxi_y[s->ind.xiy.c2C_dirichlet[i]], INSERT_VALUES);

     for(i=0; i<s->ind.xix_Ncell_Neumann; i++)
        VecSetValue(s->tempvec,s->ind.xix.C2c_neumann[s->ind.xix.c2C_neumann[i]],s->con.neumann_funct_x[s->ind.xix.c2C_neumann[i]], INSERT_VALUES);

     for(i=0; i<s->ind.xiy_Ncell_Neumann; i++)
        VecSetValue(s->tempvec,s->ind.xiy.C2c_neumann[s->ind.xiy.c2C_neumann[i]],s->con.neumann_funct_y[s->ind.xiy.c2C_neumann[i]], INSERT_VALUES);

     VecAssemblyBegin(s->con.xi_body);
     VecAssemblyEnd(s->con.xi_body);
     VecAssemblyBegin(s->con.xi_Dirichlet);
     VecAssemblyEnd(s->con.xi_Dirichlet);
     VecAssemblyBegin(s->tempvec);
     VecAssemblyEnd(s->tempvec);

     MatMult(s->NS,s->tempvec,s->con.xi_Neumann);
     VecDestroy(&(s->tempvec));

     printf("[Setting constraints] ended\n");

     // =============== Initial setting of the solid field ================

     printf("[Initial setting] started\n");
     set_initial(s);	// set initial

     VecDuplicate(s->xi,&(s->xi_old));
     VecDuplicate(s->xi,&(s->dxi_old));
     VecDuplicate(s->xi,&(s->ddxi_old));  // initialize the field_s_old

     VecDuplicate(s->xi,&(s->xi_k));
     VecDuplicate(s->xi,&(s->dxi_k));
     VecDuplicate(s->xi,&(s->ddxi_k)); // initialize the field_s_k

     fieldCopy_current_to_old(s);

     setScattering(s);

     printf("[Initial setting] ended\n");
     // ============== Preparation for solving the linear system =============

     printf("[Preparation for solving the linear system ] started\n");
     // Setup vec and mat
     VecDuplicate(s->xi, &(s->subRHS[0]));
     VecDuplicate(s->con.xi_Dirichlet, &(s->subRHS[1]));
     VecDuplicate(s->xi, &(s->subx[0]));
     VecDuplicate(s->con.xi_Dirichlet, &(s->subx[1]));
     VecDuplicate(s->xi, &(s->tempvec));

     VecCreateNest(PETSC_COMM_WORLD,2,NULL,s->subx,&(s->x));

     n = s->ind.xix_Ncell_Dirichlet +s->ind.xiy_Ncell_Dirichlet;
 	 MatCreate(PETSC_COMM_WORLD,&s->subA[3]);
 	 MatSetSizes(s->subA[3],n,n,s->ind.xi_gloNcell_Dirichlet,s->ind.xi_gloNcell_Dirichlet);
 	 MatSetFromOptions(s->subA[3]);
 	 MatSetUp(s->subA[3]);
 	 MatZeroEntries(s->subA[3]);
 	 MatAssemblyBegin(s->subA[3],MAT_FINAL_ASSEMBLY);
 	 MatAssemblyEnd(s->subA[3],MAT_FINAL_ASSEMBLY);

 	 MatDuplicate(s->DS,MAT_COPY_VALUES,&s->subA[1]); // DS
 	 MatTranspose(s->DS,MAT_INITIAL_MATRIX,&s->subA[2]);			// DS^T
 	 MatDuplicate(s->MS,MAT_DO_NOT_COPY_VALUES,&s->subA[0]);

 	printf("[Preparation for solving the linear system ] ended\n");
}


// Iteration, update the displacement
void SolidSolver(Field_S* s){

    int     i,j,n, loop=0;
    double  residue0, residue1;

    PetscScalar  	c,c1;
    PetscInt        mm,nn;
    TimeMarching_S  *timem = &(s->time);
    Param_S			*param = &(s->param);
    PetscErrorCode  ierr;
    PetscReal 		norm; // residual norm for ksp

    compute_matricesNonlinearStructure_update(s);

    // Governing matrix A
    MatCopy(s->MS,s->subA[0],DIFFERENT_NONZERO_PATTERN);

    c = param->damping[0] + 1.0/timem->dt/timem->delta;
    MatScale(s->subA[0],c);
    c = param->damping[1] + timem->dt * timem->theta / timem->delta;
    MatAXPY(s->subA[0],c,s->KLS,DIFFERENT_NONZERO_PATTERN);
    c = timem->dt * timem->theta / timem->delta;
    MatAXPY(s->subA[0],c,s->KNS,DIFFERENT_NONZERO_PATTERN);

    MatCreateNest(PETSC_COMM_WORLD,2,NULL,2,NULL,s->subA,&(s->A));
    MatNestGetISs(s->A,s->isg,NULL);

    // RHS
    c = 1.0;
    VecWAXPY(s->subRHS[0],c,s->con.xi_Neumann,s->con.xi_body);
    c = -param->damping[0]; c1 = 1.0/timem->delta -1.0;
    VecCopy(s->ddxi_old,s->tempvec);
    VecAXPBY(s->tempvec,c,c1,s->dxi_old);
    MatMultAdd(s->MS,s->tempvec,s->subRHS[0],s->subRHS[0]);

    c = -timem->dt - param->damping[1]; c1 = -timem->dt*timem->dt* (1.0-2.0*timem->theta/timem->delta);
    VecCopy(s->ddxi_old,s->tempvec);
    VecAXPBY(s->tempvec,c,c1,s->dxi_old);
    MatMultAdd(s->KLS,s->tempvec,s->subRHS[0],s->subRHS[0]);

    c = -timem->dt; c1 = -timem->dt*timem->dt* (1.0-2.0*timem->theta/timem->delta);
    VecCopy(s->ddxi_old,s->tempvec);
    VecAXPBY(s->tempvec,c,c1,s->dxi_old);
   	MatMultAdd(s->KNS,s->tempvec,s->subRHS[0],s->subRHS[0]);
   	c = -1.0;
    VecAXPY(s->subRHS[0],c,s->FS);

   	VecCopy(s->dxi_old,s->tempvec);
   	VecScale(s->tempvec,c);
   	MatMult(s->subA[2],s->tempvec,s->subRHS[1]);
   	VecAXPY(s->subRHS[1],1.0,s->con.xi_Dirichlet);

   	VecCreateNest(PETSC_COMM_WORLD,2,NULL,s->subRHS,&(s->RHS));

	//PetscViewer viewer;

    // write vec/mat
    //PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
    //PetscViewerSetType(viewer,PETSCVIEWERASCII );
    //PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);
    //PetscViewerFileSetName(viewer,"RHS.dat");
    //VecView(s->subRHS[0], viewer);
    //VecView(s->subRHS[1], viewer);
    //PetscViewerDestroy(&viewer);


	// seting the KSP solver
    KSPCreate(PETSC_COMM_WORLD,&(s->ksp));
    //KSPSetType(s->ksp, KSPGMRES);
    KSPSetInitialGuessNonzero(s->ksp,PETSC_TRUE);
	KSPSetOperators(s->ksp,s->A,s->A);
	//KSPSetTolerances(s->ksp,1e-5,PETSC_DEFAULT,PETSC_DEFAULT,PETSC_DEFAULT);
	KSPSetFromOptions(s->ksp);
	KSPGetPC(s->ksp,&(s->pc));
	PCFieldSplitSetIS(s->pc,"0",s->isg[0]);
	PCFieldSplitSetIS(s->pc,"1",s->isg[1]);

	printf("[KSPSolve] started\n");
	// solve
	KSPSolve(s->ksp,s->RHS,s->x);
	KSPGetResidualNorm(s->ksp,&norm);
	printf("[KSPSolve] ended, the residual norm is %e\n",norm);

	VecNestGetSubVec(s->x,0,&s->tempvec2);
	VecCopy(s->tempvec2,s->inc_dxi);

	//update
	c = (timem->delta - 1.0)*timem->dt;  c1 = 1.0/timem->delta/timem->dt;
	VecWAXPY(s->ddxi,c,s->ddxi_old,s->inc_dxi);
	VecScale(s->ddxi,c1);

	c = timem->dt;
	VecWAXPY(s->xi,c,s->dxi_old,s->xi_old);

	c = timem->theta*timem->dt*timem->dt;  c1 = (1.0-2.0*timem->theta)*timem->dt*timem->dt/2.0;
	VecAXPY(s->xi,c,s->ddxi);
	VecAXPY(s->xi,c1,s->ddxi_old);

	c = 1.0;
	VecWAXPY(s->dxi,c,s->dxi_old,s->inc_dxi);
	// residue
	VecNorm(s->xi,NORM_2,&residue0);
	printf("residual0 : %e\n", residue0);
	residue1 = residue0;

	fieldCopy_current_to_k (s);

	printf("[Enter While loop] started\n");
	// while loop
	while ( (residue1/residue0) > param->threshold)
	{
		loop++;
		printf("[While loop] looping\n");
		printf("%d %f\n",loop,(residue1/residue0));

		compute_matricesNonlinearStructure_update(s);

		// Governing matrix A
		MatCopy(s->MS,s->subA[0],DIFFERENT_NONZERO_PATTERN);
		c = param->damping[0] + 1.0/timem->dt/timem->delta;
		MatScale(s->subA[0],c);
		c = param->damping[1] + timem->dt * timem->theta / timem->delta;
		MatAXPY(s->subA[0],c,s->KLS,DIFFERENT_NONZERO_PATTERN);
		c = timem->dt * timem->theta / timem->delta;
		MatAXPY(s->subA[0],c,s->KNS,DIFFERENT_NONZERO_PATTERN);

		MatCreateNest(PETSC_COMM_WORLD,2,NULL,2,NULL,s->subA,&(s->A));

		// RHS
		c = 1.0;
		VecWAXPY(s->subRHS[0],c,s->con.xi_Neumann,s->con.xi_body);

		c = -1.0; c1 = 1.0/timem->dt/timem->delta;
		VecWAXPY(s->tempvec,c,s->dxi_k,s->dxi_old);
		VecScale(s->tempvec,c1);

		c = -param->damping[0]; c1 = 1.0/timem->delta -1.0;
		VecAXPY(s->tempvec,c, s->dxi_k);
		VecAXPY(s->tempvec,c1,s->ddxi_old);
		MatMultAdd(s->MS,s->tempvec,s->subRHS[0],s->subRHS[0]); // ok to here

		c = -1.0; c1 = -timem->dt*timem->theta/timem->delta - param->damping[1];
		VecWAXPY(s->tempvec,c,s->xi_old,s->xi_k);
		VecAXPY(s->tempvec,c1,s->dxi_k);
		c = -timem->dt*(1.0-timem->theta/timem->delta); c1 = -timem->dt*timem->dt* (1-2*timem->theta/timem->delta);
		VecAXPY(s->tempvec,c,s->dxi_old);
		VecAXPY(s->tempvec,c1, s->ddxi_old);
		MatMultAdd(s->KLS,s->tempvec,s->subRHS[0],s->subRHS[0]);

		c = -1.0; c1 = -timem->dt*timem->theta/timem->delta;
		VecWAXPY(s->tempvec,c,s->xi_old,s->xi_k);
		VecAXPY(s->tempvec,c1,s->dxi_k);
		c = -timem->dt*(1.0-timem->theta/timem->delta); c1 = -timem->dt*timem->dt* (1-2*timem->theta/timem->delta);
		VecAXPY(s->tempvec,c,s->dxi_old);
		VecAXPY(s->tempvec,c1,s->ddxi_old);
		MatMultAdd(s->KNS,s->tempvec,s->subRHS[0],s->subRHS[0]);

		c = -1.0;
		VecAXPY(s->subRHS[0],c,s->FS);

		VecCopy(s->dxi_k,s->tempvec);
		VecScale(s->tempvec,c);
		MatMult(s->subA[2],s->tempvec,s->subRHS[1]);
		VecAXPY(s->subRHS[1],1.0,s->con.xi_Dirichlet);
		//VecView(s->subRHS[1],PETSC_VIEWER_STDOUT_WORLD);

	   	VecCreateNest(PETSC_COMM_WORLD,2,NULL,s->subRHS,&(s->RHS));

		// seting the KSP solver
		KSPDestroy(&(s->ksp));
		KSPCreate(PETSC_COMM_WORLD,&(s->ksp));
		//KSPSetType(s->ksp, KSPFGMRES);
		KSPSetOperators(s->ksp,s->A,s->A);
		KSPSetInitialGuessNonzero(s->ksp,PETSC_TRUE); // nonzero initial guess
		KSPSetFromOptions(s->ksp);
		KSPGetPC(s->ksp,&(s->pc));
		PCFieldSplitSetIS(s->pc,"0",s->isg[0]);
		PCFieldSplitSetIS(s->pc,"1",s->isg[1]);

		// solve
		KSPSolve(s->ksp,s->RHS,s->x);
		VecNestGetSubVec(s->x,0,&s->tempvec2);
		VecCopy(s->tempvec2,s->inc_dxi);

		//update
		c = 1.0; c1  = -1.0;
		VecWAXPY(s->dxi,c,s->dxi_k,s->inc_dxi);
		VecWAXPY(s->tempvec,c1,s->dxi_old,s->dxi);

		c = (timem->delta - 1.0)*timem->dt;  c1 = 1/timem->delta/timem->dt;
		VecWAXPY(s->ddxi,c,s->ddxi_old,s->tempvec);
		VecScale(s->ddxi,c1);

		c = timem->dt;
		VecWAXPY(s->xi,c,s->dxi_old,s->xi_old);
		c = timem->theta*timem->dt*timem->dt;  c1 = (1-2*timem->theta)*timem->dt*timem->dt/2.0;
		VecAXPY(s->xi,c,s->ddxi);
		VecAXPY(s->xi,c1,s->ddxi_old);

		// residue
		c = -1.0;
		VecWAXPY(s->tempvec,c,s->xi_k,s->xi);
		VecNorm(s->tempvec,NORM_2,&residue1);
		//residue1 = 0.1*residue0;
		fieldCopy_current_to_k(s);

	}

	fieldCopy_current_to_old(s);
	printf("[Enter While loop] ended\n");

}

// output the field to day file
void SolidOutput(Field_S* s){

	VecView(s->xi,PETSC_VIEWER_STDOUT_WORLD);
}

void SolidFinalize(Field_S* s){

	printf("solid finalize");
}

// set the scattering from the global to local involved field
void setScattering(Field_S* s)
{

	Vec				local_x,local_y;
    IS 				is_xix,is_xiy;
    PetscInt		*l2g_x,*l2g_y;
    PetscInt		nx = s->ind.xix_N + s->ind.xix_ghoN, ny = s->ind.xiy_N + s->ind.xiy_ghoN;
    int				i;

    PetscMalloc1(nx,&l2g_x);
    PetscMalloc1(ny,&l2g_y);

    for(i = 0; i<nx; i++)
    	l2g_x[i] = s->ind.xix.l2g[i];

    for(i = 0; i<ny; i++)
        l2g_y[i] = s->ind.xiy.l2g[i];

	VecCreateSeq(PETSC_COMM_SELF,nx,&local_x);
    VecCreateSeq(PETSC_COMM_SELF,ny,&local_y);

    ISCreateGeneral(PETSC_COMM_SELF,nx,l2g_x,PETSC_OWN_POINTER,&is_xix);
    ISCreateGeneral(PETSC_COMM_SELF,ny,l2g_y,PETSC_OWN_POINTER,&is_xiy);

    VecScatterCreate(s->xi,is_xix,local_x,NULL,&s->scatter_x);
    VecScatterCreate(s->xi,is_xiy,local_y,NULL,&s->scatter_y);

    VecDestroy(&local_x); VecDestroy(&local_y);
    ISDestroy(&is_xix);  ISDestroy(&is_xiy);
}

void restart_file_output (Field_S* s)
{

	/*
	char  xi_file[] = "xi_data.dat";
	PetscViewer viewer;

    // write u
    PetscViewerCreate(PETSC_COMM_WORLD, &viewer);
    PetscViewerSetType(viewer,PETSCVIEWERASCII );
    PetscViewerFileSetMode(viewer, FILE_MODE_WRITE);
    PetscViewerFileSetName(viewer,xi_file);
    VecView(s->xi, viewer);
    PetscViewerDestroy(&viewer);
    */

}

void field_output (Field_S* s)
{
	char  filename_x[256] = "xi_data_";
	char  filename_y[256] = "eta_data_";
	char  ext[] = ".dat";
	char  rank_str[256];
	int		rank,i,ii,jj,cell;
	int 	nx = s->Nx;
	FILE    *fp_x, *fp_y;
	Vec		xix_local, xiy_local;
	double  *xix, *xiy;

	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	sprintf(rank_str,"%2d",rank);

	strcat(filename_x,rank_str);
	strcat(filename_x,ext);

	strcat(filename_y,rank_str);
	strcat(filename_y,ext);

	fp_x = fopen(filename_x, "w");
	fp_y = fopen(filename_y, "w");

	VecCreateSeq(PETSC_COMM_SELF,s->ind.xix_N +s->ind.xix_ghoN,&xix_local);
    VecCreateSeq(PETSC_COMM_SELF,s->ind.xiy_N +s->ind.xiy_ghoN,&xiy_local);

    VecScatterBegin	(s->scatter_x,s->xi,xix_local,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd	(s->scatter_x,s->xi,xix_local,INSERT_VALUES,SCATTER_FORWARD);

    VecScatterBegin	(s->scatter_y,s->xi,xiy_local,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd	(s->scatter_y,s->xi,xiy_local,INSERT_VALUES,SCATTER_FORWARD);

    VecGetArray(xix_local,&xix);
    VecGetArray(xiy_local,&xiy);

    fprintf(fp_x,"# %d\n",s->ind.xix_N);
    fprintf(fp_y,"# %d\n",s->ind.xiy_N);

    for(i=0;i<s->ind.xix_N ;i++){
    	cell = s->ind.xix.l2G[i];
    	ii = cell%nx;
    	jj = (int) floor((double)cell/nx);
    	fprintf(fp_x,"%d %d %.15e\n",ii,jj,xix[i]);
    }

    for(i=0;i<s->ind.xiy_N ;i++){
    	cell = s->ind.xiy.l2G[i];
    	ii = cell%nx;
    	jj = (int) floor((double)cell/nx);
    	fprintf(fp_y,"%d %d %.15e\n",ii,jj,xiy[i]);
    }

    VecRestoreArray(xix_local,&xix);
    VecRestoreArray(xiy_local,&xiy);

    VecDestroy(&xix_local);
    VecDestroy(&xiy_local);
	fclose(fp_x);
	fclose(fp_y);
}

// Copy the value in from *xi to *xi_old
void fieldCopy_current_to_old (Field_S* s)
{
	VecCopy(s->xi,s->xi_old);
	VecCopy(s->dxi,s->dxi_old);
	VecCopy(s->ddxi,s->ddxi_old);
}

// Copy the value in from *xi to *xi_k
void fieldCopy_current_to_k (Field_S* s)
{
	VecCopy(s->xi,s->xi_k);
	VecCopy(s->dxi,s->dxi_k);
	VecCopy(s->ddxi,s->ddxi_k);
}
