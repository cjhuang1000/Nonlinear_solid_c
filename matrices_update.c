#include "matrices_update.h"
#include "boun_func.h"
#include "nrutil.h"
#include "init_setting.h"

#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <petscmat.h>
#include <petscvec.h>
#include <petscis.h>

struct Grad_Dist
{
    double **center_xx, **center_yy;
    double **vertex_xy, **vertex_yx;
};

struct SSG_subgrid //stress, strain, gradient
{
    double xx[4], yy[4];
    double xy[4], yx[4];
};

struct Mcell0
{
    double KlS_xx[6][6];
    double KlS_xy[6][6];
    double KlS_yy[6][6];
    double KnlS_xx[6][6];
    double KnlS_yy[6][6];

    double FS_x[6];
    double FS_y[6];

};

/*Internal Subcell Elemental Matrices*/
boundfunc_cell  boundfunc_subcell0
    ={
        {{7,2,14,4,0,0},{2,1,4,2,0,0},{14,4,49,14,0,0},{4,2,14,7,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}},
        {{7,14,0,2,4,0},{14,49,0,4,14,0},{0,0,0,0,0,0},{2,4,0,1,2,0},{4,14,0,2,7,0},{0,0,0,0,0,0}},
        {{2,-2,4,-4,0,0},{-2,2,-4,4,0,0},{4,-4,14,-14,0,0},{-4,4,-14,14,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}},
        {{2,4,0,-2,-4,0},{4,14,0,-4,-14,0},{0,0,0,0,0,0},{-2,-4,0,2,4,0},{-4,-14,0,4,14,0},{0,0,0,0,0,0}},
        {{14,4,-14,-4,0,0},{4,2,-4,-2,0,0},{-14,-4,14,4,0,0},{-4,-2,4,2,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}},
        {{14,-14,0,4,-4,0},{-14,14,0,-4,4,0},{0,0,0,0,0,0},{4,-4,0,2,-2,0},{-4,4,0,-2,2,0},{0,0,0,0,0,0}},
        {{5,-5,0,1,-1,0},{-5,5,0,-1,1,0},{13,-13,0,5,-5,0},{-13,13,0,-5,5,0},{0,0,0,0,0,0},{0,0,0,0,0,0}},
        {{1,3,0,-1,-3,0},{-1,-3,0,1,3,0},{3,9,0,-3,-9,0},{-3,-9,0,3,9,0},{0,0,0,0,0,0},{0,0,0,0,0,0}},
        {{9,-9,0,3,-3,0},{3,-3,0,1,-1,0},{-9,9,0,-3,3,0},{-3,3,0,-1,1,0},{0,0,0,0,0,0},{0,0,0,0,0,0}},
        {{5,13,0,-5,-13,0},{1,5,0,-1,-5,0},{-5,-13,0,5,13,0},{-1,-5,0,1,5,0},{0,0,0,0,0,0},{0,0,0,0,0,0}},
        //{{5,-5,13,-13,0,0},{-5,5,-13,13,0,0},{0,0,0,0,0,0},{1,-1,5,-5,0,0},{-1,1,-5,5,0,0},{0,0,0,0,0,0}},
        //{{1,-1,3,-3,0,0},{3,-3,9,-9,0,0},{0,0,0,0,0,0},{-1,1,-3,3,0,0},{-3,3,-9,9,0,0},{0,0,0,0,0,0}},
        //{{9,3,-9,-3,0,0},{-9,-3,9,3,0,0},{0,0,0,0,0,0},{3,1,-3,-1,0,0},{-3,-1,3,1,0,0},{0,0,0,0,0,0}},
        //{{5,1,-5,-1,0,0},{13,5,-13,-5,0,0},{0,0,0,0,0,0},{-5,-1,5,1,0,0},{-13,-5,13,5,0,0},{0,0,0,0,0,0}},
        {{3,-3,9,-9,0,0},{1,-1,3,-3,0,0},{-3,3,-9,9,0,0},{-1,1,-3,3,0,0},{0,0,0,0,0,0},{0,0,0,0,0,0}},
        {{3,9,0,-3,-9,0},{-3,-9,0,3,9,0},{0,0,0,0,0,0},{1,3,0,-1,-3,0},{-1,-3,0,1,3,0},{0,0,0,0,0,0}},
        {-1, 1, -3, 3, 0, 0},
        {-3,-1,  3, 1, 0, 0},
        {-1,-3,  0, 1, 3, 0},
        {-3, 3,  0,-1, 1, 0},
    };

const boundary_function boundfunc
= {
    .coeff={ {  1, -1, -1,  1,  1, -1, -1,  1},
             { -1,  0, -1,  0,-.5, .5,-.5, .5},
             {-.5,-.5, .5, .5, -1, -1,  0,  0}},
    .sym[0]= {0,0,1,1,{0,1,2,3,4,5},{0,1,2,3,4,5},{0,1,2,3}},
    .sym[1]= {0,1,-1,1,{1,0,3,2,5,4},{2,1,0,5,4,3},{1,0,3,2}},
    .sym[2]= {1,0,1,-1,{4,5,2,3,0,1},{3,4,5,0,1,2},{2,3,0,1}},
    .sym[3]= {1,1,-1,-1,{5,4,3,2,1,0},{5,4,3,2,1,0},{3,2,1,0}}
    // sym has modified since c is starting from 0
    };

boundfunc_cutcell boudfunc_subcellint[4];
struct shapefunction shapefunc;
double m2pl,mpl,lambda,mu;



void compute_grad_update(struct Grad_Dist *ptr_dx, Field_S* ptr_f_old, Index_S* ptr_i, Grid_S* ptr_g);
struct Mcell0 compute_matricesCell_interior(struct SSG_subgrid stress, struct SSG_subgrid grad_xi);
struct Mcell0 compute_matricesCell_bound_update(int cell_i,struct SSG_subgrid stress, struct SSG_subgrid grad_xi);


void compute_matricesNonlinearStructure_update(Matrices_S* ptr_ms, Index_S* ptr_i, Grid_S* ptr_g, Solid* ptr_s, Field_S* ptr_f){

    int                 nx   = ptr_g->Nx,ny   = ptr_g->Ny;
    int                 i,j,k,cell_i;
    int                 index_xix[6], index_xiy[6];
    struct Grad_Dist    Dx;
    struct SSG_subgrid  stress;
    struct SSG_subgrid  strain;
    struct SSG_subgrid  grad_xi;
    struct Mcell0       mcell;
    struct invol        index_cell;

    Mat 				KLS_full, KNS_full; // before reduction
    Vec					FS_full, FS_temp;
    PetscScalar			temp1[36] = {0}, temp2[36] = {0}, temp3[36],temp4[36],temp5[36],temp6[36];
    PetscScalar			temp7[6], temp8[6];
    PetscErrorCode  ierr;

    m2pl = 2*ptr_s->mu + ptr_s->lambda;
    mpl = ptr_s->mu + ptr_s->lambda;
    lambda = ptr_s->lambda;
    mu = ptr_s->mu;


    /*initialize the temporary matrices @@ */

    //compute gradient

    Dx.center_xx = dmatrix(0,nx-1,0,ny-1);
    Dx.vertex_xy = dmatrix(0,nx-1,0,ny-1);
    Dx.vertex_yx = dmatrix(0,nx-1,0,ny-1);
    Dx.center_yy = dmatrix(0,nx-1,0,ny-1);

    compute_grad_update(&Dx, ptr_f, ptr_i, ptr_g);


    k = ptr_i->xix_N_before + ptr_i->xiy_N_before;

    MatCreate(PETSC_COMM_SELF,&(KLS_full));
    MatSetSizes(KLS_full,k,k,k,k);
    MatSetFromOptions(KLS_full);
    MatSetUp(KLS_full);
    MatZeroEntries(KLS_full);

    MatCreate(PETSC_COMM_SELF,&(KNS_full));
    MatSetSizes(KNS_full,k,k,k,k);
    MatSetFromOptions(KNS_full);
    MatSetUp(KNS_full);
    MatZeroEntries(KNS_full);

    VecCreate(PETSC_COMM_SELF,&FS_full);
    VecSetSizes(FS_full,k,PETSC_DECIDE);
    VecSetType(FS_full,"seq");
    VecZeroEntries(FS_full);

    // =========  for interior cells ===========
    for(i=0;i<ptr_i->cell_N_interior;i++)
    {
       cell_i = (int) ptr_i->cell_interior[i];
       index_cell = involvedIndices_grid(cell_i, nx);

       for(j=0;j<6;j++)
       {
            index_xix[j] = ptr_i->xix.G2g_before[index_cell.xix[j]];
            index_xiy[j] = ptr_i->xiy.G2g_before[index_cell.xiy[j]] + ptr_i->xix_N_before;
       }

        //calculate gradient in each subcell
       for(j=0;j<4;j++)
       {
           grad_xi.xx[j]=Dx.center_xx[index_cell.i][index_cell.j];
           grad_xi.yy[j]=Dx.center_yy[index_cell.i][index_cell.j];
       }
        grad_xi.xy[0] = Dx.vertex_xy[index_cell.i]  [index_cell.j];
        grad_xi.xy[1] = Dx.vertex_xy[index_cell.i+1][index_cell.j];
        grad_xi.xy[2] = Dx.vertex_xy[index_cell.i]  [index_cell.j+1];
        grad_xi.xy[3] = Dx.vertex_xy[index_cell.i+1][index_cell.j+1];

        grad_xi.yx[0] = Dx.vertex_yx[index_cell.i]  [index_cell.j];
        grad_xi.yx[1] = Dx.vertex_yx[index_cell.i+1][index_cell.j];
        grad_xi.yx[2] = Dx.vertex_yx[index_cell.i]  [index_cell.j+1];
        grad_xi.yx[3] = Dx.vertex_yx[index_cell.i+1][index_cell.j+1];

       //calculate the strain in each subcell
       for(j=0;j<4;j++)
       {
           strain.xx[j] = grad_xi.xx[j]+(grad_xi.xx[j]*grad_xi.xx[j] + grad_xi.yx[j]*grad_xi.yx[j])/2;
           strain.yy[j] = grad_xi.yy[j]+(grad_xi.xy[j]*grad_xi.xy[j] + grad_xi.yy[j]*grad_xi.yy[j])/2;
           strain.xy[j] = (grad_xi.xy[j] +grad_xi.yx[j] )/2+(grad_xi.xx[j]*grad_xi.xy[j] + grad_xi.yx[j]*grad_xi.yy[j])/2;
       }

       for(j=0;j<4;j++)
       {
           stress.xx[j]= m2pl*strain.xx[j]+lambda*strain.yy[j];
           stress.yy[j]= m2pl*strain.yy[j]+lambda*strain.xx[j];
           stress.xy[j]= 2*mu*strain.xy[j];
       }

       mcell = compute_matricesCell_interior(stress, grad_xi);

       for(j=0;j<6;j++)
        for(k=0;k<6;k++)
        {
        	temp1[j*6+k] = mcell.KnlS_xx[j][k];
        	temp2[j*6+k] = mcell.KnlS_yy[j][k];

        	temp3[j*6+k] = mcell.KlS_xx[j][k];
        	temp4[j*6+k] = mcell.KlS_xy[j][k];
        	temp5[j*6+k] = mcell.KlS_xy[k][j];
        	temp6[j*6+k] = mcell.KlS_yy[j][k];

        }

       MatSetValues(KNS_full,6,index_xix,6,index_xix,temp1,ADD_VALUES);
       MatSetValues(KNS_full,6,index_xiy,6,index_xiy,temp2,ADD_VALUES);
       MatSetValues(KLS_full,6,index_xix,6,index_xix,temp3,ADD_VALUES);
       MatSetValues(KLS_full,6,index_xix,6,index_xiy,temp4,ADD_VALUES);
       MatSetValues(KLS_full,6,index_xiy,6,index_xix,temp5,ADD_VALUES);
       MatSetValues(KLS_full,6,index_xiy,6,index_xiy,temp6,ADD_VALUES);

       for(j=0;j<6;j++)
       {
    	   temp7[j] = mcell.FS_x[j];
    	   temp8[j] = mcell.FS_y[j];
       }

       VecSetValues(FS_full,6,index_xix,temp7,ADD_VALUES);
       VecSetValues(FS_full,6,index_xiy,temp8,ADD_VALUES);

    }


    // =========  for boundary cells =========== @@

    for(i=0;i<ptr_i->cell_N_boundary; i++)
    //	for(i=0;i<1; i++)
    {
       cell_i = (int) ptr_i->cell_boundary[i];

       index_cell = involvedIndices_grid(cell_i, nx);

       for(j=0;j<6;j++)
       {
            index_xix[j] = ptr_i->xix.G2g_before[index_cell.xix[j]];
            index_xiy[j] = ptr_i->xiy.G2g_before[index_cell.xiy[j]]+ ptr_i->xix_N_before;
        }

       for(j=0;j<4;j++)
       {
            grad_xi.xx[j]=Dx.center_xx[index_cell.i][index_cell.j];
            grad_xi.yy[j]=Dx.center_yy[index_cell.i][index_cell.j];
       }
       grad_xi.xy[0] = Dx.vertex_xy[index_cell.i]  [index_cell.j];
       grad_xi.xy[1] = Dx.vertex_xy[index_cell.i+1][index_cell.j];
       grad_xi.xy[2] = Dx.vertex_xy[index_cell.i]  [index_cell.j+1];
       grad_xi.xy[3] = Dx.vertex_xy[index_cell.i+1][index_cell.j+1];

       grad_xi.yx[0] = Dx.vertex_yx[index_cell.i]  [index_cell.j];
       grad_xi.yx[1] = Dx.vertex_yx[index_cell.i+1][index_cell.j];
       grad_xi.yx[2] = Dx.vertex_yx[index_cell.i]  [index_cell.j+1];
       grad_xi.yx[3] = Dx.vertex_yx[index_cell.i+1][index_cell.j+1];

       //calculate the strain in each subcell
       for(j=0;j<4;j++)
       {
            strain.xx[j] = grad_xi.xx[j]+(grad_xi.xx[j]*grad_xi.xx[j] + grad_xi.yx[j]*grad_xi.yx[j])/2.0;
            strain.yy[j] = grad_xi.yy[j]+(grad_xi.xy[j]*grad_xi.xy[j] + grad_xi.yy[j]*grad_xi.yy[j])/2.0;
            strain.xy[j] = (grad_xi.xy[j] +grad_xi.yx[j] )/2+(grad_xi.xx[j]*grad_xi.xy[j] + grad_xi.yx[j]*grad_xi.yy[j])/2.0;
       }

       for(j=0;j<4;j++)
       {
            stress.xx[j]= m2pl*strain.xx[j]+lambda*strain.yy[j];
            stress.yy[j]= m2pl*strain.yy[j]+lambda*strain.xx[j];
            stress.xy[j]= 2.0*mu*strain.xy[j];
       }
       mcell = compute_matricesCell_bound_update(cell_i,stress,grad_xi);

       for(j=0;j<6;j++)
        for(k=0;k<6;k++)
        {
        	temp1[j*6+k] = mcell.KnlS_xx[j][k];
        	temp2[j*6+k] = mcell.KnlS_yy[j][k];

        	temp3[j*6+k] = mcell.KlS_xx[j][k];
        	temp4[j*6+k] = mcell.KlS_xy[j][k];
        	temp5[j*6+k] = mcell.KlS_xy[k][j];
        	temp6[j*6+k] = mcell.KlS_yy[j][k];
        }

       MatSetValues(KNS_full,6,index_xix,6,index_xix,temp1,ADD_VALUES);
       MatSetValues(KNS_full,6,index_xiy,6,index_xiy,temp2,ADD_VALUES);
       MatSetValues(KLS_full,6,index_xix,6,index_xix,temp3,ADD_VALUES);
       MatSetValues(KLS_full,6,index_xix,6,index_xiy,temp4,ADD_VALUES);
       MatSetValues(KLS_full,6,index_xiy,6,index_xix,temp5,ADD_VALUES);
       MatSetValues(KLS_full,6,index_xiy,6,index_xiy,temp6,ADD_VALUES);

       for(j=0;j<6;j++)
       {
    	   temp7[j] = mcell.FS_x[j];
    	   temp8[j] = mcell.FS_y[j];
       }

       VecSetValues(FS_full,6,index_xix,temp7,ADD_VALUES);
       VecSetValues(FS_full,6,index_xiy,temp8,ADD_VALUES);


    }

    MatAssemblyBegin(KNS_full,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(KNS_full,MAT_FINAL_ASSEMBLY);
    MatAssemblyBegin(KLS_full,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(KLS_full,MAT_FINAL_ASSEMBLY);
    VecAssemblyBegin(FS_full);
    VecAssemblyEnd(FS_full);

    //Set the reduced governing matrices
    //ISView(ptr_i->is_xi,PETSC_VIEWER_STDOUT_SELF);
    //printf("%i %i\n",ptr_i->xix_N,ptr_i->xiy_N);

    MatDestroy(&ptr_ms->KNS); MatDestroy(&ptr_ms->KLS);

    MatGetSubMatrix(KNS_full,ptr_i->is_xi,ptr_i->is_xi,MAT_INITIAL_MATRIX ,&ptr_ms->KNS);
    MatGetSubMatrix(KLS_full,ptr_i->is_xi,ptr_i->is_xi,MAT_INITIAL_MATRIX ,&ptr_ms->KLS);


    VecGetSubVector(FS_full,ptr_i->is_xi,&FS_temp);
    VecCopy(FS_temp,ptr_ms->FS);
    VecRestoreSubVector(FS_full,ptr_i->is_xi,&FS_temp);

    // free the memory
    MatDestroy(&KLS_full); MatDestroy(&KNS_full);
    VecDestroy(&FS_temp);  VecDestroy(&FS_full);

    free_dmatrix( Dx.center_xx, 0,nx-1,0,ny-1);
    free_dmatrix( Dx.vertex_xy, 0,nx-1,0,ny-1);
    free_dmatrix( Dx.vertex_yx, 0,nx-1,0,ny-1);
    free_dmatrix( Dx.center_yy, 0,nx-1,0,ny-1);

}

void compute_grad_update(struct Grad_Dist *ptr_dx, Field_S* ptr_f_old, Index_S* ptr_i, Grid_S* ptr_g)
{
    int     i,j;
    int     nx = ptr_g->Nx, ny = ptr_g->Ny,nn =ptr_g->N;
    double  *xigrid_x, *xigrid_y;
    double  dx = ptr_g->dx;
    PetscScalar *xi;


    xigrid_x = dvector(0,nn-1);
    xigrid_y = dvector(0,nn-1);
    VecGetArray(ptr_f_old->xi,&xi);

    for(i=0;i<nn;i++)
    {
        xigrid_x[i] =0.0;
        xigrid_y[i] =0.0;
    }

    for(i=0;i<ptr_i->xix_N;i++)
        xigrid_x[ptr_i->xix.l2G[i]]= xi[i];

    for(i=0;i<ptr_i->xiy_N;i++)
        xigrid_y[ptr_i->xiy.l2G[i]]= xi[i+ptr_i->xix_N];

    for(i=0;i<nx-1;i++)
        for(j=0;j<ny-1;j++)
        {
            ptr_dx->center_xx[i][j]=(xigrid_x[i+1+j*nx]-xigrid_x[i+j*nx])/dx;
            ptr_dx->center_yy[i][j]=(xigrid_y[i+(j+1)*nx]-xigrid_y[i+j*nx])/dx;

            if(i==0)  ptr_dx->vertex_yx[i][j] =(xigrid_y[1+j*nx]-xigrid_y[j*nx])/dx;
            else      ptr_dx->vertex_yx[i][j] =(xigrid_y[i+j*nx]-xigrid_y[i-1+j*nx])/dx;

            if(j==0)  ptr_dx->vertex_xy[i][j] =(xigrid_x[i+nx]-xigrid_x[i])/dx;
            else      ptr_dx->vertex_xy[i][j] =(xigrid_x[i+j*nx]-xigrid_x[i+(j-1)*nx])/dx;

        }
    free_dvector( xigrid_x, 0,nn-1);
    free_dvector( xigrid_y, 0,nn-1);
    VecRestoreArray(ptr_f_old->xi,&xi);

}

struct Mcell0 compute_matricesCell_interior(struct SSG_subgrid stress, struct SSG_subgrid grad_xi)
{
    struct Mcell0    mcell0={0};
    int              k,i,j;

    for(k=0; k<4; k++)
        for(i=0; i<6; i++)
        for(j=0; j<6; j++)
    {
        mcell0.KlS_xx[i][j] += (DSQR(1+grad_xi.xx[k])*m2pl  +DSQR(grad_xi.xy[k])*mu     )*boudfunc_subcellint[k].KS_hxxhxx[i][j]
                            +  (DSQR(1+grad_xi.xx[k])*mu    +DSQR(grad_xi.xy[k])*m2pl   )*boudfunc_subcellint[k].KS_hxyhxy[i][j]
                            +  (grad_xi.xy[k]*(1+grad_xi.xx[k])*mpl                    )*
                            (boudfunc_subcellint[k].KS_hxxhxy[i][j]+boudfunc_subcellint[k].KS_hxxhxy[j][i]) ;

        mcell0.KlS_yy[i][j] += (DSQR(1+grad_xi.yy[k])*m2pl  +DSQR(grad_xi.yx[k])*mu     )*boudfunc_subcellint[k].KS_hyyhyy[i][j]
                            +  (DSQR(1+grad_xi.yy[k])*mu    +DSQR(grad_xi.yx[k])*m2pl   )*boudfunc_subcellint[k].KS_hyxhyx[i][j]
                            +  (grad_xi.yx[k]*(1+grad_xi.yy[k])*mpl                    )*
                            (boudfunc_subcellint[k].KS_hyyhyx[i][j]+boudfunc_subcellint[k].KS_hyyhyx[j][i] ) ;

        mcell0.KlS_xy[i][j] += ( (1+grad_xi.xx[k])*grad_xi.yx[k]*m2pl       +(1+grad_xi.yy[k])*grad_xi.xy[k]*mu            )*boudfunc_subcellint[k].KS_hxxhyx[i][j]
                            +  ( (1+grad_xi.yy[k])*(1+grad_xi.xx[k])*lambda + grad_xi.xy[k]*grad_xi.yx[k]*mu               )*boudfunc_subcellint[k].KS_hxxhyy[i][j]
                            +  ( (1+grad_xi.yy[k])*(1+grad_xi.xx[k])*mu     + grad_xi.xy[k]*grad_xi.yx[k]*lambda           )*boudfunc_subcellint[k].KS_hxyhyx[i][j]
                            +  ( (1+grad_xi.yy[k])*grad_xi.xy[k]*m2pl       +(1+grad_xi.xx[k])*grad_xi.yx[k]*mu            )*boudfunc_subcellint[k].KS_hxyhyy[i][j]  ;

        mcell0.KnlS_xx[i][j]+= stress.xx[k]*boudfunc_subcellint[k].KS_hxxhxx[i][j]
                            +  stress.yy[k]*boudfunc_subcellint[k].KS_hxyhxy[i][j]
                            +  stress.xy[k]*(boudfunc_subcellint[k].KS_hxxhxy[i][j]+boudfunc_subcellint[k].KS_hxxhxy[j][i])  ;

        mcell0.KnlS_yy[i][j]+= stress.xx[k]*boudfunc_subcellint[k].KS_hyxhyx[i][j]
                            +  stress.yy[k]*boudfunc_subcellint[k].KS_hyyhyy[i][j]
                            +  stress.xy[k]*(boudfunc_subcellint[k].KS_hyyhyx[i][j]+boudfunc_subcellint[k].KS_hyyhyx[j][i])  ;
    }


    for(k=0; k<4; k++)
        for(i=0; i<6; i++)
    {
        mcell0.FS_x[i] += ((1+grad_xi.xx[k])*stress.xx[k] + stress.xy[k]*grad_xi.xy[k])*boudfunc_subcellint[k].FS_hxx[i]
                       +  ((1+grad_xi.xx[k])*stress.xy[k] + stress.yy[k]*grad_xi.xy[k])*boudfunc_subcellint[k].FS_hxy[i];

        mcell0.FS_y[i] += ((1+grad_xi.yy[k])*stress.xy[k] + stress.xx[k]*grad_xi.yx[k])*boudfunc_subcellint[k].FS_hyx[i]
                       +  ((1+grad_xi.yy[k])*stress.yy[k] + stress.xy[k]*grad_xi.yx[k])*boudfunc_subcellint[k].FS_hyy[i];
    }



    return mcell0;
}

struct Mcell0 compute_matricesCell_bound_update(int cell_i,struct SSG_subgrid stress, struct SSG_subgrid grad_xi)
{
    struct Mcell0    mcell0={0};
    int              k,i,j,n;

    for(k=0; k<4; k++)
    {

        if (shapefunc.info[cell_i][k] == -1) //inside
        {

            for(i=0; i<6; i++)
                for(j=0; j<6; j++)
            {
                mcell0.KlS_xx[i][j] += (DSQR(1+grad_xi.xx[k])*m2pl  +DSQR(grad_xi.xy[k])*mu     )*boudfunc_subcellint[k].KS_hxxhxx[i][j]
                                    +  (DSQR(1+grad_xi.xx[k])*mu    +DSQR(grad_xi.xy[k])*m2pl   )*boudfunc_subcellint[k].KS_hxyhxy[i][j]
                                    +  (grad_xi.xy[k]*(1+grad_xi.xx[k])*mpl                    )*
                                    (boudfunc_subcellint[k].KS_hxxhxy[i][j]+boudfunc_subcellint[k].KS_hxxhxy[j][i]) ;

                mcell0.KlS_yy[i][j] += (DSQR(1+grad_xi.yy[k])*m2pl  +DSQR(grad_xi.yx[k])*mu     )*boudfunc_subcellint[k].KS_hyyhyy[i][j]
                                    +  (DSQR(1+grad_xi.yy[k])*mu    +DSQR(grad_xi.yx[k])*m2pl   )*boudfunc_subcellint[k].KS_hyxhyx[i][j]
                                    +  (grad_xi.yx[k]*(1+grad_xi.yy[k])*mpl                    )*
                                    (boudfunc_subcellint[k].KS_hyyhyx[i][j]+boudfunc_subcellint[k].KS_hyyhyx[j][i]) ;

                mcell0.KlS_xy[i][j] += ( (1+grad_xi.xx[k])*grad_xi.yx[k]*m2pl       +(1+grad_xi.yy[k])*grad_xi.xy[k]*mu            )*boudfunc_subcellint[k].KS_hxxhyx[i][j]
                                    +  ( (1+grad_xi.yy[k])*(1+grad_xi.xx[k])*lambda + grad_xi.xy[k]*grad_xi.yx[k]*mu               )*boudfunc_subcellint[k].KS_hxxhyy[i][j]
                                    +  ( (1+grad_xi.yy[k])*(1+grad_xi.xx[k])*mu     + grad_xi.xy[k]*grad_xi.yx[k]*lambda           )*boudfunc_subcellint[k].KS_hxyhyx[i][j]
                                    +  ( (1+grad_xi.yy[k])*grad_xi.xy[k]*m2pl       +(1+grad_xi.xx[k])*grad_xi.yx[k]*mu            )*boudfunc_subcellint[k].KS_hxyhyy[i][j]  ;

                mcell0.KnlS_xx[i][j]+= stress.xx[k]*boudfunc_subcellint[k].KS_hxxhxx[i][j]
                                    +  stress.yy[k]*boudfunc_subcellint[k].KS_hxyhxy[i][j]
                                    +  stress.xy[k]*(boudfunc_subcellint[k].KS_hxxhxy[i][j]+boudfunc_subcellint[k].KS_hxxhxy[j][i])  ;

                mcell0.KnlS_yy[i][j]+= stress.xx[k]*boudfunc_subcellint[k].KS_hyxhyx[i][j]
                                    +  stress.yy[k]*boudfunc_subcellint[k].KS_hyyhyy[i][j]
                                    +  stress.xy[k]*(boudfunc_subcellint[k].KS_hyyhyx[i][j]+boudfunc_subcellint[k].KS_hyyhyx[j][i])  ;
            }


            for(i=0; i<6; i++)
            {
                mcell0.FS_x[i] += ((1+grad_xi.xx[k])*stress.xx[k] + stress.xy[k]*grad_xi.xy[k])*boudfunc_subcellint[k].FS_hxx[i]
                               +  ((1+grad_xi.xx[k])*stress.xy[k] + stress.yy[k]*grad_xi.xy[k])*boudfunc_subcellint[k].FS_hxy[i];

                mcell0.FS_y[i] += ((1+grad_xi.yy[k])*stress.xy[k] + stress.xx[k]*grad_xi.yx[k])*boudfunc_subcellint[k].FS_hyx[i]
                               +  ((1+grad_xi.yy[k])*stress.yy[k] + stress.xy[k]*grad_xi.yx[k])*boudfunc_subcellint[k].FS_hyy[i];
            }
        }
        else if(shapefunc.info[cell_i][k] > -1)
        {
            n = shapefunc.info[cell_i][k];
            for(i=0; i<6; i++)
                for(j=0; j<6; j++)
            {
                mcell0.KlS_xx[i][j] += (DSQR(1+grad_xi.xx[k])*m2pl  +DSQR(grad_xi.xy[k])*mu     )*shapefunc.cutcell[n].KS_hxxhxx[i][j]
                                    +  (DSQR(1+grad_xi.xx[k])*mu    +DSQR(grad_xi.xy[k])*m2pl   )*shapefunc.cutcell[n].KS_hxyhxy[i][j]
                                    +  (grad_xi.xy[k]*(1+grad_xi.xx[k])*mpl                    )*
                                    (shapefunc.cutcell[n].KS_hxxhxy[i][j]+shapefunc.cutcell[n].KS_hxxhxy[j][i]) ;

                mcell0.KlS_yy[i][j] += (DSQR(1+grad_xi.yy[k])*m2pl  +DSQR(grad_xi.yx[k])*mu     )*shapefunc.cutcell[n].KS_hyyhyy[i][j]
                                    +  (DSQR(1+grad_xi.yy[k])*mu    +DSQR(grad_xi.yx[k])*m2pl   )*shapefunc.cutcell[n].KS_hyxhyx[i][j]
                                    +  (grad_xi.yx[k]*(1+grad_xi.yy[k])*mpl                    )*
                                    (shapefunc.cutcell[n].KS_hyyhyx[i][j]+shapefunc.cutcell[n].KS_hyyhyx[j][i]) ;

                mcell0.KlS_xy[i][j] += ( (1+grad_xi.xx[k])*grad_xi.yx[k]*m2pl       +(1+grad_xi.yy[k])*grad_xi.xy[k]*mu            )*shapefunc.cutcell[n].KS_hxxhyx[i][j]
                                    +  ( (1+grad_xi.yy[k])*(1+grad_xi.xx[k])*lambda + grad_xi.xy[k]*grad_xi.yx[k]*mu               )*shapefunc.cutcell[n].KS_hxxhyy[i][j]
                                    +  ( (1+grad_xi.yy[k])*(1+grad_xi.xx[k])*mu     + grad_xi.xy[k]*grad_xi.yx[k]*lambda           )*shapefunc.cutcell[n].KS_hxyhyx[i][j]
                                    +  ( (1+grad_xi.yy[k])*grad_xi.xy[k]*m2pl       +(1+grad_xi.xx[k])*grad_xi.yx[k]*mu            )*shapefunc.cutcell[n].KS_hxyhyy[i][j]  ;

                mcell0.KnlS_xx[i][j]+= stress.xx[k]*shapefunc.cutcell[n].KS_hxxhxx[i][j]
                                    +  stress.yy[k]*shapefunc.cutcell[n].KS_hxyhxy[i][j]
                                    +  stress.xy[k]*
                                    (shapefunc.cutcell[n].KS_hxxhxy[i][j]+shapefunc.cutcell[n].KS_hxxhxy[j][i])  ;

                mcell0.KnlS_yy[i][j]+= stress.xx[k]*shapefunc.cutcell[n].KS_hyxhyx[i][j]
                                    +  stress.yy[k]*shapefunc.cutcell[n].KS_hyyhyy[i][j]
                                    +  stress.xy[k]*
                                    (shapefunc.cutcell[n].KS_hyyhyx[i][j]+shapefunc.cutcell[n].KS_hyyhyx[j][i])  ;
            }


            for(i=0; i<6; i++)
            {
                mcell0.FS_x[i] += ((1+grad_xi.xx[k])*stress.xx[k] + stress.xy[k]*grad_xi.xy[k])*shapefunc.cutcell[n].FS_hxx[i]
                               +  ((1+grad_xi.xx[k])*stress.xy[k] + stress.yy[k]*grad_xi.xy[k])*shapefunc.cutcell[n].FS_hxy[i];

                mcell0.FS_y[i] += ((1+grad_xi.yy[k])*stress.xy[k] + stress.xx[k]*grad_xi.yx[k])*shapefunc.cutcell[n].FS_hyx[i]
                               +  ((1+grad_xi.yy[k])*stress.yy[k] + stress.xy[k]*grad_xi.yx[k])*shapefunc.cutcell[n].FS_hyy[i];
            }
        }
        else
            continue;

    }

    return mcell0;
}
