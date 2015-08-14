#include "boun_func.h"

#include <petscvec.h>
#include <petscmat.h>
#include <petscis.h>

#include <math.h>
#include "nrutil.h"
#include <stdio.h>
#include <stdlib.h>

// variables
boundfunc_cell boundfunc_Mcell0
= {
    {{2,1,4,2,0,0},{1,2,2,4,0,0},{4,2,28,14,4,2},{2,4,14,28,2,4},{0,0,4,2,2,1},{0,0,2,4,1,2}},
    {{2,4,0,1,2,0},{4,28,4,2,14,2},{0,4,2,0,2,1},{1,2,0,2,4,0},{2,14,2,4,28,4},{0,2,1,0,4,2}},
    {{1,-1,2,-2,0,0},{-1,1,-2,2,0,0},{2,-2,14,-14,2,-2},{-2,2,-14,14,-2,2},{0,0,2,-2,1,-1},{0,0,-2,2,-1,1}},
    {{1,2,0,-1,-2,0},{2,14,2,-2,-14,-2},{0,2,1,0,-2,-1},{-1,-2,0,1,2,0},{-2,-14,-2,2,14,2},{0,-2,-1,0,2,1}},
    {{2,1,-2,-1,0,0},{1,2,-1,-2,0,0},{-2,-1,4,2,-2,-1},{-1,-2,2,4,-1,-2},{0,0,-2,-1,2,1},{0,0,-1,-2,1,2}},
    {{2,-2,0,1,-1,0},{-2,4,-2,-1,2,-1},{0,-2,2,0,-1,1},{1,-1,0,2,-2,0},{-1,2,-1,-2,4,-2},{0,-1,1,0,-2,2}},
    {{5,0,-5,1,0,-1},{-5,0,5,-1,0,1},{18,0,-18,18,0,-18},{-18,0,18,-18,0,18},{1,0,-1,5,0,-5},{-1,0,1,-5,0,5}},
    {{1,6,1,-1,-6,-1},{-1,-6,-1,1,6,1},{6,36,6,-6,-36,-6},{-6,-36,-6,6,36,6},{1,6,1,-1,-6,-1},{-1,-6,-1,1,6,1}},
    {{9,-6,-3,3,-2,-1},{3,6,-9,1,2,-3},{-6,4,2,6,-4,-2},{-2,-4,6,2,4,-6},{-3,2,1,-9,6,3},{-1,-2,3,-3,-6,9}},
    {{5,18,1,-5,-18,-1},{1,18,5,-1,-18,-5},{0,0,0,0,0,0},{0,0,0,0,0,0},{-5,-18,-1,5,18,1},{-1,-18,-5,1,18,5}},
    //{{5,-5,18,-18,1,-1},{0,0,0,0,0,0},{-5,5,-18,18,-1,1},{1,-1,18,-18,5,-5},{0,0,0,0,0,0},{-1,1,-18,18,-5,5}},
    //{{1,-1,6,-6,1,-1},{6,-6,36,-36,6,-6},{1,-1,6,-6,1,-1},{-1,1,-6,6,-1,1},{-6,6,-36,36,-6,6},{-1,1,-6,6,-1,1}},
    //{{9,3,-6,-2,-3,-1},{-6,6,4,-4,2,-2},{-3,-9,2,6,1,3},{3,1,6,2,-9,-3},{-2,2,-4,4,6,-6},{-1,-3,-2,-6,3,9}},
    //{{5,1,0,0,-5,-1},{18,18,0,0,-18,-18},{1,5,0,0,-1,-5},{-5,-1,0,0,5,1},{-18,-18,0,0,18,18},{-1,-5,0,0,1,5}},
    {{1,-1,3,-3,0,0},{1,-1,3,-3,0,0},{-1,1,0,0,1,-1},{-1,1,0,0,1,-1},{0,0,-3,3,-1,1},{0,0,-3,3,-1,1}},
    {{1,3,0,-1,-3,0},{-1,0,1,1,0,-1},{0,-3,-1,0,3,1},{1,3,0,-1,-3,0},{-1,0,1,1,0,-1},{0,-3,-1,0,3,1}},
    {-1, 1, -6, 6,-1, 1},
    {-1, -1, 0, 0, 1, 1},
    {-1, -6, -1, 1, 6, 1},
    {-1, 0,  1, -1, 0, 1}

};

/*Internal Subcell Elemental Matrices*/
extern boundfunc_cell  boundfunc_subcell0;
extern const boundary_function boundfunc;
extern boundfunc_cutcell boudfunc_subcellint[4];
extern struct shapefunction shapefunc;

static double dx;
static double dx2;

struct boundfunc_edge{

	double MS_xx[4][4];
	double MS_yy[4][4];
	double KS_hxxhxx[4][4];
	double KS_hyyhyy[4][4];
	double KS_hxyhxy[4][4];
	double KS_hyxhyx[4][4];
	double KS_hxxhyx[4][4];
	double KS_hxxhyy[4][4];
	double KS_hxyhyx[4][4];
	double KS_hxyhyy[4][4];
	double KS_hxxhxy[4][4];
	double KS_hyyhyx[4][4];
	double FS_hxx[4];
	double FS_hxy[4];
	double FS_hyy[4];
	double FS_hyx[4];
};

struct boundfunc_edgeNS{

	double NS_xc[4];
	double NS_yc[4];
};

/// functions
double boundFunct_MS_alpha(int k, int l,double x,double y,double alpha);
double boundFunct_MS_gamma(int k, int l,double x,double y,double gamma);

double boundFunct_KS_xx_alpha(int k, int l,double y1,double y2,double alpha, double beta);
double boundFunct_KS_xx_gamma(int k, int l,double x1,double x2,double gamma, double delta);
double boundFunct_KS_yy_alpha(int k, int l,double y1,double y2,double alpha, double beta);
double boundFunct_KS_yy_gamma(int k, int l,double x1,double x2,double gamma, double delta);
double boundFunct_KS_xy_alpha(int k, int l,double y1,double y2,double alpha, double beta);
double boundFunct_KS_xy_gamma(int k, int l,double x1,double x2,double gamma, double delta);

double boundFunct_FS_x_alpha(int k,double y1,double y2,double alpha, double beta);
double boundFunct_FS_x_gamma(int k,double x1,double x2,double gamma, double delta);
double boundFunct_FS_y_alpha(int k,double y1,double y2,double alpha, double beta);
double boundFunct_FS_y_gamma(int k,double x1,double x2,double gamma, double delta);

double boundFunct_NS_alpha(int l, double x, double y, double alpha);
double boundFunct_NS_gamma(int l, double x, double y, double gamma);

struct boundfunc_edge   compute_matricesEdge   (double start_coord[2], double end_coord[2]);
struct boundfunc_edgeNS compute_matricesEdgeNS (double start_coord[2], double end_coord[2]);
int              type_determination     (double x[4]);
boundfunc_cell   compute_matricesCell_bound(int cell_i, double sdf_value[4], char bd_type, int* ptr_bd);
boundfunc_cell   compute_matricesSubcell(double bound_value[4], char bd_type, int* ptr_type);

void reduce_vector(short int** a, int n_b, int n_a);
void reduce_system(Grid_S* ptr_g, Index_S* ptr_i, Mat* ptr_MS, Mat* ptr_KLS, Mat* ptr_NS, Mat* ptr_DS);

/// Setting up the
void set_boundfunc(double ddx){

    int     i,j,subcell_k;
    double temp1[6][6],temp2[6][6];

    dx = ddx;
    dx2 = ddx*ddx;

    /* Set the Cell Elemental Matrices with. Adding dx^2 and denominators  */
    for(i=0;i<6;i++){
        for(j=0;j<6;j++)
        {
            boundfunc_Mcell0.MS_xx[i][j]     *= (dx2/144);
            boundfunc_Mcell0.MS_yy[i][j]     *= (dx2/144);
            boundfunc_Mcell0.KS_hxxhxx[i][j] /= 24;
            boundfunc_Mcell0.KS_hyyhyy[i][j] /= 24;
            boundfunc_Mcell0.KS_hxyhxy[i][j] /= 12;
            boundfunc_Mcell0.KS_hyxhyx[i][j] /= 12;
            boundfunc_Mcell0.KS_hxxhyx[i][j] /= 96;
            boundfunc_Mcell0.KS_hxxhyy[i][j] /= 64;
            boundfunc_Mcell0.KS_hxyhyx[i][j] /= 64;
            boundfunc_Mcell0.KS_hxyhyy[i][j] /= 96;
            boundfunc_Mcell0.KS_hxxhxy[i][j] /= 16;
            boundfunc_Mcell0.KS_hyyhyx[i][j] /= 16;

            boundfunc_subcell0.MS_xx[i][j]     *= (dx2/576);
            boundfunc_subcell0.MS_yy[i][j]     *= (dx2/576);
            boundfunc_subcell0.KS_hxxhxx[i][j] /= 96;
            boundfunc_subcell0.KS_hyyhyy[i][j] /= 96;
            boundfunc_subcell0.KS_hxyhxy[i][j] /= 96;
            boundfunc_subcell0.KS_hyxhyx[i][j] /= 96;
            boundfunc_subcell0.KS_hxxhyx[i][j] /= 96;
            boundfunc_subcell0.KS_hxxhyy[i][j] /= 64;
            boundfunc_subcell0.KS_hxyhyx[i][j] /= 64;
            boundfunc_subcell0.KS_hxyhyy[i][j] /= 96;
            boundfunc_subcell0.KS_hxxhxy[i][j] /= 64;
            boundfunc_subcell0.KS_hyyhyx[i][j] /= 64;

        }

        boundfunc_Mcell0.FS_hxx[i] *= dx/8;
        boundfunc_Mcell0.FS_hxy[i] *= dx/4;
        boundfunc_Mcell0.FS_hyy[i] *= dx/8;
        boundfunc_Mcell0.FS_hyx[i] *= dx/4;

        boundfunc_subcell0.FS_hxx[i] *= dx/16;
        boundfunc_subcell0.FS_hxy[i] *= dx/16;
        boundfunc_subcell0.FS_hyy[i] *= dx/16;
        boundfunc_subcell0.FS_hyx[i] *= dx/16;
    }


    for(subcell_k=0; subcell_k<4;subcell_k++)
    {
        for(i=0;i<6;i++){
            for(j=0;j<6;j++)
            {

                boudfunc_subcellint[subcell_k].KS_hxxhxx[i][j]
                = boundfunc_subcell0.KS_hxxhxx[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].x[j]];
                boudfunc_subcellint[subcell_k].KS_hyyhyy[i][j]
                = boundfunc_subcell0.KS_hyyhyy[boundfunc.sym[subcell_k].y[i]][boundfunc.sym[subcell_k].y[j]];
                boudfunc_subcellint[subcell_k].KS_hxyhxy[i][j]
                = boundfunc_subcell0.KS_hxyhxy[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].x[j]];
                boudfunc_subcellint[subcell_k].KS_hyxhyx[i][j]
                = boundfunc_subcell0.KS_hyxhyx[boundfunc.sym[subcell_k].y[i]][boundfunc.sym[subcell_k].y[j]];

                boudfunc_subcellint[subcell_k].KS_hxxhyx[i][j]
                = boundfunc_subcell0.KS_hxxhyx[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].y[j]];
                boudfunc_subcellint[subcell_k].KS_hxxhyy[i][j]
                = boundfunc_subcell0.KS_hxxhyy[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].y[j]]
                *boundfunc.sym[subcell_k].signx*boundfunc.sym[subcell_k].signy;
                boudfunc_subcellint[subcell_k].KS_hxyhyx[i][j]
                = boundfunc_subcell0.KS_hxyhyx[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].y[j]]
                *boundfunc.sym[subcell_k].signx*boundfunc.sym[subcell_k].signy;;
                boudfunc_subcellint[subcell_k].KS_hxyhyy[i][j]
                = boundfunc_subcell0.KS_hxyhyy[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].y[j]];

                boudfunc_subcellint[subcell_k].KS_hxxhxy[i][j]
                = boundfunc_subcell0.KS_hxxhxy[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].x[j]]
                *boundfunc.sym[subcell_k].signx*boundfunc.sym[subcell_k].signy;
                boudfunc_subcellint[subcell_k].KS_hyyhyx[i][j]
                = boundfunc_subcell0.KS_hyyhyx[boundfunc.sym[subcell_k].y[i]][boundfunc.sym[subcell_k].y[j]]
                *boundfunc.sym[subcell_k].signx*boundfunc.sym[subcell_k].signy;
            }
            boudfunc_subcellint[subcell_k].FS_hxx[i] = boundfunc_subcell0.FS_hxx[boundfunc.sym[subcell_k].x[i]]* boundfunc.sym[subcell_k].signx ;
            boudfunc_subcellint[subcell_k].FS_hxy[i] = boundfunc_subcell0.FS_hxy[boundfunc.sym[subcell_k].x[i]]* boundfunc.sym[subcell_k].signy ;
            boudfunc_subcellint[subcell_k].FS_hyy[i] = boundfunc_subcell0.FS_hyy[boundfunc.sym[subcell_k].y[i]]* boundfunc.sym[subcell_k].signy ;
            boudfunc_subcellint[subcell_k].FS_hyx[i] = boundfunc_subcell0.FS_hyx[boundfunc.sym[subcell_k].y[i]]* boundfunc.sym[subcell_k].signx ;
        }
        /*
        for(i=0;i<6;i++)
            for(j=0;j<6;j++){
                temp1[i][j] = boudfunc_subcellint[subcell_k].KS_hxxhxy[j][i];
                temp2[i][j] = boudfunc_subcellint[subcell_k].KS_hyyhyx[j][i];
            }

         for(i=0;i<6;i++)
            for(j=0;j<6;j++){
                boudfunc_subcellint[subcell_k].KS_hxxhxy[i][j] += temp1[i][j];
                boudfunc_subcellint[subcell_k].KS_hyyhyx[i][j] += temp2[i][j];
            }*/
        }

}



/// first time to calculate governing matrices for reducing system: assuming undeformed state
/// Kl, M, DS, NS are constructed
void compute_matricesNonlinearStructure(Matrices_S* ptr_ms, Index_S* ptr_i, Grid_S* ptr_g, Solid* ptr_s, char* fnd){

    double          m2pl = 2*ptr_s->mu + ptr_s->lambda; // mu*2 plus lambda
    double          sdf[4]; // signed distance function
    
    int             nx   = ptr_g->Nx;
    int             i,j,k, ii, jj;
    int             cell_i,nsubcell,bdcell = 0;
    int             index_xix[6], index_xiy[6], index_stagx[2],index_stagy[2];
    
    struct invol    index_cell;
    boundfunc_cell  Mcell;

    Mat             MS_full, KLS_full, NS_full, DS_full;
    PetscScalar     temp1[36],temp2[36],temp3[36],temp4[36],temp5[36],temp6[36]; // temporary arrays for MAT construction 
    PetscScalar     temp7[12],temp8[12];

    /*initialize the matrices MS, KLS, NS, DS */

    k = ptr_i->xi_gloN;

    MatCreate(MPI_COMM_WORLD,&(MS_full));
    MatSetSizes(MS_full,PETSC_DECIDE,PETSC_DECIDE,k,k);
    MatSetFromOptions(MS_full);
    MatSetUp(MS_full);
    MatZeroEntries(MS_full);

    MatCreate(MPI_COMM_WORLD,&(KLS_full));
    MatSetSizes(KLS_full,PETSC_DECIDE,PETSC_DECIDE,k,k);
    MatSetFromOptions(KLS_full);
    MatSetUp(KLS_full);
    MatZeroEntries(KLS_full);

    MatCreate(MPI_COMM_WORLD,&(NS_full));
    MatSetSizes(NS_full,PETSC_DECIDE,PETSC_DECIDE,k,ptr_i->xi_gloNcell_Neumann);
    MatSetFromOptions(NS_full);
    MatSetUp(NS_full);
    MatZeroEntries(NS_full);

    MatCreate(MPI_COMM_WORLD,&(DS_full));
    MatSetSizes(DS_full,PETSC_DECIDE,PETSC_DECIDE,k,ptr_i->xi_gloNcell_Dirichlet);
    MatSetFromOptions(DS_full);
    MatSetUp(DS_full);
    MatZeroEntries(DS_full);

    // =========  for interior cells ===========
    for(i=0;i<ptr_i->cell_N_interior;i++)
    {
       cell_i = (int) ptr_i->cell_interior[i];
       index_cell = involvedIndices_grid(cell_i, nx);

       for(j=0;j<6;j++)
       {
            index_xix[j] = ptr_i->xix.G2g[index_cell.xix[j]];
            index_xiy[j] = ptr_i->xiy.G2g[index_cell.xiy[j]];  // y component if after x component
       }

       for(j=0;j<6;j++)
        for(k=0;k<6;k++)
        {
            temp1[j*6+k] = boundfunc_Mcell0.MS_xx[j][k];
            temp2[j*6+k] = boundfunc_Mcell0.MS_yy[j][k];

            temp3[j*6+k] = m2pl*boundfunc_Mcell0.KS_hxxhxx[j][k]
                           + ptr_s->mu*boundfunc_Mcell0.KS_hxyhxy[j][k];
            temp4[j*6+k] = ptr_s->lambda*boundfunc_Mcell0.KS_hxxhyy[j][k]
                           + ptr_s->mu*boundfunc_Mcell0.KS_hxyhyx[j][k];
            temp5[j*6+k] = ptr_s->lambda*boundfunc_Mcell0.KS_hxxhyy[k][j]
                           + ptr_s->mu*boundfunc_Mcell0.KS_hxyhyx[k][j];
            temp6[j*6+k] = m2pl*boundfunc_Mcell0.KS_hyyhyy[j][k]
                           + ptr_s->mu*boundfunc_Mcell0.KS_hyxhyx[j][k];
        }

    MatSetValues(MS_full ,6,index_xix,6,index_xix,temp1,ADD_VALUES);
    MatSetValues(MS_full ,6,index_xiy,6,index_xiy,temp2,ADD_VALUES);
    MatSetValues(KLS_full,6,index_xix,6,index_xix,temp3,ADD_VALUES);
    MatSetValues(KLS_full,6,index_xix,6,index_xiy,temp4,ADD_VALUES);
    MatSetValues(KLS_full,6,index_xiy,6,index_xix,temp5,ADD_VALUES);
    MatSetValues(KLS_full,6,index_xiy,6,index_xiy,temp6,ADD_VALUES);

    }

    //Determine how many subcells are intersecting with the boundaries
    nsubcell = 0;
    shapefunc.info = imatrix(0,ptr_g->N-1,0,3);

    for(i=0;i<ptr_i->cell_N_boundary ;i++)
    {
        cell_i = (int) ptr_i->cell_boundary[i];
        ii = cell_i%nx;
        jj = (int) floor((double)cell_i/nx);

        sdf[0] = ptr_s->boundary_value[ii][jj];
        sdf[1] = ptr_s->boundary_value[ii+1][jj];
        sdf[2] = ptr_s->boundary_value[ii][jj+1];
        sdf[3] = ptr_s->boundary_value[ii+1][jj+1];

        nsubcell+= type_determination(sdf);
    }

    shapefunc.cutcell=(boundfunc_cutcell*) malloc(nsubcell*sizeof(boundfunc_cutcell));

    for(i=0;i<ptr_g->N;i++)
        for(j=0;j<4;j++)
            shapefunc.info[i][j]=-3; // initial value: -3
                                     //   not a single corner in :-2
                                     //   the whole cell is in: -1
                                     //   cut-subcell >=0

    // =========  for boundary cells ===========
    for(i=0;i<ptr_i->cell_N_boundary; i++)
    	//for(i=0;i<1; i++)
    {
       cell_i = (int) ptr_i->cell_boundary[i];

       index_cell = involvedIndices_grid(cell_i, nx);
       ii = index_cell.i;
       jj = index_cell.j;

       for(j=0;j<6;j++)
       {
            index_xix[j] = ptr_i->xix.G2g[index_cell.xix[j]];
            index_xiy[j] = ptr_i->xiy.G2g[index_cell.xiy[j]]; // y component if after x component
       }

       sdf[0] = ptr_s->boundary_value[ii][jj];
       sdf[1] = ptr_s->boundary_value[ii+1][jj];
       sdf[2] = ptr_s->boundary_value[ii][jj+1];
       sdf[3] = ptr_s->boundary_value[ii+1][jj+1];


       Mcell = compute_matricesCell_bound(cell_i, sdf , fnd[cell_i] ,&bdcell);

       for(j=0;j<6;j++)
        for(k=0;k<6;k++)
        {
            temp1[j*6+k] = Mcell.MS_xx[j][k];
            temp2[j*6+k] = Mcell.MS_yy[j][k];

            temp3[j*6+k] = m2pl*Mcell.KS_hxxhxx[j][k]
                           + ptr_s->mu*Mcell.KS_hxyhxy[j][k];
            temp4[j*6+k] = ptr_s->lambda*Mcell.KS_hxxhyy[j][k]
                           + ptr_s->mu*Mcell.KS_hxyhyx[j][k];
            temp5[j*6+k] = ptr_s->lambda*Mcell.KS_hxxhyy[k][j]
                           + ptr_s->mu*Mcell.KS_hxyhyx[k][j];
            temp6[j*6+k] = m2pl*Mcell.KS_hyyhyy[j][k]
                           + ptr_s->mu*Mcell.KS_hyxhyx[j][k];

        }

        MatSetValues(MS_full ,6,index_xix,6,index_xix,temp1,ADD_VALUES);
        MatSetValues(MS_full ,6,index_xiy,6,index_xiy,temp2,ADD_VALUES);
        MatSetValues(KLS_full,6,index_xix,6,index_xix,temp3,ADD_VALUES);
        MatSetValues(KLS_full,6,index_xix,6,index_xiy,temp4,ADD_VALUES);
        MatSetValues(KLS_full,6,index_xiy,6,index_xix,temp5,ADD_VALUES);
        MatSetValues(KLS_full,6,index_xiy,6,index_xiy,temp6,ADD_VALUES);

       if (fnd[cell_i] == 'D')
       {
           for(j=0;j<2;j++)
            {
                index_stagx[j] = ptr_i->xix.C2c_dirichlet[index_cell.stagx_cell[j]];
                index_stagy[j] = ptr_i->xiy.C2c_dirichlet[index_cell.stagy_cell[j]];
            }
            for(j=0;j<6;j++)
                for(k=0;k<2;k++)
                {
                    temp7[j*2+k] = Mcell.DS_xc[j][k];
                    temp8[j*2+k] = Mcell.DS_yc[j][k];
                }

            MatSetValues(DS_full,6,index_xix,2,index_stagx,temp7,ADD_VALUES);
            MatSetValues(DS_full,6,index_xiy,2,index_stagy,temp8,ADD_VALUES);

       }
       else if (fnd[cell_i] == 'N')
       {
            for(j=0;j<2;j++)
            {
                index_stagx[j] = ptr_i->xix.C2c_neumann[index_cell.stagx_cell[j]];
                index_stagy[j] = ptr_i->xiy.C2c_neumann[index_cell.stagy_cell[j]];

            }
            for(j=0;j<6;j++)
                for(k=0;k<2;k++)
                {
                    temp7[j*2+k] =  Mcell.NS_xc[j][k];
                    temp8[j*2+k] =  Mcell.NS_yc[j][k];
                }

            MatSetValues(NS_full,6,index_xix,2,index_stagx,temp7,ADD_VALUES);
            MatSetValues(NS_full,6,index_xiy,2,index_stagy,temp8,ADD_VALUES);
        }
    }

    //Matrices assembly
    MatAssemblyBegin(MS_full,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(MS_full,MAT_FINAL_ASSEMBLY);

    MatAssemblyBegin(KLS_full,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(KLS_full,MAT_FINAL_ASSEMBLY);

    MatAssemblyBegin(NS_full,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(NS_full,MAT_FINAL_ASSEMBLY);

    MatAssemblyBegin(DS_full,MAT_FINAL_ASSEMBLY);
    MatAssemblyEnd(DS_full,MAT_FINAL_ASSEMBLY);

    //Reduce system

    reduce_system(ptr_g,ptr_i,&MS_full,&KLS_full,&NS_full,&DS_full);
    /*
    //Set the reduced governing matrices
    MatGetSubMatrix(MS_full,ptr_i->is_xi,ptr_i->is_xi,MAT_INITIAL_MATRIX,&ptr_ms->MS);
    MatGetSubMatrix(KLS_full,ptr_i->is_xi,ptr_i->is_xi,MAT_INITIAL_MATRIX,&ptr_ms->KLS);
    MatGetSubMatrix(KLS_full,ptr_i->is_xi,ptr_i->is_xi,MAT_INITIAL_MATRIX,&ptr_ms->KNS); // just for initialization
    MatGetSubMatrix(NS_full,ptr_i->is_xi,ptr_i->is_neu,MAT_INITIAL_MATRIX,&ptr_ms->NS);
    MatGetSubMatrix(DS_full,ptr_i->is_xi,ptr_i->is_dir,MAT_INITIAL_MATRIX,&ptr_ms->DS);
     */
    // free the memory
    MatDestroy(&MS_full); MatDestroy(&KLS_full); MatDestroy(&NS_full); MatDestroy(&DS_full);
}

// input the sdf of four corners and return how many subcells is intersecting with boundaries
int type_determination(double x[4])
{
    int         number=0,i;
    double      sdf[9] = {x[0],0,x[1],0,0,0,x[2],0,x[3]};
    int         y[4][4] = {{0,1,3,4},{1,2,4,5},{3,4,6,7},{4,5,7,8}},ssdf[9];
    _Bool       isext,isint;

    sdf[1] = (sdf[0]+sdf[2])/2;
    sdf[3] = (sdf[0]+sdf[6])/2;
    sdf[5] = (sdf[2]+sdf[8])/2;
    sdf[7] = (sdf[6]+sdf[8])/2;
    sdf[4] = (sdf[0]+sdf[2]+sdf[6]+sdf[8])/4;

    for(i=0;i<9;i++)
        ssdf[i] = (int) SIGN(1,sdf[i]);

    for(i=0;i<4;i++)
    {
        isext = ((ssdf[y[i][0]]!= -1) && (ssdf[y[i][1]]!= -1) && (ssdf[y[i][2]]!= -1) && (ssdf[y[i][3]]!= -1));
        isint = ((ssdf[y[i][0]]+ssdf[y[i][1]]+ssdf[y[i][2]]+ssdf[y[i][3]])<=-3);
        if (!isext && !isint)
            number +=1;
    }
    return number;
}

boundfunc_cell compute_matricesCell_bound(int cell_i, double sdf_value[4], char bd_type, int* ptr_bd)
{
    // ptr_bd: the number of cut-subcell
    boundfunc_cell  Mcell={0};
    int             subcell_k,i,j;
    int             sym[4][4] = {{0,1,4,3},{2,1,4,5},{6,7,4,3},{8,7,4,5}};
    double          sdf[9] = {sdf_value[0],0,sdf_value[1],0,0,0,sdf_value[2],0,sdf_value[3]};
    double          bound_value[4];
    boundfunc_cell  subcell;
    int             type,n;

    sdf[1] = (sdf[0]+sdf[2])/2;
    sdf[3] = (sdf[0]+sdf[6])/2;
    sdf[5] = (sdf[2]+sdf[8])/2;
    sdf[7] = (sdf[6]+sdf[8])/2;
    sdf[4] = (sdf[0]+sdf[2]+sdf[6]+sdf[8])/4;

    for(subcell_k=0;subcell_k<4; subcell_k++ )
    {
        for(i=0; i<4; i++)
            bound_value[i]=sdf[sym[subcell_k][i]];

        subcell=compute_matricesSubcell(bound_value,bd_type,&type);

        for(i=0; i<6; i++)
            for(j=0; j<6; j++)
            {
                Mcell.MS_xx[i][j]+=subcell.MS_xx[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].x[j]];
                Mcell.MS_yy[i][j]+=subcell.MS_yy[boundfunc.sym[subcell_k].y[i]][boundfunc.sym[subcell_k].y[j]];

                Mcell.KS_hxxhxx[i][j] +=subcell.KS_hxxhxx[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].x[j]];
                Mcell.KS_hxyhxy[i][j] +=subcell.KS_hxyhxy[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].x[j]];

                Mcell.KS_hxxhyy[i][j] +=subcell.KS_hxxhyy[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].y[j]]
                                        *boundfunc.sym[subcell_k].signx*boundfunc.sym[subcell_k].signy;
                Mcell.KS_hxyhyx[i][j] +=subcell.KS_hxyhyx[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].y[j]]
                                        *boundfunc.sym[subcell_k].signx*boundfunc.sym[subcell_k].signy;

                Mcell.KS_hyyhyy[i][j] +=subcell.KS_hyyhyy[boundfunc.sym[subcell_k].y[i]][boundfunc.sym[subcell_k].y[j]];
                Mcell.KS_hyxhyx[i][j] +=subcell.KS_hyxhyx[boundfunc.sym[subcell_k].y[i]][boundfunc.sym[subcell_k].y[j]];

                Mcell.KS_hxxhxy[i][j] +=subcell.KS_hxxhxy[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].x[j]]
                                        *boundfunc.sym[subcell_k].signx*boundfunc.sym[subcell_k].signy;
                Mcell.KS_hyyhyx[i][j] +=subcell.KS_hyyhyx[boundfunc.sym[subcell_k].y[i]][boundfunc.sym[subcell_k].y[j]]
                                        *boundfunc.sym[subcell_k].signx*boundfunc.sym[subcell_k].signy;
                Mcell.KS_hxxhyx[i][j] +=subcell.KS_hxxhyx[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].y[j]];
                Mcell.KS_hxyhyy[i][j] +=subcell.KS_hxyhyy[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].y[j]];
            }

        if(bd_type == 'N')
        {
            for(i=0; i<6; i++)
            {
                Mcell.NS_xc[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].cellx] += subcell.NS_xc[i][0];
                Mcell.NS_yc[boundfunc.sym[subcell_k].y[i]][boundfunc.sym[subcell_k].celly] += subcell.NS_yc[i][0];
            }

        }
        else if (bd_type == 'D')
        {
             for(i=0; i<6; i++)
            {
                Mcell.DS_xc[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].cellx] += subcell.DS_xc[i][0];
                Mcell.DS_yc[boundfunc.sym[subcell_k].y[i]][boundfunc.sym[subcell_k].celly] += subcell.DS_yc[i][0];
            }
        }

        if (type == 1)
        {
            n = *ptr_bd;
            for(i=0; i<6; i++)
                for(j=0; j<6; j++)
                {
                    shapefunc.cutcell[n].KS_hxxhxx[i][j]=subcell.KS_hxxhxx[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].x[j]];
                    shapefunc.cutcell[n].KS_hxyhxy[i][j]=subcell.KS_hxyhxy[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].x[j]];

                    shapefunc.cutcell[n].KS_hxxhyy[i][j]=subcell.KS_hxxhyy[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].y[j]]
                                                        *boundfunc.sym[subcell_k].signx*boundfunc.sym[subcell_k].signy;
                    shapefunc.cutcell[n].KS_hxyhyx[i][j]=subcell.KS_hxyhyx[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].y[j]]
                                                        *boundfunc.sym[subcell_k].signx*boundfunc.sym[subcell_k].signy;
                    shapefunc.cutcell[n].KS_hyyhyy[i][j]=subcell.KS_hyyhyy[boundfunc.sym[subcell_k].y[i]][boundfunc.sym[subcell_k].y[j]];
                    shapefunc.cutcell[n].KS_hyxhyx[i][j]=subcell.KS_hyxhyx[boundfunc.sym[subcell_k].y[i]][boundfunc.sym[subcell_k].y[j]];

                    shapefunc.cutcell[n].KS_hxxhxy[i][j]=subcell.KS_hxxhxy[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].x[j]]
                                                        //+ subcell.KS_hxxhxy[boundfunc.sym[subcell_k].x[j]][boundfunc.sym[subcell_k].x[i]])
                                                        *boundfunc.sym[subcell_k].signx*boundfunc.sym[subcell_k].signy;
                    shapefunc.cutcell[n].KS_hyyhyx[i][j]=subcell.KS_hyyhyx[boundfunc.sym[subcell_k].y[i]][boundfunc.sym[subcell_k].y[j]]
                                                        //+ subcell.KS_hyyhyx[boundfunc.sym[subcell_k].y[j]][boundfunc.sym[subcell_k].y[i]] )
                                                        *boundfunc.sym[subcell_k].signx*boundfunc.sym[subcell_k].signy;

                    shapefunc.cutcell[n].KS_hxxhyx[i][j]=subcell.KS_hxxhyx[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].y[j]];
                    shapefunc.cutcell[n].KS_hxyhyy[i][j]=subcell.KS_hxyhyy[boundfunc.sym[subcell_k].x[i]][boundfunc.sym[subcell_k].y[j]];

                }

            for(i=0; i<6; i++)
            {
                shapefunc.cutcell[n].FS_hxx[i] = subcell.FS_hxx[boundfunc.sym[subcell_k].x[i]]*boundfunc.sym[subcell_k].signx;
                shapefunc.cutcell[n].FS_hxy[i] = subcell.FS_hxy[boundfunc.sym[subcell_k].x[i]]*boundfunc.sym[subcell_k].signy;
                shapefunc.cutcell[n].FS_hyy[i] = subcell.FS_hyy[boundfunc.sym[subcell_k].y[i]]*boundfunc.sym[subcell_k].signy;
                shapefunc.cutcell[n].FS_hyx[i] = subcell.FS_hyx[boundfunc.sym[subcell_k].y[i]]*boundfunc.sym[subcell_k].signx;
            }
            /*shapefunc.info<0 : completely outside or inside*/
            /*shapefunc.info>=0: indes of the cut-cell*/

            shapefunc.info[cell_i][subcell_k]=*ptr_bd;
            *ptr_bd += 1;
        }
        else
            shapefunc.info[cell_i][subcell_k]=type;

    }
    return Mcell;
}


boundfunc_cell compute_matricesSubcell(double bound_value[4], char bd_type, int* ptr_type)
{
    /* We define the variable CORNERS containing the x-coordinates,
    % y-coordinates and signed-function value for each of the corners.
    % The coordinates are with respect to the left-lower corner origin, and
    %normalized by the cell size dim.h*/
    int     i,j;
    int     bound_sign[4],side,side_next,ncorner1 =4,ncorner2;
    int     ind1[4] = {0,1,2,3}, ind2[4]={0,1,3,4};
    double  corners[4][8] = {{0,0.5,0.5,0},{0,0,0.5,0.5},{bound_value[0],bound_value[1],bound_value[2],bound_value[3]},{0,1,2,3}};
    double  new_corner[2],start_coord[2],end_coord[2];
    struct boundfunc_edge      edge;
    struct boundfunc_edgeNS    edgeNS;
    boundfunc_cell             subcell={0};

    for(i=0;i<4;i++)
        bound_sign[i]=SIGN(1,bound_value[i]);

    /*not a single corner in*/
    if( (bound_sign[0]!= -1) && (bound_sign[1]!= -1) && (bound_sign[2]!= -1) && (bound_sign[3]!= -1) )
    {
        *ptr_type = -2;
        return subcell;
    }

    /* no more than one corner on the border*/
    if( (bound_sign[0]+bound_sign[1]+bound_sign[2]+bound_sign[3]) <= -3 )
    {
        *ptr_type = -1;
        subcell = boundfunc_subcell0;
        return subcell;
    }

    *ptr_type = 1;

    // Compute the boundary corners
    for(side=3;side>-1;side--)
    {
        side_next = (side+1)%4;
        if( corners[2][side]*corners[2][side_next]<0)
        {
            ncorner1 +=1;

            for(i=0;i<2;i++)
                new_corner[i] = corners[i][side]+(corners[i][side_next]-corners[i][side])
                            *corners[2][side]/(corners[2][side]-corners[2][side_next]);

            for(i=0; i<4;i++)
                for(j=ncorner1-1; j>side+1;j--)
                    corners[i][j] = corners[i][j-1];

            corners[0][side+1] = new_corner[0];
            corners[1][side+1] = new_corner[1];
            corners[2][side+1] = 0;
            corners[3][side+1] = side;
        }
    }
    ncorner2 = ncorner1;

    // Delete the external corners
    for(side=ncorner1-1;side>-1;side--)
    {
        if (corners[2][side] >0)
        {
            ncorner2 -= 1;
            for(j=side; j<ncorner2;j++)
                for(i=0; i<4; i++)
                    corners[i][j] = corners[i][j+1];
        }
    }

    // ncorner2 - corner number after update

    /*Compute for boundary subcells contribution, edge-by-edge*/
    // Add the contribution of each edge into the line integral
    for(side=0; side<ncorner2;side++)
    {
        side_next = (side+1)%ncorner2;
        start_coord[0] = corners[0][side];
        start_coord[1] = corners[1][side];
        end_coord[0]   = corners[0][side_next];
        end_coord[1]   = corners[1][side_next];

        if ( corners[1][side] != corners[1][side_next])
        {
            edge = compute_matricesEdge( start_coord, end_coord);

            for(i=0; i<4; i++)
            {
                for(j=0; j<4; j++)
                {
                    subcell.MS_xx[ind1[i]][ind1[j]]     += edge.MS_xx[i][j];
                    subcell.MS_yy[ind2[i]][ind2[j]]     += edge.MS_yy[i][j];

                    subcell.KS_hxxhxx[ind1[i]][ind1[j]] += edge.KS_hxxhxx[i][j];
                    subcell.KS_hyyhyy[ind2[i]][ind2[j]] += edge.KS_hyyhyy[i][j];
                    subcell.KS_hxyhxy[ind1[i]][ind1[j]] += edge.KS_hxyhxy[i][j];
                    subcell.KS_hyxhyx[ind2[i]][ind2[j]] += edge.KS_hyxhyx[i][j];

                    subcell.KS_hxxhyx[ind1[i]][ind2[j]] += edge.KS_hxxhyx[i][j];
                    subcell.KS_hxxhyy[ind1[i]][ind2[j]] += edge.KS_hxxhyy[i][j];
                    subcell.KS_hxyhyx[ind1[i]][ind2[j]] += edge.KS_hxyhyx[i][j];
                    subcell.KS_hxyhyy[ind1[i]][ind2[j]] += edge.KS_hxyhyy[i][j];

                    subcell.KS_hxxhxy[ind1[i]][ind1[j]] += edge.KS_hxxhxy[i][j];
                    subcell.KS_hyyhyx[ind2[i]][ind2[j]] += edge.KS_hyyhyx[i][j];
                }
                subcell.FS_hxx[ind1[i]] += edge.FS_hxx[i];
                subcell.FS_hxy[ind1[i]] += edge.FS_hxy[i];
                subcell.FS_hyy[ind2[i]] += edge.FS_hyy[i];
                subcell.FS_hyx[ind2[i]] += edge.FS_hyx[i];
            }
        }

        if ((corners[2][side]==0) && (corners[2][side_next]==0))
        {
            edgeNS = compute_matricesEdgeNS( start_coord, end_coord);

            /*if (bd_type == 'F')
            {
                subcell.coord = [ corners(1:2,side) , corners(1:2,side_next) ];
                subcell.TS0.xc( : , end+1 ) = - edge.NS.xc;
                subcell.TS0.yc( : , end+1 ) = - edge.NS.yc;
            }*/
            if (bd_type == 'N')
            {
                for(i=0; i<4; i++)
                {
                    subcell.NS_xc[ind1[i]][0] += edgeNS.NS_xc[i];
                    subcell.NS_yc[ind2[i]][0] += edgeNS.NS_yc[i];
                }
            }
            else if (bd_type == 'D')
            {
                for(i=0; i<4; i++)
                {
                    subcell.DS_xc[ind1[i]][0] -= edgeNS.NS_xc[i];
                    subcell.DS_yc[ind2[i]][0] -= edgeNS.NS_yc[i];
                }
            }
        }
    }
    return subcell;
}

struct boundfunc_edge compute_matricesEdge(double start_coord[2], double end_coord[2])
{
    struct boundfunc_edge  edge={0};
    double          alpha, beta,delta;
    int             i,j;

    alpha    = ( end_coord[0]-start_coord[0] ) / ( end_coord[1]-start_coord[1]);

    if (abs(alpha) <= 1)
        beta  = start_coord[0] - alpha*start_coord[1];
    else
        delta = start_coord[1] - start_coord[0]/alpha;

    //Calculation using alpha-formula
    if ( abs(alpha) <= 1)
    {
        for(i=0;i<4;i++)
        {
           for (j=0;j<4;j++)
          {
            //MS.xx coefficients
            edge.MS_xx[i][j] = boundFunct_MS_alpha(  i,  j,  end_coord[0],  end_coord[1], alpha)
                             - boundFunct_MS_alpha(  i,  j,start_coord[0],start_coord[1], alpha);
            //MS.yy coefficients
            edge.MS_yy[i][j] = boundFunct_MS_alpha(4+i,4+j,  end_coord[0],  end_coord[1], alpha)
                             - boundFunct_MS_alpha(4+i,4+j,start_coord[0],start_coord[1], alpha);

            //KS coefficients
            edge.KS_hxxhxx[i][j] = boundFunct_KS_xx_alpha(  i,  j,  start_coord[1],  end_coord[1], alpha, beta);
            edge.KS_hyyhyy[i][j] = boundFunct_KS_yy_alpha(4+i,4+j,  start_coord[1],  end_coord[1], alpha, beta);
            edge.KS_hxyhxy[i][j] = boundFunct_KS_yy_alpha(  i,  j,  start_coord[1],  end_coord[1], alpha, beta);
            edge.KS_hyxhyx[i][j] = boundFunct_KS_xx_alpha(4+i,4+j,  start_coord[1],  end_coord[1], alpha, beta);

            edge.KS_hxxhyx[i][j] = boundFunct_KS_xx_alpha(  i,4+j,  start_coord[1],  end_coord[1], alpha, beta);
            edge.KS_hxxhyy[i][j] = boundFunct_KS_xy_alpha(j+4,  i,  start_coord[1],  end_coord[1], alpha, beta);
            edge.KS_hxyhyx[i][j] = boundFunct_KS_xy_alpha(  i,j+4,  start_coord[1],  end_coord[1], alpha, beta);
            edge.KS_hxyhyy[i][j] = boundFunct_KS_yy_alpha( j+4, i,  start_coord[1],  end_coord[1], alpha, beta);

            edge.KS_hxxhxy[i][j] = boundFunct_KS_xy_alpha(   i,  j,  start_coord[1],  end_coord[1], alpha, beta);
            edge.KS_hyyhyx[i][j] = boundFunct_KS_xy_alpha( j+4,i+4,  start_coord[1],  end_coord[1], alpha, beta);
          }

         edge.FS_hxx[i] = boundFunct_FS_x_alpha(  i,   start_coord[1],  end_coord[1], alpha, beta);
         edge.FS_hxy[i] = boundFunct_FS_y_alpha(  i,   start_coord[1],  end_coord[1], alpha, beta);
         edge.FS_hyy[i] = boundFunct_FS_y_alpha(4+i,   start_coord[1],  end_coord[1], alpha ,beta);
         edge.FS_hyx[i] = boundFunct_FS_x_alpha(4+i,   start_coord[1],  end_coord[1], alpha ,beta);
        }
    }

    //Calculation using gamma-formula
    else
    {
        for(i=0;i<4;i++)
        {
           for (j=0;j<4;j++)
          {
            //MS.xx coefficients
            edge.MS_xx[i][j] = boundFunct_MS_gamma(  i,  j,  end_coord[0],  end_coord[1], 1/alpha)
                             - boundFunct_MS_gamma(  i,  j,start_coord[0],start_coord[1], 1/alpha);
            //MS.yy coefficients
            edge.MS_yy[i][j] = boundFunct_MS_gamma(4+i,4+j,  end_coord[0],  end_coord[1], 1/alpha)
                             - boundFunct_MS_gamma(4+i,4+j,start_coord[0],start_coord[1], 1/alpha);
            //KS coefficients
            edge.KS_hxxhxx[i][j] = boundFunct_KS_xx_gamma(  i,  j,  start_coord[0],  end_coord[0], 1/alpha, delta);
            edge.KS_hyyhyy[i][j] = boundFunct_KS_yy_gamma(4+i,4+j,  start_coord[0],  end_coord[0], 1/alpha, delta);
            edge.KS_hxyhxy[i][j] = boundFunct_KS_yy_gamma(  i,  j,  start_coord[0],  end_coord[0], 1/alpha, delta);
            edge.KS_hyxhyx[i][j] = boundFunct_KS_xx_gamma(4+i,4+j,  start_coord[0],  end_coord[0], 1/alpha, delta);

            edge.KS_hxxhyx[i][j] = boundFunct_KS_xx_gamma(  i,4+j,  start_coord[0],  end_coord[0], 1/alpha, delta);
            edge.KS_hxxhyy[i][j] = boundFunct_KS_xy_gamma(j+4,  i,  start_coord[0],  end_coord[0], 1/alpha, delta);
            edge.KS_hxyhyx[i][j] = boundFunct_KS_xy_gamma(  i,j+4,  start_coord[0],  end_coord[0], 1/alpha, delta);
            edge.KS_hxyhyy[i][j] = boundFunct_KS_yy_gamma( j+4, i,  start_coord[0],  end_coord[0], 1/alpha, delta);

            edge.KS_hxxhxy[i][j] = boundFunct_KS_xy_gamma(   i,  j,  start_coord[0],  end_coord[0], 1/alpha, delta);
            edge.KS_hyyhyx[i][j] = boundFunct_KS_xy_gamma( j+4,i+4,  start_coord[0],  end_coord[0], 1/alpha, delta);

          }
         edge.FS_hxx[i] = boundFunct_FS_x_gamma(  i,   start_coord[0],  end_coord[0], 1/alpha, delta);
         edge.FS_hxy[i] = boundFunct_FS_y_gamma(  i,   start_coord[0],  end_coord[0], 1/alpha, delta);
         edge.FS_hyy[i] = boundFunct_FS_y_gamma(4+i,   start_coord[0],  end_coord[0], 1/alpha ,delta);
         edge.FS_hyx[i] = boundFunct_FS_x_gamma(4+i,   start_coord[0],  end_coord[0], 1/alpha ,delta);
        }

    }
    return edge;
}

struct boundfunc_edgeNS compute_matricesEdgeNS(double start_coord[2], double end_coord[2])
{

    struct boundfunc_edgeNS  edge={0};
    double            alpha;
    int               i;

    alpha    = ( end_coord[0]-start_coord[0] ) / ( end_coord[1]-start_coord[1]);
        //Calculation using alpha-formula
    if ( abs(alpha) <= 1)
    {
        for(i=0;i<4;i++)
        {
         //Calculation of the NS coefficients
         edge.NS_xc[i] = SIGN( 1, end_coord[1]-start_coord[1] ) *
                           ( boundFunct_NS_alpha(   i,  end_coord[0],  end_coord[1], alpha)
                            -boundFunct_NS_alpha(   i,start_coord[0],start_coord[1], alpha) );
         edge.NS_yc[i] = SIGN( 1, end_coord[1]-start_coord[1] ) *
                           ( boundFunct_NS_alpha( 4+i,  end_coord[0],  end_coord[1], alpha)
                            -boundFunct_NS_alpha( 4+i,start_coord[0],start_coord[1], alpha) );
        }
    }
    else
    {
        for(i=0;i<4;i++)
        {
         //Calculation of the NS coefficients
         edge.NS_xc[i] = SIGN( 1,end_coord[0]-start_coord[0] ) *
                           ( boundFunct_NS_gamma(   i,  end_coord[0],  end_coord[1], 1/alpha)
                            -boundFunct_NS_gamma(   i,start_coord[0],start_coord[1], 1/alpha) );

         edge.NS_yc[i] = SIGN( 1,end_coord[0]-start_coord[0] ) *
                           ( boundFunct_NS_gamma( 4+i,  end_coord[0],  end_coord[1], 1/alpha)
                            -boundFunct_NS_gamma( 4+i,start_coord[0],start_coord[1], 1/alpha) );
        }
    }

    return edge;
}

double boundFunct_MS_alpha(int k, int l,double x,double y,double alpha){

     double z;
     double y3 = y*y*y, y4 = y3*y, y5= y4*y, y6 = y5*y;

     z = dx2* boundfunc.coeff[0][k]*boundfunc.coeff[0][l]*(
         (2*(x*x*x)/6       + (boundfunc.coeff[1][k]+boundfunc.coeff[1][l])*(x*x)/2  + boundfunc.coeff[1][k]*boundfunc.coeff[1][l]*x)
         * (2*(y3)/6        + (boundfunc.coeff[2][k]+boundfunc.coeff[2][l])*(y*y)/2  + boundfunc.coeff[2][k]*boundfunc.coeff[2][l]*y)
         + (-alpha)*(2*(x*x)/2 + (boundfunc.coeff[1][k]+boundfunc.coeff[1][l])* x    + boundfunc.coeff[1][k]*boundfunc.coeff[1][l]  )
         * (2*(y4)/24 + (boundfunc.coeff[2][k]+boundfunc.coeff[2][l])*(y3)/6         + boundfunc.coeff[2][k]*boundfunc.coeff[2][l]*(y*y)/2 )
         + alpha*alpha*(2* x                                                         + (boundfunc.coeff[1][k]+boundfunc.coeff[1][l])                                                        )
         * (2*(y5)/120+ (boundfunc.coeff[2][k]+boundfunc.coeff[2][l])*(y4)/24        + boundfunc.coeff[2][k]*boundfunc.coeff[2][l]*(y3)/6 )
         - alpha*alpha*alpha*2
         * (2*(y6)/720+ (boundfunc.coeff[2][k]+boundfunc.coeff[2][l])*(y5)/120       + boundfunc.coeff[2][k]*boundfunc.coeff[2][l]*(y4)/24) );


    return z;
}

double boundFunct_MS_gamma(int k, int l,double x,double y,double gamma){

     double z;
     double x3 = x*x*x, x4 = x3*x, x5= x4*x, x6 = x5*x;

     z = dx2 * boundfunc.coeff[0][k] * boundfunc.coeff[0][l]* gamma
         *(   (2*(x4)/24 + (boundfunc.coeff[1][k]+boundfunc.coeff[1][l])*(x3)/6  + boundfunc.coeff[1][k]*boundfunc.coeff[1][l]*(x*x)/2 )
            * (2*(y*y)/2 + (boundfunc.coeff[2][k]+boundfunc.coeff[2][l])*y + boundfunc.coeff[2][k]*boundfunc.coeff[2][l])
            + (-gamma)  *(2*(x5)/120+ (boundfunc.coeff[1][k]+boundfunc.coeff[1][l])*(x4)/24 + boundfunc.coeff[1][k]*boundfunc.coeff[1][l]*(x3)/6 )
            * (2* y      + (boundfunc.coeff[2][k]+boundfunc.coeff[2][l])                                                )
            + gamma*gamma*(2*(x6)/720+ (boundfunc.coeff[1][k]+boundfunc.coeff[1][l])*(x5)/120+ boundfunc.coeff[1][k]*boundfunc.coeff[1][l]*(x4)/24)
            * 2);
     return z;
}

double boundFunct_KS_xx_alpha(int k, int l,double y1,double y2,double alpha, double beta){

    double z;
    double y13 = y1*y1*y1, y14 = y13*y1, y23= y2*y2*y2, y24 = y23*y2;

    z = boundfunc.coeff[0][k] * boundfunc.coeff[0][l]
        *( alpha *(y24-y14)/4
         + (beta + alpha* (boundfunc.coeff[2][k]+boundfunc.coeff[2][l])) *(y23-y13)/3
         + ( boundfunc.coeff[2][k]*boundfunc.coeff[2][l]*alpha
         + beta*( boundfunc.coeff[2][k] + boundfunc.coeff[2][l]))        *(y2*y2-y1*y1)/2
         + beta*boundfunc.coeff[2][k] * boundfunc.coeff[2][l]            *(y2-y1) );

    return z;
}

double boundFunct_KS_xx_gamma(int k, int l,double x1,double x2,double gamma, double delta){

    double z;
    double x13 = x1*x1*x1, x14 = x13*x1, x23= x2*x2*x2, x24 = x23*x2;

    z = boundfunc.coeff[0][k] * boundfunc.coeff[0][l]
        *( gamma*gamma*gamma                                                         *(x24-x14)/4
           + gamma*gamma*( boundfunc.coeff[2][k] + boundfunc.coeff[2][l] + 2*delta ) *(x23-x13)/3
           + gamma*(delta + boundfunc.coeff[2][k] )* (delta + boundfunc.coeff[2][l]) *(x2*x2-x1*x1)/2 );

    return z;
}

double boundFunct_KS_yy_alpha(int k, int l,double y1,double y2,double alpha, double beta){

    double z;
    double y13 = y1*y1*y1, y14 = y13*y1, y23= y2*y2*y2, y24 = y23*y2;

    z = boundfunc.coeff[0][k] * boundfunc.coeff[0][l]
        *(   alpha*alpha*alpha                                                         *(y24-y14)/12
           + alpha*alpha*  ( beta + (boundfunc.coeff[1][k]+boundfunc.coeff[1][l])/2 )  *(y23-y13)/3
           + alpha*(beta + boundfunc.coeff[1][k])*(beta + boundfunc.coeff[1][l])       *(y2*y2-y1*y1)/2
           + (beta*beta*beta/3 + (boundfunc.coeff[1][k]+boundfunc.coeff[1][l])/2*beta*beta + boundfunc.coeff[1][k]*boundfunc.coeff[1][l]*beta ) *(y2-y1) );

    return z;
}

double boundFunct_KS_yy_gamma(int k, int l,double x1,double x2,double gamma, double delta){

    double z;
    double x13 = x1*x1*x1, x14 = x13*x1, x23= x2*x2*x2, x24 = x23*x2;

    z = boundfunc.coeff[0][k] * boundfunc.coeff[0][l]* gamma
        *( (x24-x14)/12
            + (boundfunc.coeff[1][k] + boundfunc.coeff[1][l]) *(x23-x13)/6
            +  boundfunc.coeff[1][k]*boundfunc.coeff[1][l]    *(x2*x2-x1*x1)/2 );

    return z;
}

double boundFunct_KS_xy_alpha(int k, int l,double y1,double y2,double alpha, double beta){

    double z;
    double y13 = y1*y1*y1, y14 = y13*y1, y23= y2*y2*y2, y24 = y23*y2;

    z = boundfunc.coeff[0][k] * boundfunc.coeff[0][l]/2
        *( alpha*alpha                                                                 *(y24-y14)/4
           +  alpha* (alpha*boundfunc.coeff[2][l] +2*(beta + boundfunc.coeff[1][k]) )  *(y23-y13)/3
           + (2*boundfunc.coeff[2][l]*alpha+ beta + boundfunc.coeff[1][k])*(beta + boundfunc.coeff[1][k])      *(y2*y2-y1*y1)/2
           + boundfunc.coeff[2][l]*(beta + boundfunc.coeff[1][k])*(beta + boundfunc.coeff[1][k])               *(y2-y1)  );

    return z;
}

double boundFunct_KS_xy_gamma(int k, int l,double x1,double x2,double gamma, double delta){

    double z;
    double x13 = x1*x1*x1, x14 = x13*x1, x23= x2*x2*x2, x24 = x23*x2;

    z = boundfunc.coeff[0][k] * boundfunc.coeff[0][l]/2
        *( gamma*gamma                                                                *(x24-x14)/4
           + gamma* (2*gamma*boundfunc.coeff[1][k] + delta + boundfunc.coeff[2][l] )  *(x23-x13)/3
           + (2*boundfunc.coeff[1][k]*gamma*( delta +boundfunc.coeff[2][l]) + boundfunc.coeff[1][k]*boundfunc.coeff[1][k]*gamma*gamma)    *(x2*x2-x1*x1)/2
           + boundfunc.coeff[1][k]*boundfunc.coeff[1][k]*gamma*(delta + boundfunc.coeff[2][l])                                            *(x2-x1) );

    return z;
}

double boundFunct_FS_x_alpha(int k,double y1,double y2,double alpha, double beta){

    double z;
    double y12 = y1*y1, y13 = y12*y1, y22 = y2*y2, y23= y22*y2;

    z = dx * boundfunc.coeff[0][k]
        *( alpha                                    *(y23-y13)/3
           + (alpha*boundfunc.coeff[2][k] + beta )  *(y22-y12)/2
           + beta*boundfunc.coeff[2][k]             *(y2-y1) );

    return z;
}

double boundFunct_FS_x_gamma(int k,double x1,double x2,double gamma, double delta){

    double z;
    double x12 = x1*x1, x13 = x12*x1, x22= x2*x2, x23 = x22*x2;

    z = dx * boundfunc.coeff[0][k]
        *( gamma*gamma                          *(x23-x13)/3
        + gamma*(delta + boundfunc.coeff[2][k]) *(x22-x12)/2 );

    return z;
}

double boundFunct_FS_y_alpha(int k,double y1,double y2,double alpha, double beta){

    double z;
    double y12 = y1*y1, y13 = y12*y1, y22 = y2*y2, y23= y22*y2;

    z = dx * boundfunc.coeff[0][k]/2
        *( alpha*alpha                              *(y23-y13)/3
           + alpha*(boundfunc.coeff[1][k] + beta )  *(y22-y12)
           + (boundfunc.coeff[1][k] + beta )*(boundfunc.coeff[1][k] + beta )      *(y2-y1) );

    return z;
}

double boundFunct_FS_y_gamma(int k,double x1,double x2,double gamma, double delta){

    double z;

    z = dx * boundfunc.coeff[0][k]* gamma
        *( ( (x2+boundfunc.coeff[1][k])*(x2+boundfunc.coeff[1][k])*(x2+boundfunc.coeff[1][k])
            -(x1+boundfunc.coeff[1][k])*(x1+boundfunc.coeff[1][k])*(x1+boundfunc.coeff[1][k]) )/6 );

    return z;
}

double boundFunct_NS_alpha(int l, double x, double y, double alpha){

    double z;
    z = dx * boundfunc.coeff[0][l] * sqrt(alpha*alpha+1)
        *(   (x + boundfunc.coeff[1][l] )  *  ( (y*y)/2 + boundfunc.coeff[2][l]* y)
             + (-alpha)  *( (y*y*y)/6 + boundfunc.coeff[2][l]*(y*y)/2 ) );

    return z;
}

double boundFunct_NS_gamma(int l, double x, double y, double gamma){

    double z;
    z = dx * boundfunc.coeff[0][l] * sqrt(gamma*gamma+1)
        *(  ( y + boundfunc.coeff[2][l] )  *  ( (x*x)/2 + boundfunc.coeff[1][l]* x )
              + (-gamma)  *( (x*x*x)/6 + boundfunc.coeff[1][l]*(x*x)/2 ) );
    return z;
}

void reduce_system(Grid_S* ptr_g, Index_S* ptr_i, Mat* ptr_MS, Mat* ptr_KLS, Mat* ptr_NS, Mat* ptr_DS)
// reduce_system.m
{
    /*% REDUCE_SYSTEM reduces the initial system of equation  A * t[u, v, p] = F by removing:
    - irrelevant unknowns in u, v or p which have no-influence on the domain
    - by getting rid of the arbitrary solid-body-motion*/

    int         glon_xi =0, glon_cell_neu =0,glon_cell_dir =0; // global number
    int         n_xix =0, n_xix_gho =0, n_cellx_neu =0,n_cellx_dir =0;       // local number of cells
    int         n_xiy =0, n_xiy_gho =0, n_celly_neu =0,n_celly_dir =0;
    int         ns_xixy =0, ns_cellxy_neu =0, ns_cellxy_dir =0;
    int         i,j,temp;
    int         *index_keep_xi;            // g -> g'
    int         *index_keep_stag_Neumann;  // c -> c'
    int         *index_keep_stag_Dirichlet;// c -> c'
    Vec          MS_diag_vec, KLS_diag_vec,MS_diag_locv, KLS_diag_locv;
    PetscScalar *MS_diag,*KLS_diag;
    PetscReal   *NS_norm, *DS_norm;
    PetscInt	*is_xi, *is_dir, *is_neu;
    IS 			ixy;
    VecScatter  scatter_MS,scatter_KS;
    
    VecCreate(PETSC_COMM_WORLD,&MS_diag_vec);
    VecSetSizes(MS_diag_vec,PETSC_DECIDE,ptr_i->xi_gloN);
    VecSetFromOptions(MS_diag_vec);
    VecDuplicate(MS_diag_vec,&KLS_diag_vec);

    MatGetDiagonal(*ptr_MS,MS_diag_vec);  // can only be applied on the parallel mat and vec
    MatGetDiagonal(*ptr_KLS,KLS_diag_vec);

    VecCreateSeq(PETSC_COMM_SELF,ptr_i->xi_gloN,&MS_diag_locv); // create a sequential loc vec
    VecDuplicate(MS_diag_locv,&KLS_diag_locv);

    ISCreateStride(PETSC_COMM_SELF,ptr_i->xi_gloN,0,1,&ixy);
    VecScatterCreate(MS_diag_vec,ixy,MS_diag_locv,ixy,&scatter_MS);
    VecScatterBegin(scatter_MS,MS_diag_vec,MS_diag_locv,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(scatter_MS,MS_diag_vec,MS_diag_locv,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterCreate(KLS_diag_vec,ixy,KLS_diag_locv,ixy,&scatter_KS);
    VecScatterBegin(scatter_KS,KLS_diag_vec,KLS_diag_locv,INSERT_VALUES,SCATTER_FORWARD);
    VecScatterEnd(scatter_KS,KLS_diag_vec,KLS_diag_locv,INSERT_VALUES,SCATTER_FORWARD);

    VecGetArray(MS_diag_locv, &MS_diag);
    VecGetArray(KLS_diag_locv,&KLS_diag); // each processor will get all the info
    
    PetscMalloc1(ptr_i->xi_gloNcell_Dirichlet,&DS_norm);
    PetscMalloc1(ptr_i->xi_gloNcell_Neumann  ,&NS_norm);

    MatGetColumnNorms(*ptr_DS ,NORM_INFINITY,DS_norm);
    MatGetColumnNorms(*ptr_NS ,NORM_INFINITY,NS_norm);

    index_keep_xi              = (int*) calloc(sizeof(int),ptr_i->xi_gloN);
    index_keep_stag_Neumann    = (int*) calloc(sizeof(int),ptr_i->xi_gloNcell_Neumann);
    index_keep_stag_Dirichlet  = (int*) calloc(sizeof(int),ptr_i->xi_gloNcell_Dirichlet );

    //VecView(MS_diag_vec,PETSC_VIEWER_STDOUT_SELF); //
    //VecView(KLS_diag_vec,PETSC_VIEWER_STDOUT_SELF); //

    // ok to here
    // =================  xix, xiy ======================
    //List of the indices of variables which are not outside of the domain
    for(i=0; i<ptr_i->xi_gloN;i++)
    {
        if( ( MS_diag[i] != 0 ) || ( KLS_diag[i] != 0 ))
        {
            index_keep_xi[i] = glon_xi;
            glon_xi +=1;
        }
        else
        	index_keep_xi[i] = -1;
    }

    temp = ptr_i->xix_N + ptr_i->xiy_N;
    MPI_Scan(&temp,&ns_xixy,1,MPIU_INT,MPI_SUM,PETSC_COMM_WORLD); // get the starting index of xi x,y on each processor
    ns_xixy -= temp;

    for(i= 0 ; i< ptr_i->xix_N ;i++)
    {
    	if(index_keep_xi[ptr_i->xix.l2g[i]] != -1)
    	{
    		ptr_i->xix.l2g[n_xix] = ptr_i->xix.l2g[i];
    		ptr_i->xix.l2G[n_xix] = ptr_i->xix.l2G[i];
    		n_xix ++;  // new xix_N
    	}
    }

    for(i= ptr_i->xix_N ; i< ptr_i->xix_N + ptr_i->xix_ghoN ;i++)
    {
        if(index_keep_xi[ptr_i->xix.l2g[i]] != -1)
        {
        	ptr_i->xix.l2g[n_xix + n_xix_gho] = ptr_i->xix.l2g[i];
        	ptr_i->xix.l2G[n_xix + n_xix_gho] = ptr_i->xix.l2G[i];
        	n_xix_gho ++;  // new xix_N
        }
    }

    for(i= 0 ; i< ptr_i->xiy_N ;i++)
    {
    	if(index_keep_xi[ptr_i->xiy.l2g[i]] != -1)
    	{
    		ptr_i->xiy.l2g[n_xiy] = ptr_i->xiy.l2g[i];
    		ptr_i->xiy.l2G[n_xiy] = ptr_i->xiy.l2G[i];
    		n_xiy ++;  // new xix_N
    	}
    }

    for(i= ptr_i->xiy_N ; i< ptr_i->xiy_N + ptr_i->xiy_ghoN ;i++)
    {
        if(index_keep_xi[ptr_i->xiy.l2g[i]] != -1)
        {
        	ptr_i->xiy.l2g[n_xiy + n_xiy_gho] = ptr_i->xiy.l2g[i];
        	ptr_i->xiy.l2G[n_xiy + n_xiy_gho] = ptr_i->xiy.l2G[i];
        	n_xiy_gho ++;  // new xix_N
        }
    }
    /*
    //temp = n_xix + n_xiy;
    //MPI_Scan(&temp,&ns_xixy_new,1,MPIU_INT,MPI_SUM,PETSC_COMM_WORLD); // get the starting index of xi x,y on each processor
    //ns_xixy_new -= temp;

    ptr_i->xix.G2g_before = (short int *) malloc(sizeof(short int)*ptr_g->N);
    ptr_i->xiy.G2g_before = (short int *) malloc(sizeof(short int)*ptr_g->N);

    for(i= 0 ; i< ptr_g->N ;i++)
    {
        ptr_i->xix.G2g_before[i] = ptr_i->xix.G2g[i];
        ptr_i->xiy.G2g_before[i] = ptr_i->xiy.G2g[i];
        ptr_i->xix.G2g[i] = index_keep_xi[ptr_i->xix.G2g[i]];
        ptr_i->xiy.G2g[i] = index_keep_xi[ptr_i->xiy.G2g[i]];
    }

    PetscMalloc1(n_xix + n_xiy,&is_xi);   // local index
    temp = 0;
    for(i= ns_xixy ; i< ns_xixy + ptr_i->xix_N + ptr_i->xiy_N ;i++)
    {
    	if(index_keep_xi[i] != -1)
    	{
    		is_xi[temp] = i;
    		temp++;
    	}
    }

    //List of the indices of cells acting on Dirichlet boundary
    //========================= dir ============================

    for(j=0;j< ptr_i->xi_gloNcell_Dirichlet; j++)
    {
        if (DS_norm[j] != 0)
        {
             index_keep_stag_Dirichlet[j] = glon_cell_dir;
             glon_cell_dir +=1;
        }
        else
        	index_keep_stag_Dirichlet[j] = -1;
    }

    temp = ptr_i->xix_Ncell_Dirichlet + ptr_i->xiy_Ncell_Dirichlet;
    MPI_Scan(&temp,&ns_cellxy_dir,1,MPIU_INT,MPI_SUM,PETSC_COMM_WORLD); // get the starting index of xi x,y on each processor
    ns_cellxy_dir -= temp;


    for(i= 0 ; i< ptr_i->xix_Ncell_Dirichlet ;i++)
    {
    	if(index_keep_stag_Dirichlet[ ptr_i->xix.C2c_dirichlet[ptr_i->xix.c2C_dirichlet[i]] ] != -1)
    	{
    		ptr_i->xix.c2C_dirichlet[n_cellx_dir] = ptr_i->xix.c2C_dirichlet[i];
    		n_cellx_dir ++;  // new xix_N
    	}
    }

    for(i= 0 ; i< ptr_i->xiy_Ncell_Dirichlet ;i++)
    {
    	if(index_keep_stag_Dirichlet[ ptr_i->xiy.C2c_dirichlet[ptr_i->xiy.c2C_dirichlet[i]] ] != -1)
    	{
    		ptr_i->xiy.c2C_dirichlet[n_celly_dir] = ptr_i->xiy.c2C_dirichlet[i];
    		n_celly_dir ++;  // new xix_N
    	}
    }

    for(i= 0 ; i< ptr_g->N ;i++)
    {
    	ptr_i->xix.C2c_dirichlet[i] = index_keep_stag_Dirichlet[ptr_i->xix.C2c_dirichlet[i]];
        ptr_i->xiy.C2c_dirichlet[i] = index_keep_stag_Dirichlet[ptr_i->xiy.C2c_dirichlet[i]];
    }

    PetscMalloc1(n_cellx_dir + n_celly_dir,&is_dir);

    temp = 0;
    for(i= ns_cellxy_dir ; i< ns_cellxy_dir + ptr_i->xix_Ncell_Dirichlet + ptr_i->xiy_Ncell_Dirichlet ;i++)
    {
    	if(index_keep_stag_Dirichlet[i] != -1)
    	{
    		is_dir[temp] = i;
    		temp++;
    	}
    }

    //List of the indices of cells acting on Neumann boundary
    // ========================= neu ===========================

    for(j=0;j< ptr_i->xix_gloNcell_Neumann; j++)
     {
         if (NS_norm[j] != 0)
         {
        	 index_keep_stag_Neumann[j] = glon_cell_neu;
             glon_cell_neu +=1;
         }
         else
        	 index_keep_stag_Neumann[j] = -1;
     }

     temp = ptr_i->xix_Ncell_Neumann + ptr_i->xiy_Ncell_Neumann;
     MPI_Scan(&temp,&ns_cellxy_dir,1,MPIU_INT,MPI_SUM,PETSC_COMM_WORLD); // get the starting index of xi x,y on each processor
     ns_cellxy_neu -= temp;


     for(i= 0 ; i< ptr_i->xix_Ncell_Neumann ;i++)
     {
     	if(index_keep_stag_Neumann[ ptr_i->xix.C2c_neumann[ptr_i->xix.c2C_neumann[i]] ] != -1)
     	{
     		ptr_i->xix.c2C_neumann[n_cellx_neu] = ptr_i->xix.c2C_neumann[i];
     		n_cellx_neu ++;  // new xix_N
     	}
     }

     for(i= 0 ; i< ptr_i->xiy_Ncell_Neumann ;i++)
     {
     	if(index_keep_stag_Neumann[ ptr_i->xiy.C2c_neumann[ptr_i->xiy.c2C_neumann[i]] ] != -1)
     	{
     		ptr_i->xiy.c2C_neumann[n_celly_neu] = ptr_i->xiy.c2C_neumann[i];
     		n_celly_neu ++;  // new xix_N
     	}
     }

     for(i= 0 ; i< ptr_g->N ;i++)
     {
     	ptr_i->xix.C2c_neumann[i] = index_keep_stag_Neumann[ptr_i->xix.C2c_neumann[i]];
        ptr_i->xiy.C2c_neumann[i] = index_keep_stag_Neumann[ptr_i->xiy.C2c_neumann[i]];
     }

     PetscMalloc1(n_cellx_neu + n_celly_neu,&is_neu);

     temp = 0;
     for(i= ns_cellxy_neu ; i< ns_cellxy_neu + ptr_i->xix_Ncell_Neumann + ptr_i->xiy_Ncell_Neumann ;i++)
     {
     	if(index_keep_stag_Neumann[i] != -1)
     	{
     		is_neu[temp] = i;
     		temp++;
     	}
     }

    // --------Update the l2g, l2G and c2C arrays ------------

    //ptr_i->xix.keep       = (short int *) malloc(sizeof(short int)*n_xix);
    //ptr_i->xiy.keep       = (short int *) malloc(sizeof(short int)*n_xiy);


	//for(i=0;i<n_xix;i++)
    //{
    //    ptr_i->xix.keep[i] = index_keep_xix[i];
    //}

	//for(i=0;i<n_xiy;i++)
    //{
    //    ptr_i->xiy.keep[i] = index_keep_xiy[i];
    //}

    reduce_vector(&(ptr_i->xix.l2G),ptr_i->xix_N + ptr_i->xix_ghoN,n_xix + n_xix_gho);
    reduce_vector(&(ptr_i->xix.l2g),ptr_i->xix_N + ptr_i->xix_ghoN,n_xix + n_xix_gho);

    reduce_vector(&(ptr_i->xiy.l2G),ptr_i->xiy_N + ptr_i->xiy_ghoN,n_xiy + n_xiy_gho);
    reduce_vector(&(ptr_i->xiy.l2g),ptr_i->xiy_N + ptr_i->xiy_ghoN,n_xiy + n_xiy_gho);

    reduce_vector(&(ptr_i->xix.c2C_neumann),ptr_i->xix_Ncell_Neumann,n_cellx_neu);
    reduce_vector(&(ptr_i->xiy.c2C_neumann),ptr_i->xiy_Ncell_Neumann,n_celly_neu);

    reduce_vector(&(ptr_i->xix.c2C_dirichlet),ptr_i->xix_Ncell_Dirichlet,n_cellx_dir);
    reduce_vector(&(ptr_i->xiy.c2C_dirichlet),ptr_i->xiy_Ncell_Dirichlet,n_celly_dir);


	MPI_Reduce(&n_xix, &ptr_i->xix_gloN, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ptr_i->xix_gloN, 1, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Reduce(&n_xiy, &ptr_i->xi_gloN, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ptr_i->xi_gloN, 1, MPI_INT, 0, MPI_COMM_WORLD);
	ptr_i->xi_gloN += ptr_i->xix_gloN;

	MPI_Reduce(&n_cellx_neu, &ptr_i->xix_gloNcell_Neumann, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ptr_i->xix_gloNcell_Neumann, 1, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Reduce(&n_celly_neu, &ptr_i->xi_gloNcell_Neumann, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ptr_i->xi_gloNcell_Neumann, 1, MPI_INT, 0, MPI_COMM_WORLD);
	ptr_i->xi_gloNcell_Neumann += ptr_i->xix_gloNcell_Neumann;

	MPI_Reduce(&n_cellx_dir, &ptr_i->xix_gloNcell_Dirichlet, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ptr_i->xix_gloNcell_Dirichlet, 1, MPI_INT, 0, MPI_COMM_WORLD);

	MPI_Reduce(&n_celly_dir, &ptr_i->xi_gloNcell_Dirichlet, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);
	MPI_Bcast(&ptr_i->xi_gloNcell_Dirichlet, 1, MPI_INT, 0, MPI_COMM_WORLD);
	ptr_i->xi_gloNcell_Dirichlet += ptr_i->xix_gloNcell_Dirichlet;

	ptr_i->xi_gloN_before = ptr_i->xi_gloN;
	ptr_i->xix_gloN_before = ptr_i->xix_gloN;

    ptr_i->xix_N_before = ptr_i->xix_N;
    ptr_i->xiy_N_before = ptr_i->xiy_N;
    ptr_i->xix_N = n_xix;
    ptr_i->xiy_N = n_xiy;
    ptr_i->xix_Ncell_Neumann = n_cellx_neu;
    ptr_i->xiy_Ncell_Neumann = n_celly_neu;
    ptr_i->xix_Ncell_Dirichlet = n_cellx_dir;
    ptr_i->xiy_Ncell_Dirichlet = n_celly_dir;

    // to here
    ISCreateGeneral(MPI_COMM_WORLD,n_xix+n_xiy,is_xi,PETSC_COPY_VALUES,&ptr_i->is_xi);
    ISCreateGeneral(MPI_COMM_WORLD,n_cellx_neu + n_celly_neu,is_neu,PETSC_COPY_VALUES,&ptr_i->is_neu);
    ISCreateGeneral(MPI_COMM_WORLD,n_cellx_dir + n_celly_dir,is_dir,PETSC_COPY_VALUES,&ptr_i->is_dir);

    free(index_keep_xi );
    free(index_keep_stag_Neumann );
    free(index_keep_stag_Dirichlet );

    VecRestoreArray(MS_diag_vec, &MS_diag);
    VecRestoreArray(KLS_diag_vec,&KLS_diag);
    PetscFree(MS_diag); 
    PetscFree(KLS_diag);
    PetscFree(NS_norm); 
    PetscFree(DS_norm);
    PetscFree(is_xi);
    PetscFree(is_dir);
    PetscFree(is_neu);
    VecDestroy(&MS_diag_vec);
    VecDestroy(&KLS_diag_vec);
    */
    VecScatterDestroy(&scatter_MS);
    VecScatterDestroy(&scatter_KS);
}

// shrink the vector size from length n_b to n_a
void reduce_vector(short int** a, int n_b, int n_a)
{
    int i;
    short int* temp1, temp2;

    temp1 = (short int *) malloc(sizeof(short int)*n_a);
    for(i=0;i<n_a;i++)
        temp1[i] = (*a)[i];

    temp2 = *a;
    *a = temp1;

    free(temp2);

}
