/* parameters that is set by the users*/
#include "user_param.h"

#include <string.h>
#include <stdio.h>
#include <math.h>
#include "nrutil.h"

void set_grid 		(Field_S* s);
void set_constraint	(Field_S* s,char bc);
void set_processor	(Field_S* s, AppCtx* ptr_u);

void user_param (Field_S* s, AppCtx* ptr_u)
{

	int		solid_resolution	= 10;
	double	domain[4]			={-1.0,1.0,-1.0,1.0};
	double  char_length			=1.0-1e-15;

	double  density				=1.0;
	double	nu					=0.25;
	double	E					=2.5;
	char	shape[]				="Disk"; //"Beam","Disk","Elliptic"

	char	boundarycondition	='N'; // N: all Neumann, F:all FSI, D: all Dirichlet
									  // M: all mixed, P: Neumann+Dirichlet
    double  solid_threshold     =0.1;
    double  damping[2]			= {0.0, 0.0};

	double	timestep			=1e-2;
	double	finaltime			=1e-2;
	double	Newmark_theta		=1.0/4;
	int		plot_per_steps		=1000;

	char	initial[]			="specified"; // rest, specified

	int		i,j;

	/* setting grid information*/
	s->char_lengh	= char_length;
	s->domain_x[0]	= domain[0];
	s->domain_x[1]	= domain[1];
	s->domain_y[0]	= domain[2];
	s->domain_y[1]	= domain[3];

	s->domain_size[0] = domain[1]-domain[0];
	s->domain_size[1] = domain[3]-domain[2];

	s->res = solid_resolution;

	if (s->domain_size[0] > s->domain_size[1])
		s->dx = s->domain_size[1]/solid_resolution;
	else
		s->dx = s->domain_size[0]/solid_resolution;
	set_grid(s);

	/* setting structure information*/
	s->param.length			= char_length;
	s->param.rho				= density;
	s->param.lambda			= E*nu/(1.0+nu)/(1.0-2.0*nu);
	s->param.mu				= E/2.0/(1.0+nu);
	s->param.threshold		= solid_threshold;
	strcpy( s->param.type, shape);
	strcpy( s->param.initial, initial);
	s->param.boundcond 		= boundarycondition;
	s->param.damping[0] 		= damping[0];
	s->param.damping[1] 		= damping[1];
	s->param.boundary_value 	= dmatrix(0,s->Nx-1,0,s->Ny-1);
	s->param.boundary_sign  	= imatrix(0,s->Nx-1,0,s->Ny-1);

	set_sdf(s, 0,s->Nx, 0, s->Ny, s->x_grid, s->y_grid, s->param.boundary_value);				// set sdf

	for(i=0; i<s->Nx; i++)
		for(j=0; j<s->Ny; j++)
			s->param.boundary_sign[i][j] = (int)((s->param.boundary_value[i][j] > 0) - (s->param.boundary_value[i][j] < 0));

	/* setting time marching information*/
	s->time.initial	= 0;
	s->time.dt		= timestep;
	s->time.nstep	= (int) ceil(finaltime/timestep);
	s->time.nstep_output	= plot_per_steps;
	s->time.delta	= 0.5;
	s->time.theta	= Newmark_theta;

	/* setting constraint of*/
	set_constraint(s,boundarycondition);

	//
	PetscMalloc1(s->N,&s->v2p);
	set_processor(s, ptr_u);
}


void set_grid(Field_S* s){

	int i;
	int n;
	double dx = s->dx;

	s->Nx = (int) floor((double)s->domain_size[0]/dx)+1;
	s->Ny = (int) floor((double)s->domain_size[1]/dx)+1;
	s->N  = s->Nx*s->Ny;

	s->x_grid = dvector(0,s->Nx-1);
	s->y_grid = dvector(0,s->Ny-1);

	/* setting the x- y- coordinates of the grid*/
	for(i=0; i<s->Nx; i++)
		s->x_grid[i] = s->domain_x[0] + i*dx;

	for(i=0; i<s->Ny; i++)
		s->y_grid[i] = s->domain_y[0] + i*dx;

}

/* compute the sdf of the solid*/
void set_sdf(Field_S* s, int mi, int mf, int ni, int nf, double* x, double* y, double** sdf){

	int		i,j;
	double	center_x,center_y,r,hl_ratio,Lx,Ly;

	// Calculate the siged distance function on each grid point
	// Need to specify the shape parameter
	if(strcmp(s->param.type, "Beam")==0){

		center_x = 0;
		center_y = 0;
		hl_ratio = 1;
		Lx = s->param.length;
		Ly = Lx*hl_ratio;

		for(i=mi;i<mf;i++)
			for(j=ni;j<nf;j++)
				sdf[i][j] = - DMIN(
				DMIN( (x[i]-center_x)+ Lx/2,  Lx/2 - (x[i]-center_x) ) ,
                DMIN( (y[j]-center_y)+ Ly/2,  Ly/2 - (y[j]-center_y) ) );

	}
	else if (strcmp(s->param.type, "Disk")==0){

		center_x = 0.0;
		center_y = 0.0;
		r = 0.3-1e-15;

		for(i=mi;i<mf;i++)
			for(j=ni;j<nf;j++)
				sdf[i][j] = -r + sqrt((x[i]-center_x)*(x[i]-center_x)
				+(y[j] - center_y)*(y[j] - center_y));
	}
	else
		printf("Please enter a new shape of solid");
}

void set_constraint(Field_S* s,char bc){

	int		m = s->Nx,n = s->Ny;
	int		i,j;
	int		toln = m*n;

	// body force
	s->con.body_funct_x = dvector(0,toln-1);
	s->con.body_funct_y = dvector(0,toln-1);

	for(i=0; i<m; i++)
		for(j=0;j<n;j++){
			s->con.body_funct_x[j*m+i]=0;
			s->con.body_funct_y[j*m+i]=0;
		}

	// fsi force
	s->con.fsi_func_x = dvector(0,toln-1);
	s->con.fsi_func_y = dvector(0,toln-1);

	for(i=0; i<m; i++)
		for(j=0;j<n;j++){
			s->con.fsi_func_x[j*m+i]=0;
			s->con.fsi_func_y[j*m+i]=0;
		}

	// Neumann
	s->con.neumann_funct_x = dvector(0,toln-1);
	s->con.neumann_funct_y = dvector(0,toln-1);

	for(i=0; i<m; i++)
		for(j=0;j<n;j++){
			// uni-direction
			if ((s->y_grid[j]) > 0.3)
				s->con.neumann_funct_x[j*m+i]=1.0;
			else if ((s->y_grid[j]) < -0.3)
				s->con.neumann_funct_x[j*m+i]=-1.0;
			else
				s->con.neumann_funct_x[j*m+i]=0;

			// zero force
			s->con.neumann_funct_y[j*m+i]=0;
			//s->con.fsi_func_y[j*m+i]=0;
		}

	//Dirichlet
	s->con.dirichlet_dxi_x = dvector(0,toln-1);
	s->con.dirichlet_dxi_y = dvector(0,toln-1);

	for(i=0; i<m; i++)
		for(j=0;j<n;j++){
			s->con.dirichlet_dxi_x[j*m+i]=0;
			s->con.dirichlet_dxi_y[j*m+i]=0;
		}

	//Surface specification
	s->con.fsineumanndirichlet=(char *)malloc((size_t) ((m*n)*sizeof(char)));

	switch (bc)
	{
		case 'N':
			for(i=0; i<m; i++)
				for(j=0;j<n;j++){
					s->con.fsineumanndirichlet[j*m+i]='N';
			}
			break;

		case 'P':
			for(i=0; i<m; i++)
				for(j=0;j<n;j++){
					if (((s->y_grid[j]+s->dx/2) > 0.3) || ((s->y_grid[j]+s->dx/2) < -0.3))
						s->con.fsineumanndirichlet[j*m+i]='N';
					else
						s->con.fsineumanndirichlet[j*m+i]='D';
			}
			break;

		case 'D':
			for(i=0; i<m; i++)
				for(j=0;j<n;j++){
					s->con.fsineumanndirichlet[j*m+i]='D';
			}
			break;

	}
}

// assign cells to each processor
void set_processor(Field_S* s, AppCtx* ptr_u){

	int i,j,m = s->Nx;
	int m_mid = s->Nx/2, n_mid = s->Ny/2;


	// in the test, rank must be 2
	for(i=0; i<s->Nx; i++)
			for(j=0;j<s->Ny;j++)
			{
				//ptr_u->v2p[j*m+i] = -1;
				//if (i/m_mid == ptr_u->rank)
				s->v2p[j*m+i]=i/m_mid;
			}

	// in the test, rank must be 4
	/*for(i=0; i<s->Nx; i++)
			for(j=0;j<s->Ny;j++)
			{
				ptr_u->v2p[j*m+i]=(j/m_mid)*2 + i/m_mid;
			}
	*/

}
