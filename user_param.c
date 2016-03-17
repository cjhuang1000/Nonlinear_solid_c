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

	int		solid_resolution	= 200;
	double	domain[4]			={-1.0,1.0,-1.0,1.0};
	double  char_length			=1.0-1e-10;

	double  density				=1.0;
	double	nu					=0.2;
	double	E					=12e7;
	char	shape[]				="Beam"; //"Beam","Disk","Elliptic"

	char	boundarycondition	='P'; // N: all Neumann, F:all FSI, D: all Dirichlet
									  // M: all mixed, P: Neumann+Dirichlet
    double  solid_threshold     =1e-5;
    double  damping[2]			= {0.0, 0.0};

	double	timestep			=15e-5;
	double	finaltime			=15e-5*2;
	double	Newmark_theta		=1.0/4;
	int		plot_per_steps		=1000;

	char	initial[]			="rest"; // rest, specified

	int		i,j;

	printf("[%s] \n", __func__);
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
	s->param.rho			= density;
	s->param.lambda			= E*nu/(1.0+nu)/(1.0-2.0*nu);
	s->param.mu				= E/2.0/(1.0+nu);
	s->param.threshold		= solid_threshold;
	strcpy( s->param.type, shape);
	strcpy( s->param.initial, initial);
	s->param.boundcond 		= boundarycondition;
	s->param.damping[0] 		= damping[0];
	s->param.damping[1] 		= damping[1];
	s->param.boundary_value 	= dmatrix(0,(s->Nx-1)*2,0,(s->Ny-1)*2);
	s->param.boundary_sign  	= imatrix(0,(s->Nx-1)*2,0,(s->Ny-1)*2);

	set_sdf(s,0,s->Nx,0, s->Ny, s->x_grid, s->y_grid, s->param.boundary_value,'s');				// set sdf

	for(i=0; i<=(s->Nx-1)*2; i++)
		for(j=0; j<=(s->Ny-1)*2; j++)
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
void set_sdf(Field_S* s, int mx, int nx, int my, int ny, double* x, double* y, double** sdf, char type){

	int		i,j, rank, sdf_nx, sdf_ny;
	double	center_x,center_y,r,hl_ratio,Lx,Ly, temp1, temp2;
	double  *refined_grid_x, *refined_grid_y;

	sdf_nx = (nx-1)*2+1;
	sdf_ny = (ny-1)*2+1;

	refined_grid_x = dvector(0,sdf_nx-1);
	refined_grid_y = dvector(0,sdf_ny-1);

	for(i=0;i<nx;i++){
		refined_grid_x[i*2] 	= x[i];
		refined_grid_x[i*2+1] 	= (x[i]+x[i+1])/2.0;
	}
	for(i=0;i<ny;i++){
		refined_grid_y[i*2] 	= y[i];
		refined_grid_y[i*2+1] 	= (y[i]+y[i+1])/2.0;
	}
	refined_grid_x[sdf_nx-1]=x[nx-1];
	refined_grid_y[sdf_ny-1]=y[ny-1];

	// Calculate the siged distance function on each grid point
	// Need to specify the shape parameter
	if(strcmp(s->param.type, "Beam")==0){

		center_x = 0.0;
		center_y = 0.0;
		hl_ratio = 0.1;
		Lx = s->param.length;
		Ly = Lx*hl_ratio;
		MPI_Comm_rank(MPI_COMM_WORLD,&rank);

		switch(type)
		{
		case 's':
			for(i=0;i<sdf_nx;i++){
				for(j=0;j<sdf_ny;j++){
					temp1 = DMIN( (refined_grid_x[i]-center_x)+ Lx/2,  Lx/2 - (refined_grid_x[i]-center_x) );
					temp2 = DMIN( (refined_grid_y[j]-center_y)+ Ly/2,  Ly/2 - (refined_grid_y[j]-center_y) ) ;
					sdf[i][j] = -DMIN(temp1, temp2);
				}
			}
			break;
		case 'f':
			for(i=0;i<nx;i++){
				for(j=0;j<ny;j++){
					temp1 = DMIN( (x[i]-center_x)+ Lx/2,  Lx/2 - (x[i]-center_x) );
					temp2 = DMIN( (y[j]-center_y)+ Ly/2,  Ly/2 - (y[j]-center_y) ) ;
					sdf[i][j] = -DMIN(temp1, temp2);
				}
			}
		break;
		}


	}
	else if (strcmp(s->param.type, "Disk")==0){

		center_x = 0.0;
		center_y = 0.0;
		r = 0.3-1e-15;

		switch(type)
		{
		case 's':
			for(i=0;i<sdf_nx;i++)
				for(j=0;j<sdf_ny;j++)
					sdf[i][j] = -r + sqrt((refined_grid_x[i]-center_x)*(refined_grid_x[i]-center_x)
					+(refined_grid_y[j] - center_y)*(refined_grid_y[j] - center_y));
			break;
		case 'f':
			for(i=0;i<nx;i++)
				for(j=0;j<ny;j++)
					sdf[i][j] = -r + sqrt((x[i]-center_x)*(x[i]-center_x)
					+(y[j] - center_y)*(y[j] - center_y));
			break;
		}
	}
	else
		printf("Please enter a new shape of solid");

	free_dvector(refined_grid_x,0,sdf_nx-1);
	free_dvector(refined_grid_y,0,sdf_ny-1);
}

void set_constraint(Field_S* s,char bc){

	int		m = s->Nx,n = s->Ny;
	int		i,j;
	int		toln = m*n;
	double  strength, ss, x, y;

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

	strength = 2.85e4;
	ss = 1.0/0.1;

	for(i=0; i<m; i++)
		for(j=0;j<n;j++){
			// uni-direction
			s->con.neumann_funct_x[j*m+i]=0.0;

			x = s->x_grid[i];
			y = s->y_grid[j]+ s->dx/2;

			// zero force
			if ((x > -0.45) && (x < 0.45))
				s->con.neumann_funct_y[j*m+i]=-strength/2.0;
			else
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
		//'N'*(x>-0.49)  + 'D'*(x<=-0.49);
		case 'P':
			for(i=0; i<m; i++)
				for(j=0;j<n;j++){
					if ((s->x_grid[i]+s->dx/2) > -0.45)
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


	// uni processor
	//for(i=0; i<s->Nx; i++)
	//		for(j=0;j<s->Ny;j++)
	//		{
	//			s->v2p[j*m+i]=0;
	//		}

	// 2 processes
	// +-------+-------+
	// |   0   |   1   |
	// +-------+-------+

	// in the test, rank must be 2
	for(i=0; i<s->Nx; i++)
			for(j=0;j<s->Ny;j++)
			{
				s->v2p[j*m+i]=i/m_mid;
			}

	// 4 processes
	// +-------+-------+
	// |   2   |   3   |
	// +-------+-------+
	// |   0   |   1   |
	// +-------+-------+
	// in the test, rank must be 4
	/*for(i=0; i<s->Nx; i++)
			for(j=0;j<s->Ny;j++)
			{
				ptr_u->v2p[j*m+i]=(j/m_mid)*2 + i/m_mid;
			}
	*/

}
