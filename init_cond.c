#include "init_cond.h"
#include <string.h>
#include <petscvec.h>

void set_initial(Field_S* s)
{

    int       i,j;
    PetscInt  n;

    // initialize field_s
    VecCreate(PETSC_COMM_WORLD,&(s->xi));
    VecSetSizes(s->xi,s->ind.xix_N + s->ind.xiy_N,s->ind.xi_gloN);
    VecSetFromOptions(s->xi);

    VecDuplicate(s->xi,&s->dxi);
    VecDuplicate(s->xi,&s->ddxi);
    VecDuplicate(s->xi,&s->inc_dxi);

    /* set initial condition*/
    if (strcmp(s->param.initial,"rest")==0)
    {
        VecZeroEntries(s->xi);
        VecZeroEntries(s->dxi);
        VecZeroEntries(s->ddxi);
    }

    if (strcmp(s->param.initial,"specified")==0)
    {
        VecSet(s->xi,0.06);
        VecSet(s->dxi,0.0);
        VecSet(s->ddxi,0.0);
    }
    VecZeroEntries(s->inc_dxi);

    VecAssemblyBegin(s->xi);
    VecAssemblyEnd(s->xi);
    VecAssemblyBegin(s->dxi);
    VecAssemblyEnd(s->dxi);
    VecAssemblyBegin(s->ddxi);
    VecAssemblyEnd(s->ddxi);
    VecAssemblyBegin(s->inc_dxi);
    VecAssemblyEnd(s->inc_dxi);

    VecCreate(PETSC_COMM_WORLD,&s->traction);
    VecSetSizes(s->traction,s->ind.xix_Ncell_Neumann + s->ind.xiy_Ncell_Neumann,s->ind.xi_gloNcell_Neumann);
    VecSetFromOptions(s->traction);

    VecSet(s->traction,0.0);

}
