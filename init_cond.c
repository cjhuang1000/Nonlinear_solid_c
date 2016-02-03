#include "init_cond.h"
#include <string.h>
#include <petscvec.h>

void set_initial(Index_S* ptr_i, Solid* ptr_s, Field_S* ptr_fs)
{

    int       i,j;
    PetscInt  n;

    // initialize field_s
    VecCreate(PETSC_COMM_WORLD,&ptr_fs->xi);
    VecSetSizes(ptr_fs->xi,ptr_i->xix_N + ptr_i->xiy_N,ptr_i->xi_gloN);
    VecSetFromOptions(ptr_fs->xi);

    VecDuplicate(ptr_fs->xi,&ptr_fs->dxi);
    VecDuplicate(ptr_fs->xi,&ptr_fs->ddxi);
    VecDuplicate(ptr_fs->xi,&ptr_fs->inc_dxi);

    /* set initial condition*/
    if (strcmp(ptr_s->initial,"rest")==0)
    {
        VecZeroEntries(ptr_fs->xi);
        VecZeroEntries(ptr_fs->dxi);
        VecZeroEntries(ptr_fs->ddxi);
    }

    if (strcmp(ptr_s->initial,"specified")==0)
    {
        VecSet(ptr_fs->xi,0.06);
        VecSet(ptr_fs->dxi,0.0);
        VecSet(ptr_fs->ddxi,0.0);
    }
    VecZeroEntries(ptr_fs->inc_dxi);

    VecAssemblyBegin(ptr_fs->xi);
    VecAssemblyEnd(ptr_fs->xi);
    VecAssemblyBegin(ptr_fs->dxi);
    VecAssemblyEnd(ptr_fs->dxi);
    VecAssemblyBegin(ptr_fs->ddxi);
    VecAssemblyEnd(ptr_fs->ddxi);
    VecAssemblyBegin(ptr_fs->inc_dxi);
    VecAssemblyEnd(ptr_fs->inc_dxi);

    VecCreate(PETSC_COMM_WORLD,&ptr_fs->traction);
    VecSetSizes(ptr_fs->traction,ptr_i->xix_Ncell_Neumann + ptr_i->xiy_Ncell_Neumann,ptr_i->xi_gloNcell_Neumann);
    VecSetFromOptions(ptr_fs->traction);

    VecSet(ptr_fs->traction,0.0);

}
