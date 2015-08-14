#ifndef INIT_SETTING_H
#define INIT_SETTING_H 1


// set index, bound function, inital
#include "Field_s.h"
#include <petsc.h>
#include <petscao.h>

struct invol
{
	int xix[6];
	int xiy[6];
	int stagx_cell [2];
	int stagy_cell [2];
	int i,j;
};


void set_index(Index_S* ind, Grid_S* g, AppCtx* ptr_u, int **bs, char *fsineumanndirichlet);
struct invol involvedIndices_grid(int cell, int m);

#endif  /* INIT_SETTING_H */
