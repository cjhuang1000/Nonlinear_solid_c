#ifndef INIT_SETTING_H
#define INIT_SETTING_H 1
// DONE

// set index, bound function
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


void set_index(Field_S* s, AppCtx* ptr_u);
struct invol involvedIndices_grid(int cell, int m);

#endif  /* INIT_SETTING_H */
