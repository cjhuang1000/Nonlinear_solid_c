#ifndef MATRICES_UPDATE_H
#define MATRICES_UPDATE_H 1

#include "Field_s.h"
#include <petscmat.h>

void compute_matricesNonlinearStructure_update(Matrices_S* ptr_ms, Index_S* ptr_i, Grid_S* ptr_g,
		Solid* ptr_s, AppCtx* ptr_u, Field_S* ptr_f);

#endif  /* MATRICES_UPDATE_H */
