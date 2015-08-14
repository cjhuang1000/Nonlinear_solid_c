#ifndef USER_PARAM_H
#define USER_PARAM_H 1


#include "Field_s.h"
#include <petsc.h>
#include <petscsys.h>

void user_param (Grid_S* ptr_g, Solid* ptr_s, TimeMarching* ptr_t,Constraint_S* ptr_c, AppCtx* ptr_u);

#endif  /* USER_PARAM_H */
