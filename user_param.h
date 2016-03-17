#ifndef USER_PARAM_H
#define USER_PARAM_H 1

//DONE

#include "Field_s.h"
#include <petsc.h>
#include <petscsys.h>

void user_param (Field_S* s, AppCtx* ptr_u);
void set_sdf(Field_S* s, int mx, int nx, int my, int ny, double* x, double* y, double** sdf, char type);

#endif  /* USER_PARAM_H */
