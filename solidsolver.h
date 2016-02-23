/*
 * solidsolver.h
 *
 *  Created on: Feb 17, 2016
 *      Author: peggyhuang
 */

#ifndef SOLIDSOLVER_H_
#define SOLIDSOLVER_H_

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#include <petscsys.h>
#include <petsc.h>
#include <petscvec.h>
#include <petscmat.h>
#include <petscksp.h>
#include <petscis.h>


#include "Field_s.h"
#include "user_param.h"
#include "init_setting.h"
#include "boun_func.h"
#include "init_cond.h"

void SolidInitialize(Field_S* s, AppCtx* ptr_u);
void SolidSolver(Field_S* s);
void SolidOutput(Field_S* s);
void SolidFinalize(Field_S* s);

#endif /* SOLIDSOLVER_H_ */
