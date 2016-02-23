/*
 * interface.h
 *
 *  Created on: Feb 18, 2016
 *      Author: peggyhuang
 */

#ifndef INTERFACE_H_
#define INTERFACE_H_

#include "interface_marker.h"
#include "Field_s.h"
#include "fluid.h"

void InterfaceInitialize(Lag_marker* mk, Field_S* s, Field_F* f, AppCtx* ptr_u);
void InterfaceUpdate	(Lag_marker* mk, Field_S* s, Field_F* f);
void InterfaceTraction	(Lag_marker* mk, Field_S* s, Field_F* f);
void InterfaceOutput	(Lag_marker* mk);
void InterfaceFinalize	(Lag_marker* mk);

#endif /* INTERFACE_H_ */
