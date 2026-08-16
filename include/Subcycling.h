#pragma once

#ifndef _PDS_SUBCYCLING_H
#define _PDS_SUBCYCLING_H

#include "Typedefs.h"

/*
 * These helpers are the deliberately small interface between the
 * drift-mode subcycling driver and the existing force calculation.  When
 * subcycling is inactive they preserve the normal LocalSegForces behavior.
 */
int  SubcyclingIsActive(Home_t *home);
int  SubcyclingUseBaseForces(Home_t *home);
int  SubcyclingSelectSegmentPair(Home_t *home, Node_t *node1, Node_t *node2,
                                 Node_t *node3, Node_t *node4);
void SubcyclingAddCachedForces(Home_t *home);

#endif  /* _PDS_SUBCYCLING_H */
