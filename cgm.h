#ifndef CGM_H_INCLUDED__
#define CGM_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>

extern void
cgm_CRS(double *val,
    const int *col,
    const int *ptr,
    double *bvec,
    double *xvec,
    const int ndata,
    const double eps,
    const int i_max,
    const int m);


#endif //CGM_H_INCLUDED__

