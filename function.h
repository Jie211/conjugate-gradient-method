#ifndef FUNCTION_H_INCLUDED__
#define FUNCTION_H_INCLUDED__

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>
#include <sys/time.h>

extern void
bxcreate(const double *a,
         double *b,
         double *x,
         const int n_data);

extern void 
CSR(const double *a,
    double *val,
    int *col,
    int *ptr,
    const int n_data,
    const int m);

extern void
display_mat(const double *a, 
            const int n_data);

extern void
display_vec(const double *a, 
            const int n_data);

extern void
list_vec_i(const int *a, 
            const int n_data);

extern void
list_vec_d(const double *a, 
            const int n_data);

extern void
CSR(const double *a,
    double *val,
    int *col,
    int *ptr,
    const int n_data,
    const int m);

extern void 
usage(const char *a,const char *b,const char *x, int *ndata, int *m);

extern void 
getmatrix(const char *file,double *a,int ndata,int m);

extern void
getvector(const char *file,double *v,int ndata);

extern void
gethead(const char *file1, const char *file2, const char *file3, int *N, int *NZ);

extern void
getdata(const char *file1, const char *file2, const char *file3,int *col,int *ptr,double *val,double *b,double *x, int N, int NZ);

extern double
gettimeofday_sec();
#endif //FUNCTION_H_INCLUDED__

