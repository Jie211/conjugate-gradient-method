#include <stdio.h>
#include <stdlib.h>
#include <omp.h>
#include "function.h"
#include "cgm.h"
int main(int argc, char const* argv[])
{
  int N;
  int M;
  int MAX;
  double *b;
  double *x;

  double *val;
  int *col;
  int *ptr;

  if(argc!=6)
  {
    printf("Usage \n");
    printf("%s colval ptr bx maxloop thread\n", argv[0]);
    exit(0);
  }
  printf("--OpenMP--\n");
#ifdef EBUG
  printf("OpenMP Max Threads %d\n", omp_get_max_threads());
  printf("OpenMP Threads set to %d\n", atoi(argv[5]));
#endif
  omp_set_num_threads(atoi(argv[5]));
  printf("--OpenMP over--\n");
  printf("--check error--\n");
#ifdef EBUG
  printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
  printf("get headpart and check headerror\n");
#endif
  gethead(argv[1],argv[2],argv[3],&N,&M); 
#ifdef EBUG
  printf("get N = %d M = %d\n", N, M);
#endif
  b=(double *)malloc(sizeof(double)*N);
  x=(double *)malloc(sizeof(double)*N);
  if(!b || !x)
  {
    printf("malloc b x error");
    if(b)
      free(b);
    if(x)
      free(x);
    exit(1);
  }


  val=(double *)malloc(sizeof(double)*M);
  col=(int *)malloc(sizeof(int)*M);
  ptr=(int *)malloc(sizeof(int)*(N+1));
  if(!val || !col || !ptr)
  {
    printf("malloc val col ptr error");
    if(b)
      free(b);
    if(x)
      free(x);
    if(val)
      free(val);
    if(col)
      free(col);
    if(ptr)
      free(ptr);
    exit(1);
  }
#ifdef EBUG  
  printf("read data\n");
#endif
  getdata(argv[1],argv[2],argv[3],col,ptr,val,b,x,N,M);
#ifdef EBUG
  printf("read data over\n");
  printf("@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@@\n");
#endif
  printf("--check error over--\n");
  int alpha=atoi(argv[4]);
  MAX = alpha;
  printf("--result--\n");
  cgm_CRS(val, col, ptr, b, x, N, 1e-10, MAX, M);

  free(b);
  free(x);
  free(val);
  free(col);
  free(ptr);

  return 0;
}
