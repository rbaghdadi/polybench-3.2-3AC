/**
 * adi.c: This file is part of the PolyBench/C 3.2 test suite.
 *
 *
 * Contact: Louis-Noel Pouchet <pouchet@cse.ohio-state.edu>
 * Web address: http://polybench.sourceforge.net
 */
#include <stdio.h>
#include <unistd.h>
#include <string.h>
#include <math.h>

/* Include polybench common header. */
#include <polybench.h>

/* Include benchmark-specific header. */
/* Default data type is double, default size is 10x1024x1024. */
#include "adi.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_2D(X,N,N,n,n),
		 DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
		 DATA_TYPE POLYBENCH_2D(B,N,N,n,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++)
      {
	X[i][j] = ((DATA_TYPE) i*(j+1) + 1) / n;
	A[i][j] = ((DATA_TYPE) i*(j+2) + 2) / n;
	B[i][j] = ((DATA_TYPE) i*(j+3) + 3) / n;
      }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_2D(X,N,N,n,n))

{
  int i, j;

  for (i = 0; i < n; i++)
    for (j = 0; j < n; j++) {
      fprintf(stderr, DATA_PRINTF_MODIFIER, X[i][j]);
      if ((i * N + j) % 20 == 0) fprintf(stderr, "\n");
    }
  fprintf(stderr, "\n");
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_adi(int tsteps,
		int n,
		DATA_TYPE POLYBENCH_2D(X,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(A,N,N,n,n),
		DATA_TYPE POLYBENCH_2D(B,N,N,n,n))
{
  int t, i1, i2;

DATA_TYPE temp1, temp2, temp3, temp4, temp5, temp6, temp7, temp8, temp9,
	  temp10, temp11, temp12, temp26, temp27, temp28, temp29, temp30,
	  temp31;

int temp13, temp14, temp15, temp16, temp17,
	  temp18, temp19, temp20, temp21, temp22, temp24, temp23, temp25;

#pragma scop
  for (t = 0; t < _PB_TSTEPS; t++)
    {
      for (i1 = 0; i1 < _PB_N; i1++)
	for (i2 = 1; i2 < _PB_N; i2++)
	  {
	    temp4 = X[i1][i2-1] * A[i1][i2];
	    temp1 = temp4 / B[i1][i2-1];
	    temp26 = X[i1][i2] - temp1;
	    X[i1][i2] = temp26; 
	    temp5 = A[i1][i2] * A[i1][i2];
	    temp13 = i2-1;
	    temp2 = temp5 / B[i1][temp13];
	    temp27 = B[i1][i2] - temp2;
	    B[i1][i2] = temp27; 
	  }

      for (i1 = 0; i1 < _PB_N; i1++)
      {
	temp28 = X[i1][_PB_N-1] / B[i1][_PB_N-1]; 
	X[i1][_PB_N-1] = temp28; 
      }

      for (i1 = 0; i1 < _PB_N; i1++)
	for (i2 = 0; i2 < _PB_N-2; i2++)
	  {
	    temp14 = 1-i2;
	    temp15 = _PB_N-i2-3;
	    temp6 = X[i1][_PB_N-temp14] * A[i1][temp15];
	    temp16 = _PB_N-2-i2;
	    temp3 = (X[i1][temp16] - temp6);
	    temp17 = _PB_N-i2-2;
	    temp18 = _PB_N-3-i2;
	    temp29 = temp3 / B[i1][temp18];
	    X[i1][temp17] = temp29; 
	  }

      for (i1 = 1; i1 < _PB_N; i1++)
	for (i2 = 0; i2 < _PB_N; i2++) {
	  temp19 = i1-1;
	  temp7 = X[temp19][i2] * A[i1][i2]; 
	  temp20 = i1-1;
	  temp8 = temp7 / B[temp20][i2];
	  temp31 =  X[i1][i2] - temp8;
	  X[i1][i2] = temp31;
	  temp9 = A[i1][i2] * A[i1][i2];
	  temp21 = i1-1;
	  temp10 = temp9 / B[temp21][i2];
	  temp30 = B[i1][i2] - temp10;
	  B[i1][i2] = temp30;
	}

      for (i2 = 0; i2 < _PB_N; i2++)
      {
	temp22 = X[_PB_N-1][i2] / B[_PB_N-1][i2];
	X[_PB_N-1][i2] = temp22; 
      }

      for (i1 = 0; i1 < _PB_N-2; i1++)
	for (i2 = 0; i2 < _PB_N; i2++)
	{
	  temp22 = _PB_N-i1-3;
	  temp23 = _PB_N-3-i1;
	  temp11 = X[temp22][i2] * A[temp23][i2];
	  temp24 = _PB_N-2-i1;
	  temp12 = (X[temp24][i2] - temp11);
	  temp25 = _PB_N-2-i1;
	  temp23 = temp12 / B[temp25][i2];
	  X[temp25][i2] = temp23; 
	}
    }
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;
  int tsteps = TSTEPS;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(X, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(A, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(B, DATA_TYPE, N, N, n, n);


  /* Initialize array(s). */
  init_array (n, POLYBENCH_ARRAY(X), POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_adi (tsteps, n, POLYBENCH_ARRAY(X),
	      POLYBENCH_ARRAY(A), POLYBENCH_ARRAY(B));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(X)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(X);
  POLYBENCH_FREE_ARRAY(A);
  POLYBENCH_FREE_ARRAY(B);

  return 0;
}
