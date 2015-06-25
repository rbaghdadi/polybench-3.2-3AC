/**
 * durbin.c: This file is part of the PolyBench/C 3.2 test suite.
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
/* Default data type is double, default size is 4000. */
#include "durbin.h"


/* Array initialization. */
static
void init_array (int n,
		 DATA_TYPE POLYBENCH_2D(y,N,N,n,n),
		 DATA_TYPE POLYBENCH_2D(sum,N,N,n,n),
		 DATA_TYPE POLYBENCH_1D(alpha,N,n),
		 DATA_TYPE POLYBENCH_1D(beta,N,n),
		 DATA_TYPE POLYBENCH_1D(r,N,n))
{
  int i, j;

  for (i = 0; i < n; i++)
    {
      alpha[i] = i;
      beta[i] = (i+1)/n/2.0;
      r[i] = (i+1)/n/4.0;
      for (j = 0; j < n; j++) {
	y[i][j] = ((DATA_TYPE) i*j) / n;
	sum[i][j] = ((DATA_TYPE) i*j) / n;
      }
    }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int n,
		 DATA_TYPE POLYBENCH_1D(out,N,n))

{
  int i;

  for (i = 0; i < n; i++) {
    fprintf (stderr, DATA_PRINTF_MODIFIER, out[i]);
    if (i % 20 == 0) fprintf (stderr, "\n");
  }
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_durbin(int n,
		   DATA_TYPE POLYBENCH_2D(y,N,N,n,n),
		   DATA_TYPE POLYBENCH_2D(sum,N,N,n,n),
		   DATA_TYPE POLYBENCH_1D(alpha,N,n),
		   DATA_TYPE POLYBENCH_1D(beta,N,n),
		   DATA_TYPE POLYBENCH_1D(r,N,n),
		   DATA_TYPE POLYBENCH_1D(out,N,n))
{
  int i, k;

  DATA_TYPE temp1, temp2, temp3, temp4, temp5,
                   temp12, temp14, temp15, temp16;
  int temp6, temp7, temp8, temp11, temp10, temp9, temp13;

#pragma scop
  y[0][0] = r[0];
  beta[0] = 1;
  alpha[0] = r[0];
  for (k = 1; k < _PB_N; k++)
    {
      temp6 = k-1;
      temp2 = alpha[temp6] * beta[temp6];
      temp1 = alpha[temp6] * temp2;
      temp14 = beta[temp6] - temp1;
      beta[k] = temp14;
      sum[0][k] = r[k];
      temp12 =  k - 1;
      for (i = 0; i <=temp12; i++)
	{
	  temp7 = k-i;
	  temp8 = temp7-1;
	  temp3 = r[temp8] * y[i][temp6];
	  temp11 = i+1;
	  temp15 = sum[i][k] + temp3;
	  sum[temp11][k] = temp15; 
	}
      temp5 = -sum[k][k];
      alpha[k] = temp5 * beta[k];
      temp13 = k-1;
      for (i = 0; i <= temp13; i++)
	{
	  temp9 = k-i;
	  temp10 = temp9-1;
	  temp4 = alpha[k] * y[temp10][temp6];
	  temp16 = y[i][temp6] + temp4;
	  y[i][k] = temp16; 
	}
      y[k][k] = alpha[k];
    }
  for (i = 0; i < _PB_N; i++)
    out[i] = y[i][_PB_N-1];
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int n = N;

  /* Variable declaration/allocation. */
  POLYBENCH_2D_ARRAY_DECL(y, DATA_TYPE, N, N, n, n);
  POLYBENCH_2D_ARRAY_DECL(sum, DATA_TYPE, N, N, n, n);
  POLYBENCH_1D_ARRAY_DECL(alpha, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(beta, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(r, DATA_TYPE, N, n);
  POLYBENCH_1D_ARRAY_DECL(out, DATA_TYPE, N, n);


  /* Initialize array(s). */
  init_array (n,
	      POLYBENCH_ARRAY(y),
	      POLYBENCH_ARRAY(sum),
	      POLYBENCH_ARRAY(alpha),
	      POLYBENCH_ARRAY(beta),
	      POLYBENCH_ARRAY(r));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_durbin (n,
		 POLYBENCH_ARRAY(y),
		 POLYBENCH_ARRAY(sum),
		 POLYBENCH_ARRAY(alpha),
		 POLYBENCH_ARRAY(beta),
		 POLYBENCH_ARRAY(r),
		 POLYBENCH_ARRAY(out));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(n, POLYBENCH_ARRAY(out)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(y);
  POLYBENCH_FREE_ARRAY(sum);
  POLYBENCH_FREE_ARRAY(alpha);
  POLYBENCH_FREE_ARRAY(beta);
  POLYBENCH_FREE_ARRAY(r);
  POLYBENCH_FREE_ARRAY(out);

  return 0;
}
