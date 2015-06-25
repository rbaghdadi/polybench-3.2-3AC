/**
 * fdtd-apml.c: This file is part of the PolyBench/C 3.2 test suite.
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
/* Default data type is double, default size is 256x256x256. */
#include "fdtd-apml.h"


/* Array initialization. */
static
void init_array (int cz,
		 int cxm,
		 int cym,
		 DATA_TYPE *mui,
		 DATA_TYPE *ch,
		 DATA_TYPE POLYBENCH_2D(Ax,CZ+1,CYM+1,cz+1,cym+1),
		 DATA_TYPE POLYBENCH_2D(Ry,CZ+1,CYM+1,cz+1,cym+1),
		 DATA_TYPE POLYBENCH_3D(Ex,CZ+1,CYM+1,CXM+1,cz+1,cym+1,cxm+1),
		 DATA_TYPE POLYBENCH_3D(Ey,CZ+1,CYM+1,CXM+1,cz+1,cym+1,cxm+1),
		 DATA_TYPE POLYBENCH_3D(Hz,CZ+1,CYM+1,CXM+1,cz+1,cym+1,cxm+1),
		 DATA_TYPE POLYBENCH_1D(czm,CZ+1,cz+1),
		 DATA_TYPE POLYBENCH_1D(czp,CZ+1,cz+1),
		 DATA_TYPE POLYBENCH_1D(cxmh,CXM+1,cxm+1),
		 DATA_TYPE POLYBENCH_1D(cxph,CXM+1,cxm+1),
		 DATA_TYPE POLYBENCH_1D(cymh,CYM+1,cym+1),
		 DATA_TYPE POLYBENCH_1D(cyph,CYM+1,cym+1))
{
  int i, j, k;
  *mui = 2341;
  *ch = 42;
  for (i = 0; i <= cz; i++)
    {
      czm[i] = ((DATA_TYPE) i + 1) / cxm;
      czp[i] = ((DATA_TYPE) i + 2) / cxm;
    }
  for (i = 0; i <= cxm; i++)
    {
      cxmh[i] = ((DATA_TYPE) i + 3) / cxm;
      cxph[i] = ((DATA_TYPE) i + 4) / cxm;
    }
  for (i = 0; i <= cym; i++)
    {
      cymh[i] = ((DATA_TYPE) i + 5) / cxm;
      cyph[i] = ((DATA_TYPE) i + 6) / cxm;
    }

  for (i = 0; i <= cz; i++)
    for (j = 0; j <= cym; j++)
      {
	Ry[i][j] = ((DATA_TYPE) i*(j+1) + 10) / cym;
	Ax[i][j] = ((DATA_TYPE) i*(j+2) + 11) / cym;
	for (k = 0; k <= cxm; k++)
	  {
	    Ex[i][j][k] = ((DATA_TYPE) i*(j+3) + k + 1) / cxm;
	    Ey[i][j][k] = ((DATA_TYPE) i*(j+4) + k + 2) / cym;
	    Hz[i][j][k] = ((DATA_TYPE) i*(j+5) + k + 3) / cz;
	  }
      }
}


/* DCE code. Must scan the entire live-out data.
   Can be used also to check the correctness of the output. */
static
void print_array(int cz,
		 int cxm,
		 int cym,
		 DATA_TYPE POLYBENCH_3D(Bza,CZ+1,CYM+1,CXM+1,cz+1,cym+1,cxm+1),
		 DATA_TYPE POLYBENCH_3D(Ex,CZ+1,CYM+1,CXM+1,cz+1,cym+1,cxm+1),
		 DATA_TYPE POLYBENCH_3D(Ey,CZ+1,CYM+1,CXM+1,cz+1,cym+1,cxm+1),
		 DATA_TYPE POLYBENCH_3D(Hz,CZ+1,CYM+1,CXM+1,cz+1,cym+1,cxm+1))
{
  int i, j, k;

  for (i = 0; i <= cz; i++)
    for (j = 0; j <= cym; j++)
      for (k = 0; k <= cxm; k++) {
	fprintf(stderr, DATA_PRINTF_MODIFIER, Bza[i][j][k]);
	fprintf(stderr, DATA_PRINTF_MODIFIER, Ex[i][j][k]);
	fprintf(stderr, DATA_PRINTF_MODIFIER, Ey[i][j][k]);
	fprintf(stderr, DATA_PRINTF_MODIFIER, Hz[i][j][k]);
	if ((i * cxm + j) % 20 == 0) fprintf(stderr, "\n");
      }
  fprintf(stderr, "\n");
}


/* Main computational kernel. The whole function will be timed,
   including the call and return. */
static
void kernel_fdtd_apml(int cz,
		      int cxm,
		      int cym,
		      DATA_TYPE mui,
		      DATA_TYPE ch,
		      DATA_TYPE POLYBENCH_2D(Ax,CZ+1,CYM+1,cz+1,cym+1),
		      DATA_TYPE POLYBENCH_2D(Ry,CZ+1,CYM+1,cz+1,cym+1),
		      DATA_TYPE POLYBENCH_2D(clf,CYM+1,CXM+1,cym+1,cxm+1),
		      DATA_TYPE POLYBENCH_2D(tmp,CYM+1,CXM+1,cym+1,cxm+1),
		      DATA_TYPE POLYBENCH_3D(Bza,CZ+1,CYM+1,CXM+1,cz+1,cym+1,cxm+1),
		      DATA_TYPE POLYBENCH_3D(Ex,CZ+1,CYM+1,CXM+1,cz+1,cym+1,cxm+1),
		      DATA_TYPE POLYBENCH_3D(Ey,CZ+1,CYM+1,CXM+1,cz+1,cym+1,cxm+1),
		      DATA_TYPE POLYBENCH_3D(Hz,CZ+1,CYM+1,CXM+1,cz+1,cym+1,cxm+1),
		      DATA_TYPE POLYBENCH_1D(czm,CZ+1,cz+1),
		      DATA_TYPE POLYBENCH_1D(czp,CZ+1,cz+1),
		      DATA_TYPE POLYBENCH_1D(cxmh,CXM+1,cxm+1),
		      DATA_TYPE POLYBENCH_1D(cxph,CXM+1,cxm+1),
		      DATA_TYPE POLYBENCH_1D(cymh,CYM+1,cym+1),
		      DATA_TYPE POLYBENCH_1D(cyph,CYM+1,cym+1))
{
  int iz, iy, ix;

  DATA_TYPE temp1, temp2, temp3, temp4, temp5, temp6, tmp1,
  temp7, temp8, temp9, temp10, temp11, temp12, temp13, temp14, temp15, temp16, temp17, temp18, temp19, temp20, temp21, temp22, temp23, temp24,    
  temp25, temp26, temp27, temp28, temp29, temp30, temp31, temp32, temp33, temp34, temp35, temp36, temp37, temp38, temp39, temp40, temp41, temp42, 
  temp43, temp44, temp45, temp46, temp47, temp48, temp49, temp50, temp51, temp52, temp53, temp54, temp55, temp56, temp57, temp58; 
  

  int temp60, temp59, temp61, temp62;

#pragma scop
  for (iz = 0; iz < _PB_CZ; iz++)
    {
      for (iy = 0; iy < _PB_CYM; iy++)
	{
	  for (ix = 0; ix < _PB_CXM; ix++)
	    {
	      temp59 = iy+1;
	      temp1 = Ex[iz][iy][ix] - Ex[iz][temp59][ix]; 
	      temp60 = ix+1;
	      temp2 = Ey[iz][iy][temp60] - Ey[iz][iy][ix];
	      clf[iz][iy] = temp1 + temp2;
	      temp3 = (cymh[iy] / cyph[iy]);
	      temp4 = (ch / cyph[iy]);
	      temp5 = temp4 * clf[iz][iy];
	      temp7 =  Bza[iz][iy][ix] - temp5;
	      tmp1 = temp3 * temp7;
	      temp6 = (cxmh[ix] /cxph[ix]);
	      temp8 = mui * czp[iz] / cxph[ix];
	      temp9 = mui * czm[iz] ;
	      temp10 = temp9/ cxph[ix];
	      temp11 = (temp10) * Bza[iz][iy][ix];
	      temp12 = temp6 * Hz[iz][iy][ix];
	      temp13 = (temp8) * tmp1;
	      temp14 = temp13 - temp11;
	      Hz[iz][iy][ix] = temp12 + temp14;
	      Bza[iz][iy][ix] = tmp1;
	    }
	    temp61 = iy+1;
	    temp15 = Ex[iz][iy][_PB_CXM] - Ex[iz][temp61][_PB_CXM];
	    temp16 = Ry[iz][iy] - Ey[iz][iy][_PB_CXM];
	  clf[iz][iy] = temp15 + temp16;
	  temp17 = cymh[iy] / cyph[iy];
	  temp18 = (temp17) * Bza[iz][iy][_PB_CXM];
	  temp19 = ch / cyph[iy];
	  temp20 = (temp19) * clf[iz][iy];
	  tmp1 = temp18 - temp20;
	  temp21 = cxmh[_PB_CXM] / cxph[_PB_CXM];
	  temp22 = mui * czp[iz];
	  temp23 = temp22 / cxph[_PB_CXM];
	  temp24 = (temp23) * tmp1;
	  temp25 = mui * czm[iz];
	  temp26 = temp25 / cxph[_PB_CXM];
	  temp27 = (temp21) * Hz[iz][iy][_PB_CXM];
	  temp28 = (temp26) * Bza[iz][iy][_PB_CXM];
	  temp29 = temp24 - temp28;
	  Hz[iz][iy][_PB_CXM]= temp27 + temp29;
	  Bza[iz][iy][_PB_CXM] = tmp1;
	  for (ix = 0; ix < _PB_CXM; ix++)
	    {
	      temp62 = ix+1;
	      temp30 = Ey[iz][_PB_CYM][temp62] - Ey[iz][_PB_CYM][ix];
	      temp31 = Ex[iz][_PB_CYM][ix] - Ax[iz][ix];
	      clf[iz][iy] = temp31 + temp30;
	      temp31 = cymh[_PB_CYM] / cyph[iy];
	      temp32 = ch / cyph[iy];
	      temp33 = (temp32) * clf[iz][iy];
	      temp34 = (temp31) * Bza[iz][iy][ix];
	      tmp1 = temp34 - temp33;
	      temp35 = cxmh[ix] / cxph[ix];
	      temp36 = (temp35) * Hz[iz][_PB_CYM][ix];
	      temp37 = czp[iz] / cxph[ix];
	      temp38 = mui * temp37;
	      temp39 = (temp38) * tmp1;
	      temp40 = mui * czm[iz];
	      temp41 = temp40 / cxph[ix];
	      temp42 = (temp41) * Bza[iz][_PB_CYM][ix];
	      temp43 = temp36 + temp39;
	      Hz[iz][_PB_CYM][ix] = temp43 - temp42;
	      Bza[iz][_PB_CYM][ix] = tmp1;
	    }
	    temp44 = Ax[iz][_PB_CXM] + Ry[iz][_PB_CYM];
	    temp45 = temp44 - Ey[iz][_PB_CYM][_PB_CXM];
	  clf[iz][iy] = Ex[iz][_PB_CYM][_PB_CXM] - temp45;
	  temp46 = cymh[_PB_CYM] / cyph[_PB_CYM];
	  temp47 = (temp46) * Bza[iz][_PB_CYM][_PB_CXM];
	  temp48 = ch / cyph[_PB_CYM];
	  temp49 = (temp48) * clf[iz][iy];
	  tmp1 = temp47 - temp49;
	  temp50 = cxmh[_PB_CXM] / cxph[_PB_CXM];
	  temp51 = (temp50) * Hz[iz][_PB_CYM][_PB_CXM];
	  temp52 = czp[iz] / cxph[_PB_CXM];
	  temp53 = mui * temp52;
	  temp54 = (temp53) * tmp1;
	  temp55 = czm[iz] / cxph[_PB_CXM];
	  temp56 = mui * temp55;
	  temp57 = (temp56) * Bza[iz][_PB_CYM][_PB_CXM];
	  temp58 = temp54 - temp57;
	  Hz[iz][_PB_CYM][_PB_CXM] = temp51 + temp58;
	  Bza[iz][_PB_CYM][_PB_CXM] = tmp1;
	}
    }
#pragma endscop

}


int main(int argc, char** argv)
{
  /* Retrieve problem size. */
  int cz = CZ;
  int cym = CYM;
  int cxm = CXM;

  /* Variable declaration/allocation. */
  DATA_TYPE mui;
  DATA_TYPE ch;
  POLYBENCH_2D_ARRAY_DECL(Ax,DATA_TYPE,CZ+1,CYM+1,cz+1,cym+1);
  POLYBENCH_2D_ARRAY_DECL(Ry,DATA_TYPE,CZ+1,CYM+1,cz+1,cym+1);
  POLYBENCH_2D_ARRAY_DECL(clf,DATA_TYPE,CYM+1,CXM+1,cym+1,cxm+1);
  POLYBENCH_2D_ARRAY_DECL(tmp,DATA_TYPE,CYM+1,CXM+1,cym+1,cxm+1);
  POLYBENCH_3D_ARRAY_DECL(Bza,DATA_TYPE,CZ+1,CYM+1,CXM+1,cz+1,cym+1,cxm+1);
  POLYBENCH_3D_ARRAY_DECL(Ex,DATA_TYPE,CZ+1,CYM+1,CXM+1,cz+1,cym+1,cxm+1);
  POLYBENCH_3D_ARRAY_DECL(Ey,DATA_TYPE,CZ+1,CYM+1,CXM+1,cz+1,cym+1,cxm+1);
  POLYBENCH_3D_ARRAY_DECL(Hz,DATA_TYPE,CZ+1,CYM+1,CXM+1,cz+1,cym+1,cxm+1);
  POLYBENCH_1D_ARRAY_DECL(czm,DATA_TYPE,CZ+1,cz+1);
  POLYBENCH_1D_ARRAY_DECL(czp,DATA_TYPE,CZ+1,cz+1);
  POLYBENCH_1D_ARRAY_DECL(cxmh,DATA_TYPE,CXM+1,cxm+1);
  POLYBENCH_1D_ARRAY_DECL(cxph,DATA_TYPE,CXM+1,cxm+1);
  POLYBENCH_1D_ARRAY_DECL(cymh,DATA_TYPE,CYM+1,cym+1);
  POLYBENCH_1D_ARRAY_DECL(cyph,DATA_TYPE,CYM+1,cym+1);

  /* Initialize array(s). */
  init_array (cz, cxm, cym, &mui, &ch,
  	      POLYBENCH_ARRAY(Ax),
  	      POLYBENCH_ARRAY(Ry),
  	      POLYBENCH_ARRAY(Ex),
  	      POLYBENCH_ARRAY(Ey),
  	      POLYBENCH_ARRAY(Hz),
  	      POLYBENCH_ARRAY(czm),
  	      POLYBENCH_ARRAY(czp),
  	      POLYBENCH_ARRAY(cxmh),
  	      POLYBENCH_ARRAY(cxph),
  	      POLYBENCH_ARRAY(cymh),
  	      POLYBENCH_ARRAY(cyph));

  /* Start timer. */
  polybench_start_instruments;

  /* Run kernel. */
  kernel_fdtd_apml (cz, cxm, cym, mui, ch,
  		    POLYBENCH_ARRAY(Ax),
  		    POLYBENCH_ARRAY(Ry),
  		    POLYBENCH_ARRAY(clf),
  		    POLYBENCH_ARRAY(tmp),
  		    POLYBENCH_ARRAY(Bza),
  		    POLYBENCH_ARRAY(Ex),
  		    POLYBENCH_ARRAY(Ey),
  		    POLYBENCH_ARRAY(Hz),
  		    POLYBENCH_ARRAY(czm),
  		    POLYBENCH_ARRAY(czp),
  		    POLYBENCH_ARRAY(cxmh),
  		    POLYBENCH_ARRAY(cxph),
  		    POLYBENCH_ARRAY(cymh),
  		    POLYBENCH_ARRAY(cyph));

  /* Stop and print timer. */
  polybench_stop_instruments;
  polybench_print_instruments;

  /* Prevent dead-code elimination. All live-out data must be printed
     by the function call in argument. */
  polybench_prevent_dce(print_array(cz, cxm, cym,
  				    POLYBENCH_ARRAY(Bza),
  				    POLYBENCH_ARRAY(Ex),
  				    POLYBENCH_ARRAY(Ey),
  				    POLYBENCH_ARRAY(Hz)));

  /* Be clean. */
  POLYBENCH_FREE_ARRAY(Ax);
  POLYBENCH_FREE_ARRAY(Ry);
  POLYBENCH_FREE_ARRAY(clf);
  POLYBENCH_FREE_ARRAY(tmp);
  POLYBENCH_FREE_ARRAY(Bza);
  POLYBENCH_FREE_ARRAY(Ex);
  POLYBENCH_FREE_ARRAY(Ey);
  POLYBENCH_FREE_ARRAY(Hz);
  POLYBENCH_FREE_ARRAY(czm);
  POLYBENCH_FREE_ARRAY(czp);
  POLYBENCH_FREE_ARRAY(cxmh);
  POLYBENCH_FREE_ARRAY(cxph);
  POLYBENCH_FREE_ARRAY(cymh);
  POLYBENCH_FREE_ARRAY(cyph);

  return 0;
}

