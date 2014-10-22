/*
 * chol.c
 *
 * Code generation for function 'chol'
 *
 * C source code generated on: Wed Oct 22 22:07:03 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "magCali.h"
#include "chol.h"
#include "eig.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void b_chol(creal_T A[9])
{
  int32_T j;
  int32_T exitg1;
  int32_T info;
  boolean_T exitg2;
  int32_T jj;
  creal_T c;
  int32_T ix;
  int32_T iy;
  int32_T k;
  real_T ajj;
  int32_T jp1;
  int32_T jp1j;
  int32_T ia;
  real_T A_re;
  real_T A_im;
  static const creal_T dc1 = { 1.0, 0.0 };

  j = 0;
  do {
    exitg1 = 0;
    if (j + 1 < 4) {
      if (A[j + 3 * j].im != 0.0) {
        exitg1 = 1;
      } else {
        j++;
      }
    } else {
      info = 0;
      j = 0;
      exitg2 = FALSE;
      while ((exitg2 == 0U) && (j + 1 < 4)) {
        jj = j + j * 3;
        c.re = 0.0;
        c.im = 0.0;
        if (j < 1) {
        } else {
          ix = j;
          iy = j;
          for (k = 1; k <= j; k++) {
            c.re += A[ix].re * A[iy].re + A[ix].im * A[iy].im;
            c.im += A[ix].re * A[iy].im - A[ix].im * A[iy].re;
            ix += 3;
            iy += 3;
          }
        }

        ajj = A[jj].re - c.re;
        if ((A[jj].im == 0.0) && (ajj > 0.0)) {
          ajj = sqrt(ajj);
          A[jj].re = ajj;
          A[jj].im = 0.0;
          if (j + 1 < 3) {
            jp1 = j + 2;
            jp1j = jj + 1;
            for (k = 0; k + 1 <= j; k++) {
              A[j + 3 * k].im = -A[j + 3 * k].im;
            }

            if (j == 0) {
            } else {
              ix = j + 1;
              jj = jp1 + 3 * (j - 1);
              while (jp1 <= jj) {
                c.re = -A[ix - 1].re - 0.0 * A[ix - 1].im;
                c.im = -A[ix - 1].im + 0.0 * A[ix - 1].re;
                iy = jp1j;
                k = (jp1 - j) + 1;
                for (ia = jp1; ia <= k; ia++) {
                  A_re = A[ia - 1].re * c.re - A[ia - 1].im * c.im;
                  A_im = A[ia - 1].re * c.im + A[ia - 1].im * c.re;
                  A[iy].re += A_re;
                  A[iy].im += A_im;
                  iy++;
                }

                ix += 3;
                jp1 += 3;
              }
            }

            for (k = 0; k + 1 <= j; k++) {
              A[j + 3 * k].im = -A[j + 3 * k].im;
            }

            c = eml_div(dc1, ajj);
            jj = (jp1j - j) + 2;
            while (jp1j + 1 <= jj) {
              A_re = A[jp1j].re;
              A_im = A[jp1j].im;
              A[jp1j].re = c.re * A_re - c.im * A_im;
              A[jp1j].im = c.re * A_im + c.im * A_re;
              jp1j++;
            }
          }

          j++;
        } else {
          A[jj].re = ajj;
          A[jj].im = 0.0;
          info = j + 1;
          exitg2 = TRUE;
        }
      }

      if (info == 0) {
        jj = 3;
      } else {
        jj = info - 1;
      }

      for (j = 1; j + 1 <= jj; j++) {
        for (k = 1; k <= j; k++) {
          A[(k + 3 * j) - 1].re = 0.0;
          A[(k + 3 * j) - 1].im = 0.0;
        }
      }

      exitg1 = 1;
    }
  } while (exitg1 == 0U);
}

void chol(real_T A[16])
{
  int32_T info;
  int32_T j;
  boolean_T exitg1;
  int32_T jj;
  real_T ajj;
  int32_T ix;
  int32_T iy;
  int32_T jmax;
  int32_T jp1;
  real_T c;
  int32_T i2;
  int32_T ia;
  info = 0;
  j = 0;
  exitg1 = FALSE;
  while ((exitg1 == 0U) && (j + 1 < 5)) {
    jj = j + (j << 2);
    ajj = 0.0;
    if (j < 1) {
    } else {
      ix = j;
      iy = j;
      for (jmax = 1; jmax <= j; jmax++) {
        ajj += A[ix] * A[iy];
        ix += 4;
        iy += 4;
      }
    }

    ajj = A[jj] - ajj;
    if (ajj > 0.0) {
      ajj = sqrt(ajj);
      A[jj] = ajj;
      if (j + 1 < 4) {
        jp1 = j + 2;
        jmax = jj + 1;
        if (j == 0) {
        } else {
          ix = j + 1;
          jj = jp1 + ((j - 1) << 2);
          while (jp1 <= jj) {
            c = -A[ix - 1];
            iy = jmax;
            i2 = (jp1 - j) + 2;
            for (ia = jp1; ia <= i2; ia++) {
              A[iy] += A[ia - 1] * c;
              iy++;
            }

            ix += 4;
            jp1 += 4;
          }
        }

        ajj = 1.0 / ajj;
        jj = (jmax - j) + 3;
        while (jmax + 1 <= jj) {
          A[jmax] *= ajj;
          jmax++;
        }
      }

      j++;
    } else {
      A[jj] = ajj;
      info = j + 1;
      exitg1 = TRUE;
    }
  }

  if (info == 0) {
    jmax = 4;
  } else {
    jmax = info - 1;
  }

  for (j = 1; j + 1 <= jmax; j++) {
    for (jj = 1; jj <= j; jj++) {
      A[(jj + (j << 2)) - 1] = 0.0;
    }
  }
}

/* End of code generation (chol.c) */
