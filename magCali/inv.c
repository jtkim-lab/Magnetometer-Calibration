/*
 * inv.c
 *
 * Code generation for function 'inv'
 *
 * C source code generated on: Wed Oct 22 22:07:03 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "magCali.h"
#include "inv.h"
#include "eig.h"
#include "magCali_rtwutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static void inv3x3(creal_T x[9], creal_T y[9]);
static void invNxN(const real_T x[16], real_T y[16]);

/* Function Definitions */
static void inv3x3(creal_T x[9], creal_T y[9])
{
  int32_T p1;
  int32_T p2;
  int32_T p3;
  real_T absx11;
  real_T absx21;
  real_T absx31;
  creal_T t1;
  int32_T itmp;
  creal_T b_x;
  creal_T t2;
  p1 = 0;
  p2 = 3;
  p3 = 6;
  absx11 = rt_hypotd_snf(fabs(x[0].re), fabs(x[0].im));
  absx21 = rt_hypotd_snf(fabs(x[1].re), fabs(x[1].im));
  absx31 = rt_hypotd_snf(fabs(x[2].re), fabs(x[2].im));
  if ((absx21 > absx11) && (absx21 > absx31)) {
    p1 = 3;
    p2 = 0;
    t1 = x[0];
    x[0] = x[1];
    x[1] = t1;
    t1 = x[3];
    x[3] = x[4];
    x[4] = t1;
    t1 = x[6];
    x[6] = x[7];
    x[7] = t1;
  } else {
    if (absx31 > absx11) {
      p1 = 6;
      p3 = 0;
      t1 = x[0];
      x[0] = x[2];
      x[2] = t1;
      t1 = x[3];
      x[3] = x[5];
      x[5] = t1;
      t1 = x[6];
      x[6] = x[8];
      x[8] = t1;
    }
  }

  x[1] = b_eml_div(x[1], x[0]);
  x[2] = b_eml_div(x[2], x[0]);
  absx11 = x[1].re * x[3].im + x[1].im * x[3].re;
  x[4].re -= x[1].re * x[3].re - x[1].im * x[3].im;
  x[4].im -= absx11;
  absx11 = x[2].re * x[3].im + x[2].im * x[3].re;
  x[5].re -= x[2].re * x[3].re - x[2].im * x[3].im;
  x[5].im -= absx11;
  absx11 = x[1].re * x[6].im + x[1].im * x[6].re;
  x[7].re -= x[1].re * x[6].re - x[1].im * x[6].im;
  x[7].im -= absx11;
  absx11 = x[2].re * x[6].im + x[2].im * x[6].re;
  x[8].re -= x[2].re * x[6].re - x[2].im * x[6].im;
  x[8].im -= absx11;
  if (rt_hypotd_snf(fabs(x[5].re), fabs(x[5].im)) > rt_hypotd_snf(fabs(x[4].re),
       fabs(x[4].im))) {
    itmp = p2;
    p2 = p3;
    p3 = itmp;
    t1 = x[1];
    x[1] = x[2];
    x[2] = t1;
    t1 = x[4];
    x[4] = x[5];
    x[5] = t1;
    t1 = x[7];
    x[7] = x[8];
    x[8] = t1;
  }

  x[5] = b_eml_div(x[5], x[4]);
  absx11 = x[5].re * x[7].im + x[5].im * x[7].re;
  x[8].re -= x[5].re * x[7].re - x[5].im * x[7].im;
  x[8].im -= absx11;
  b_x.re = (x[5].re * x[1].re - x[5].im * x[1].im) - x[2].re;
  b_x.im = (x[5].re * x[1].im + x[5].im * x[1].re) - x[2].im;
  t1 = b_eml_div(b_x, x[8]);
  b_x.re = -(x[1].re + (x[7].re * t1.re - x[7].im * t1.im));
  b_x.im = -(x[1].im + (x[7].re * t1.im + x[7].im * t1.re));
  t2 = b_eml_div(b_x, x[4]);
  b_x.re = (1.0 - (x[3].re * t2.re - x[3].im * t2.im)) - (x[6].re * t1.re - x[6]
    .im * t1.im);
  b_x.im = (0.0 - (x[3].re * t2.im + x[3].im * t2.re)) - (x[6].re * t1.im + x[6]
    .im * t1.re);
  y[p1] = b_eml_div(b_x, x[0]);
  y[p1 + 1] = t2;
  y[p1 + 2] = t1;
  b_x.re = -x[5].re;
  b_x.im = -x[5].im;
  t1 = b_eml_div(b_x, x[8]);
  b_x.re = 1.0 - (x[7].re * t1.re - x[7].im * t1.im);
  b_x.im = 0.0 - (x[7].re * t1.im + x[7].im * t1.re);
  t2 = b_eml_div(b_x, x[4]);
  b_x.re = -((x[3].re * t2.re - x[3].im * t2.im) + (x[6].re * t1.re - x[6].im *
              t1.im));
  b_x.im = -((x[3].re * t2.im + x[3].im * t2.re) + (x[6].re * t1.im + x[6].im *
              t1.re));
  y[p2] = b_eml_div(b_x, x[0]);
  y[p2 + 1] = t2;
  y[p2 + 2] = t1;
  if (x[8].im == 0.0) {
    t1.re = 1.0 / x[8].re;
    t1.im = 0.0;
  } else if (x[8].re == 0.0) {
    t1.re = 0.0;
    t1.im = -(1.0 / x[8].im);
  } else {
    absx31 = fabs(x[8].re);
    absx11 = fabs(x[8].im);
    if (absx31 > absx11) {
      absx11 = x[8].im / x[8].re;
      absx21 = x[8].re + absx11 * x[8].im;
      t1.re = (1.0 + absx11 * 0.0) / absx21;
      t1.im = (0.0 - absx11) / absx21;
    } else if (absx11 == absx31) {
      absx11 = x[8].re > 0.0 ? 0.5 : -0.5;
      absx21 = x[8].im > 0.0 ? 0.5 : -0.5;
      t1.re = absx11 / absx31;
      t1.im = (-0.0 - absx21) / absx31;
    } else {
      absx11 = x[8].re / x[8].im;
      absx21 = x[8].im + absx11 * x[8].re;
      t1.re = absx11 / absx21;
      t1.im = (absx11 * 0.0 - 1.0) / absx21;
    }
  }

  b_x.re = -x[7].re * t1.re - -x[7].im * t1.im;
  b_x.im = -x[7].re * t1.im + -x[7].im * t1.re;
  t2 = b_eml_div(b_x, x[4]);
  b_x.re = -((x[3].re * t2.re - x[3].im * t2.im) + (x[6].re * t1.re - x[6].im *
              t1.im));
  b_x.im = -((x[3].re * t2.im + x[3].im * t2.re) + (x[6].re * t1.im + x[6].im *
              t1.re));
  y[p3] = b_eml_div(b_x, x[0]);
  y[p3 + 1] = t2;
  y[p3 + 2] = t1;
}

static void invNxN(const real_T x[16], real_T y[16])
{
  real_T A[16];
  int32_T i1;
  int8_T ipiv[4];
  int32_T j;
  int32_T jj;
  int32_T jp1j;
  int32_T pipk;
  int32_T ix;
  real_T smax;
  int32_T jA;
  real_T s;
  int32_T i;
  int8_T p[4];
  for (i1 = 0; i1 < 16; i1++) {
    y[i1] = 0.0;
    A[i1] = x[i1];
  }

  for (i1 = 0; i1 < 4; i1++) {
    ipiv[i1] = (int8_T)(1 + i1);
  }

  for (j = 0; j < 3; j++) {
    jj = j * 5;
    jp1j = jj + 2;
    pipk = 1;
    ix = jj;
    smax = fabs(A[jj]);
    for (jA = 2; jA <= 4 - j; jA++) {
      ix++;
      s = fabs(A[ix]);
      if (s > smax) {
        pipk = jA;
        smax = s;
      }
    }

    if (A[(jj + pipk) - 1] != 0.0) {
      if (pipk - 1 != 0) {
        ipiv[j] = (int8_T)(j + pipk);
        ix = j;
        pipk = (j + pipk) - 1;
        for (jA = 0; jA < 4; jA++) {
          smax = A[ix];
          A[ix] = A[pipk];
          A[pipk] = smax;
          ix += 4;
          pipk += 4;
        }
      }

      i1 = (jp1j - j) + 2;
      for (i = jp1j; i <= i1; i++) {
        A[i - 1] /= A[jj];
      }
    }

    jA = jj + 5;
    pipk = jj + 4;
    for (jj = 1; jj <= 3 - j; jj++) {
      smax = A[pipk];
      if (A[pipk] != 0.0) {
        ix = jp1j;
        i1 = (jA - j) + 3;
        for (i = jA; i + 1 <= i1; i++) {
          A[i] += A[ix - 1] * -smax;
          ix++;
        }
      }

      pipk += 4;
      jA += 4;
    }
  }

  for (i1 = 0; i1 < 4; i1++) {
    p[i1] = (int8_T)(1 + i1);
  }

  for (jA = 0; jA < 3; jA++) {
    if (ipiv[jA] > 1 + jA) {
      pipk = p[ipiv[jA] - 1];
      p[ipiv[jA] - 1] = p[jA];
      p[jA] = (int8_T)pipk;
    }
  }

  for (jA = 0; jA < 4; jA++) {
    y[jA + ((p[jA] - 1) << 2)] = 1.0;
    for (j = jA; j + 1 < 5; j++) {
      if (y[j + ((p[jA] - 1) << 2)] != 0.0) {
        for (i = j + 1; i + 1 < 5; i++) {
          y[i + ((p[jA] - 1) << 2)] -= y[j + ((p[jA] - 1) << 2)] * A[i + (j << 2)];
        }
      }
    }
  }

  for (j = 0; j < 4; j++) {
    pipk = j << 2;
    for (jA = 3; jA > -1; jA += -1) {
      jj = jA << 2;
      if (y[jA + pipk] != 0.0) {
        y[jA + pipk] /= A[jA + jj];
        for (i = 0; i + 1 <= jA; i++) {
          y[i + pipk] -= y[jA + pipk] * A[i + jj];
        }
      }
    }
  }
}

void b_inv(const creal_T x[9], creal_T y[9])
{
  creal_T b_x[9];
  memcpy(&b_x[0], &x[0], 9U * sizeof(creal_T));
  inv3x3(b_x, y);
}

void inv(const real_T x[16], real_T y[16])
{
  invNxN(x, y);
}

/* End of code generation (inv.c) */
