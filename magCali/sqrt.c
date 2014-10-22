/*
 * sqrt.c
 *
 * Code generation for function 'sqrt'
 *
 * C source code generated on: Wed Oct 22 22:07:03 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "magCali.h"
#include "sqrt.h"
#include "eig.h"
#include "magCali_rtwutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void b_sqrt(creal_T *x)
{
  real_T absxi;
  real_T absxr;
  if (x->im == 0.0) {
    if (x->re < 0.0) {
      absxi = 0.0;
      absxr = sqrt(fabs(x->re));
    } else {
      absxi = sqrt(x->re);
      absxr = 0.0;
    }
  } else if (x->re == 0.0) {
    if (x->im < 0.0) {
      absxi = sqrt(-x->im / 2.0);
      absxr = -absxi;
    } else {
      absxi = sqrt(x->im / 2.0);
      absxr = absxi;
    }
  } else if (rtIsNaN(x->re) || rtIsNaN(x->im)) {
    absxi = rtNaN;
    absxr = rtNaN;
  } else if (rtIsInf(x->im)) {
    absxi = rtInf;
    absxr = x->im;
  } else if (rtIsInf(x->re)) {
    if (x->re < 0.0) {
      absxi = 0.0;
      absxr = rtInf;
    } else {
      absxi = rtInf;
      absxr = 0.0;
    }
  } else {
    absxr = fabs(x->re);
    absxi = fabs(x->im);
    if ((absxr > 4.4942328371557893E+307) || (absxi > 4.4942328371557893E+307))
    {
      absxr *= 0.5;
      absxi *= 0.5;
      absxi = rt_hypotd_snf(absxr, absxi);
      if (absxi > absxr) {
        absxi = sqrt(absxi) * sqrt(1.0 + absxr / absxi);
      } else {
        absxi = sqrt(absxi) * 1.4142135623730951;
      }
    } else {
      absxi = sqrt((rt_hypotd_snf(absxr, absxi) + absxr) * 0.5);
    }

    if (x->re > 0.0) {
      absxr = 0.5 * (x->im / absxi);
    } else {
      if (x->im < 0.0) {
        absxr = -absxi;
      } else {
        absxr = absxi;
      }

      absxi = 0.5 * (x->im / absxr);
    }
  }

  x->re = absxi;
  x->im = absxr;
}

/* End of code generation (sqrt.c) */
