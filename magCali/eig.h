/*
 * eig.h
 *
 * Code generation for function 'eig'
 *
 * C source code generated on: Wed Oct 22 22:07:03 2014
 *
 */

#ifndef __EIG_H__
#define __EIG_H__
/* Include files */
#include <math.h>
#include <stddef.h>
#include <stdlib.h>
#include <string.h>
#include "rt_nonfinite.h"

#include "rtwtypes.h"
#include "magCali_types.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
extern void b_eig(const creal_T A[9], creal_T V[9], creal_T D[9]);
extern creal_T b_eml_div(const creal_T x, const creal_T y);
extern void eig(const real_T A[36], creal_T V[36], creal_T D[36]);
extern creal_T eml_div(const creal_T x, real_T y);
#endif
/* End of code generation (eig.h) */
