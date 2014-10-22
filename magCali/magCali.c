/*
 * magCali.c
 *
 * Code generation for function 'magCali'
 *
 * C source code generated on: Wed Oct 22 22:07:03 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "magCali.h"
#include "norm.h"
#include "eig.h"
#include "sqrt.h"
#include "inv.h"
#include "chol.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */

/* Function Definitions */
void magCali(const real_T data[6000], creal_T A_i[9], creal_T B[3])
{
  real_T D[20000];
  int32_T i;
  int32_T i0;
  real_T S[100];
  int32_T ind;
  real_T L_S22[16];
  real_T L_S22_i[16];
  real_T b_L_S22[16];
  real_T S22_i[16];
  real_T A[36];
  real_T b_S[24];
  real_T c_S[36];
  real_T norm_v1;
  static const real_T b_A[36] = { 0.0, 0.5, 0.5, 0.0, 0.0, 0.0, 0.5, 0.0, 0.5,
    0.0, 0.0, 0.0, 0.5, 0.5, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.25, 0.0, 0.0,
    0.0, 0.0, 0.0, 0.0, -0.25, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, -0.25 };

  creal_T L[36];
  creal_T V[36];
  creal_T b_max;
  creal_T v1[6];
  real_T norm_VV2;
  creal_T b_S22_i[24];
  creal_T v2[4];
  creal_T v[100];
  creal_T Q[9];
  creal_T L_Q[9];
  creal_T L_Q_i[9];
  creal_T b_L_Q[9];
  creal_T b_v[3];
  creal_T QB[3];
  creal_T b_Q[9];
  real_T norm_VV3;
  creal_T b_L_Q_i[9];
  real_T ar;
  real_T ai;

  /*  magnetometer calibration */
  /*   */
  /*  Jungtaek Kim */
  /*  jungtaek.kim@jt-inc.net */
  /*   */
  memset(&D[0], 0, 20000U * sizeof(real_T));
  for (i = 0; i < 2000; i++) {
    D[i] = data[i] * data[i];
    D[2000 + i] = data[2000 + i] * data[2000 + i];
    D[4000 + i] = data[4000 + i] * data[4000 + i];
    D[6000 + i] = 2.0 * data[2000 + i] * data[4000 + i];
    D[8000 + i] = 2.0 * data[i] * data[4000 + i];
    D[10000 + i] = 2.0 * data[i] * data[2000 + i];
    D[12000 + i] = 2.0 * data[i];
    D[14000 + i] = 2.0 * data[2000 + i];
    D[16000 + i] = 2.0 * data[4000 + i];
    D[18000 + i] = 1.0;
  }

  for (i0 = 0; i0 < 10; i0++) {
    for (i = 0; i < 10; i++) {
      S[i0 + 10 * i] = 0.0;
      for (ind = 0; ind < 2000; ind++) {
        S[i0 + 10 * i] += D[ind + 2000 * i0] * D[ind + 2000 * i];
      }
    }
  }

  for (i0 = 0; i0 < 4; i0++) {
    for (i = 0; i < 4; i++) {
      L_S22[i + (i0 << 2)] = S[(i + 10 * (6 + i0)) + 6];
    }
  }

  chol(L_S22);
  inv(L_S22, L_S22_i);
  for (i0 = 0; i0 < 4; i0++) {
    for (i = 0; i < 4; i++) {
      b_L_S22[i + (i0 << 2)] = L_S22[i0 + (i << 2)];
    }
  }

  inv(b_L_S22, L_S22);
  for (i0 = 0; i0 < 4; i0++) {
    for (i = 0; i < 4; i++) {
      S22_i[i0 + (i << 2)] = 0.0;
      for (ind = 0; ind < 4; ind++) {
        S22_i[i0 + (i << 2)] += L_S22[i0 + (ind << 2)] * L_S22_i[ind + (i << 2)];
      }
    }
  }

  /*      S22_i = pinv(S22);%I4\S22; */
  for (i0 = 0; i0 < 6; i0++) {
    for (i = 0; i < 4; i++) {
      b_S[i0 + 6 * i] = 0.0;
      for (ind = 0; ind < 4; ind++) {
        b_S[i0 + 6 * i] += S[i0 + 10 * (6 + ind)] * S22_i[ind + (i << 2)];
      }
    }
  }

  for (i0 = 0; i0 < 6; i0++) {
    for (i = 0; i < 6; i++) {
      norm_v1 = 0.0;
      for (ind = 0; ind < 4; ind++) {
        norm_v1 += b_S[i0 + 6 * ind] * S[(ind + 10 * i) + 6];
      }

      c_S[i0 + 6 * i] = S[i0 + 10 * i] - norm_v1;
    }
  }

  for (i0 = 0; i0 < 6; i0++) {
    for (i = 0; i < 6; i++) {
      A[i0 + 6 * i] = 0.0;
      for (ind = 0; ind < 6; ind++) {
        A[i0 + 6 * i] += b_A[i0 + 6 * ind] * c_S[ind + 6 * i];
      }
    }
  }

  eig(A, V, L);
  b_max = L[0];
  ind = 0;
  for (i = 0; i < 6; i++) {
    if (b_max.re < L[i + 6 * i].re) {
      b_max = L[i + 6 * i];
      ind = i;
    }
  }

  memcpy(&v1[0], &V[6 * ind], 6U * sizeof(creal_T));
  if (V[6 * ind].re < 0.0) {
    for (i0 = 0; i0 < 6; i0++) {
      v1[i0].re = -V[i0 + 6 * ind].re;
      v1[i0].im = -V[i0 + 6 * ind].im;
    }
  }

  norm_v1 = norm(v1);
  for (i0 = 0; i0 < 6; i0++) {
    norm_VV2 = v1[i0].im;
    if (v1[i0].im == 0.0) {
      v1[i0].re /= norm_v1;
      v1[i0].im = 0.0;
    } else if (v1[i0].re == 0.0) {
      v1[i0].re = 0.0;
      v1[i0].im = norm_VV2 / norm_v1;
    } else {
      v1[i0].re /= norm_v1;
      v1[i0].im = norm_VV2 / norm_v1;
    }
  }

  for (i0 = 0; i0 < 4; i0++) {
    for (i = 0; i < 6; i++) {
      norm_v1 = 0.0;
      for (ind = 0; ind < 4; ind++) {
        norm_v1 += S22_i[i0 + (ind << 2)] * S[(ind + 10 * i) + 6];
      }

      b_S22_i[i0 + (i << 2)].re = norm_v1;
      b_S22_i[i0 + (i << 2)].im = 0.0;
    }
  }

  for (i0 = 0; i0 < 4; i0++) {
    v2[i0].re = 0.0;
    v2[i0].im = 0.0;
    for (i = 0; i < 6; i++) {
      v2[i0].re += b_S22_i[i0 + (i << 2)].re * v1[i].re - 0.0 * v1[i].im;
      v2[i0].im += b_S22_i[i0 + (i << 2)].re * v1[i].im + 0.0 * v1[i].re;
    }
  }

  for (i0 = 0; i0 < 100; i0++) {
    v[i0].re = 0.0;
    v[i0].im = 0.0;
  }

  v[0] = v1[0];
  v[1] = v1[1];
  v[2] = v1[2];
  v[3] = v1[3];
  v[4] = v1[4];
  v[5] = v1[5];
  v[6].re = -v2[0].re;
  v[6].im = -v2[0].im;
  v[7].re = -v2[1].re;
  v[7].im = -v2[1].im;
  v[8].re = -v2[2].re;
  v[8].im = -v2[2].im;
  v[9].re = -v2[3].re;
  v[9].im = -v2[3].im;
  Q[0] = v[0];
  Q[3] = v[5];
  Q[6] = v[4];
  Q[1] = v[5];
  Q[4] = v[1];
  Q[7] = v[3];
  Q[2] = v[4];
  Q[5] = v[3];
  Q[8] = v[2];
  memcpy(&L_Q[0], &Q[0], 9U * sizeof(creal_T));
  b_chol(L_Q);
  b_inv(L_Q, L_Q_i);
  for (i0 = 0; i0 < 3; i0++) {
    for (i = 0; i < 3; i++) {
      b_L_Q[i + 3 * i0].re = L_Q[i0 + 3 * i].re;
      b_L_Q[i + 3 * i0].im = -L_Q[i0 + 3 * i].im;
    }
  }

  b_inv(b_L_Q, L_Q);

  /*      Q_i = pinv(Q);%I3\Q; */
  for (i0 = 0; i0 < 3; i0++) {
    for (i = 0; i < 3; i++) {
      b_L_Q[i0 + 3 * i].re = 0.0;
      b_L_Q[i0 + 3 * i].im = 0.0;
      for (ind = 0; ind < 3; ind++) {
        b_L_Q[i0 + 3 * i].re += L_Q[i0 + 3 * ind].re * L_Q_i[ind + 3 * i].re -
          L_Q[i0 + 3 * ind].im * L_Q_i[ind + 3 * i].im;
        b_L_Q[i0 + 3 * i].im += L_Q[i0 + 3 * ind].re * L_Q_i[ind + 3 * i].im +
          L_Q[i0 + 3 * ind].im * L_Q_i[ind + 3 * i].re;
      }
    }
  }

  b_v[0] = v[6];
  b_v[1] = v[7];
  b_v[2] = v[8];
  for (i0 = 0; i0 < 3; i0++) {
    norm_v1 = 0.0;
    norm_VV2 = 0.0;
    for (i = 0; i < 3; i++) {
      norm_v1 += b_L_Q[i0 + 3 * i].re * b_v[i].re - b_L_Q[i0 + 3 * i].im * b_v[i]
        .im;
      norm_VV2 += b_L_Q[i0 + 3 * i].re * b_v[i].im + b_L_Q[i0 + 3 * i].im *
        b_v[i].re;
    }

    B[i0].re = -norm_v1;
    B[i0].im = -norm_VV2;
  }

  for (i0 = 0; i0 < 3; i0++) {
    QB[i0].re = 0.0;
    QB[i0].im = 0.0;
    for (i = 0; i < 3; i++) {
      QB[i0].re += Q[i0 + 3 * i].re * B[i].re - Q[i0 + 3 * i].im * B[i].im;
      QB[i0].im += Q[i0 + 3 * i].re * B[i].im + Q[i0 + 3 * i].im * B[i].re;
    }

    b_v[i0].re = B[i0].re;
    b_v[i0].im = -B[i0].im;
  }

  b_max.re = 0.0;
  b_max.im = 0.0;
  for (i = 0; i < 3; i++) {
    b_max.re += b_v[i].re * QB[i].re - b_v[i].im * QB[i].im;
    b_max.im += b_v[i].re * QB[i].im + b_v[i].im * QB[i].re;
  }

  b_max.re -= v[9].re;
  b_max.im -= v[9].im;
  b_sqrt(&b_max);
  memcpy(&b_Q[0], &Q[0], 9U * sizeof(creal_T));
  b_eig(b_Q, L_Q, Q);
  norm_v1 = b_norm(*(creal_T (*)[3])&L_Q[0]);
  norm_VV2 = b_norm(*(creal_T (*)[3])&L_Q[3]);
  norm_VV3 = b_norm(*(creal_T (*)[3])&L_Q[6]);
  for (i0 = 0; i0 < 3; i0++) {
    if (L_Q[i0].im == 0.0) {
      L_Q_i[i0].re = L_Q[i0].re / norm_v1;
      L_Q_i[i0].im = 0.0;
    } else if (L_Q[i0].re == 0.0) {
      L_Q_i[i0].re = 0.0;
      L_Q_i[i0].im = L_Q[i0].im / norm_v1;
    } else {
      L_Q_i[i0].re = L_Q[i0].re / norm_v1;
      L_Q_i[i0].im = L_Q[i0].im / norm_v1;
    }

    if (L_Q[3 + i0].im == 0.0) {
      L_Q_i[3 + i0].re = L_Q[3 + i0].re / norm_VV2;
      L_Q_i[3 + i0].im = 0.0;
    } else if (L_Q[3 + i0].re == 0.0) {
      L_Q_i[3 + i0].re = 0.0;
      L_Q_i[3 + i0].im = L_Q[3 + i0].im / norm_VV2;
    } else {
      L_Q_i[3 + i0].re = L_Q[3 + i0].re / norm_VV2;
      L_Q_i[3 + i0].im = L_Q[3 + i0].im / norm_VV2;
    }

    if (L_Q[6 + i0].im == 0.0) {
      L_Q_i[6 + i0].re = L_Q[6 + i0].re / norm_VV3;
      L_Q_i[6 + i0].im = 0.0;
    } else if (L_Q[6 + i0].re == 0.0) {
      L_Q_i[6 + i0].re = 0.0;
      L_Q_i[6 + i0].im = L_Q[6 + i0].im / norm_VV3;
    } else {
      L_Q_i[6 + i0].re = L_Q[6 + i0].re / norm_VV3;
      L_Q_i[6 + i0].im = L_Q[6 + i0].im / norm_VV3;
    }

    for (i = 0; i < 3; i++) {
      b_L_Q[i0 + 3 * i].re = 0.0;
      b_L_Q[i0 + 3 * i].im = 0.0;
      for (ind = 0; ind < 3; ind++) {
        b_L_Q[i0 + 3 * i].re += L_Q_i[i0 + 3 * ind].re * Q[ind + 3 * i].re -
          L_Q_i[i0 + 3 * ind].im * Q[ind + 3 * i].im;
        b_L_Q[i0 + 3 * i].im += L_Q_i[i0 + 3 * ind].re * Q[ind + 3 * i].im +
          L_Q_i[i0 + 3 * ind].im * Q[ind + 3 * i].re;
      }
    }
  }

  for (i0 = 0; i0 < 3; i0++) {
    for (i = 0; i < 3; i++) {
      b_L_Q_i[i0 + 3 * i].re = 0.0;
      b_L_Q_i[i0 + 3 * i].im = 0.0;
      for (ind = 0; ind < 3; ind++) {
        b_L_Q_i[i0 + 3 * i].re += b_L_Q[i0 + 3 * ind].re * L_Q_i[i + 3 * ind].re
          - b_L_Q[i0 + 3 * ind].im * -L_Q_i[i + 3 * ind].im;
        b_L_Q_i[i0 + 3 * i].im += b_L_Q[i0 + 3 * ind].re * -L_Q_i[i + 3 * ind].
          im + b_L_Q[i0 + 3 * ind].im * L_Q_i[i + 3 * ind].re;
      }
    }
  }

  for (i0 = 0; i0 < 3; i0++) {
    for (i = 0; i < 3; i++) {
      ar = b_L_Q_i[i + 3 * i0].re * 0.569;
      ai = b_L_Q_i[i + 3 * i0].im * 0.569;
      if (b_max.im == 0.0) {
        if (ai == 0.0) {
          A_i[i + 3 * i0].re = ar / b_max.re;
          A_i[i + 3 * i0].im = 0.0;
        } else if (ar == 0.0) {
          A_i[i + 3 * i0].re = 0.0;
          A_i[i + 3 * i0].im = ai / b_max.re;
        } else {
          A_i[i + 3 * i0].re = ar / b_max.re;
          A_i[i + 3 * i0].im = ai / b_max.re;
        }
      } else if (b_max.re == 0.0) {
        if (ar == 0.0) {
          A_i[i + 3 * i0].re = ai / b_max.im;
          A_i[i + 3 * i0].im = 0.0;
        } else if (ai == 0.0) {
          A_i[i + 3 * i0].re = 0.0;
          A_i[i + 3 * i0].im = -(ar / b_max.im);
        } else {
          A_i[i + 3 * i0].re = ai / b_max.im;
          A_i[i + 3 * i0].im = -(ar / b_max.im);
        }
      } else {
        norm_VV3 = fabs(b_max.re);
        norm_v1 = fabs(b_max.im);
        if (norm_VV3 > norm_v1) {
          norm_v1 = b_max.im / b_max.re;
          norm_VV2 = b_max.re + norm_v1 * b_max.im;
          A_i[i + 3 * i0].re = (ar + norm_v1 * ai) / norm_VV2;
          A_i[i + 3 * i0].im = (ai - norm_v1 * ar) / norm_VV2;
        } else if (norm_v1 == norm_VV3) {
          norm_v1 = b_max.re > 0.0 ? 0.5 : -0.5;
          norm_VV2 = b_max.im > 0.0 ? 0.5 : -0.5;
          A_i[i + 3 * i0].re = (ar * norm_v1 + ai * norm_VV2) / norm_VV3;
          A_i[i + 3 * i0].im = (ai * norm_v1 - ar * norm_VV2) / norm_VV3;
        } else {
          norm_v1 = b_max.re / b_max.im;
          norm_VV2 = b_max.im + norm_v1 * b_max.re;
          A_i[i + 3 * i0].re = (norm_v1 * ar + ai) / norm_VV2;
          A_i[i + 3 * i0].im = (norm_v1 * ai - ar) / norm_VV2;
        }
      }
    }
  }
}

/* End of code generation (magCali.c) */
