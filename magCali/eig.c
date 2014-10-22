/*
 * eig.c
 *
 * Code generation for function 'eig'
 *
 * C source code generated on: Wed Oct 22 22:07:03 2014
 *
 */

/* Include files */
#include "rt_nonfinite.h"
#include "magCali.h"
#include "eig.h"
#include "sqrt.h"
#include "magCali_rtwutil.h"

/* Type Definitions */

/* Named Constants */

/* Variable Declarations */

/* Variable Definitions */

/* Function Declarations */
static void b_eml_matlab_zggev(creal_T A[9], real_T *info, creal_T alpha1[3],
  creal_T beta1[3], creal_T V[9]);
static void b_eml_matlab_zhgeqz(creal_T A[9], int32_T ilo, int32_T ihi, creal_T
  Z[9], real_T *info, creal_T alpha1[3], creal_T beta1[3]);
static real_T b_eml_matlab_zlanhs(const creal_T A[9], int32_T ilo, int32_T ihi);
static void b_eml_matlab_zlartg(const creal_T f, const creal_T g, real_T *cs,
  creal_T *sn);
static void b_eml_matlab_ztgevc(const creal_T A[9], creal_T V[9]);
static int32_T div_s32_floor(int32_T numerator, int32_T denominator);
static void eml_matlab_zggbal(creal_T A[36], int32_T *ilo, int32_T *ihi, int32_T
  rscale[6]);
static void eml_matlab_zggev(creal_T A[36], real_T *info, creal_T alpha1[6],
  creal_T beta1[6], creal_T V[36]);
static void eml_matlab_zhgeqz(creal_T A[36], int32_T ilo, int32_T ihi, creal_T
  Z[36], real_T *info, creal_T alpha1[6], creal_T beta1[6]);
static real_T eml_matlab_zlangeM(const creal_T x[9]);
static real_T eml_matlab_zlanhs(const creal_T A[36], int32_T ilo, int32_T ihi);
static void eml_matlab_zlartg(const creal_T f, const creal_T g, real_T *cs,
  creal_T *sn, creal_T *r);
static void eml_matlab_ztgevc(const creal_T A[36], creal_T V[36]);

/* Function Definitions */
static void b_eml_matlab_zggev(creal_T A[9], real_T *info, creal_T alpha1[3],
  creal_T beta1[3], creal_T V[9])
{
  real_T anrm;
  int32_T i;
  int32_T ii;
  boolean_T ilascl;
  real_T anrmto;
  real_T cfromc;
  real_T ctoc;
  boolean_T notdone;
  real_T cfrom1;
  real_T cto1;
  real_T mul;
  creal_T b_A[9];
  int32_T rscale[3];
  int32_T ilo;
  int32_T ihi;
  int32_T exitg2;
  int32_T j;
  boolean_T exitg5;
  int32_T nzcount;
  int32_T jj;
  boolean_T exitg6;
  boolean_T guard2 = FALSE;
  creal_T atmp;
  int32_T exitg1;
  boolean_T exitg3;
  boolean_T exitg4;
  boolean_T guard1 = FALSE;
  static const int8_T iv1[9] = { 1, 0, 0, 0, 1, 0, 0, 0, 1 };

  creal_T s;
  *info = 0.0;
  anrm = eml_matlab_zlangeM(A);
  if (!((!rtIsInf(anrm)) && (!rtIsNaN(anrm)))) {
    for (i = 0; i < 3; i++) {
      alpha1[i].re = rtNaN;
      alpha1[i].im = 0.0;
      beta1[i].re = rtNaN;
      beta1[i].im = 0.0;
    }

    for (ii = 0; ii < 9; ii++) {
      V[ii].re = rtNaN;
      V[ii].im = 0.0;
    }
  } else {
    ilascl = FALSE;
    anrmto = anrm;
    if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
      anrmto = 6.7178761075670888E-139;
      ilascl = TRUE;
    } else {
      if (anrm > 1.4885657073574029E+138) {
        anrmto = 1.4885657073574029E+138;
        ilascl = TRUE;
      }
    }

    if (ilascl) {
      cfromc = anrm;
      ctoc = anrmto;
      notdone = TRUE;
      while (notdone) {
        cfrom1 = cfromc * 2.0041683600089728E-292;
        cto1 = ctoc / 4.9896007738368E+291;
        if ((fabs(cfrom1) > fabs(ctoc)) && (ctoc != 0.0)) {
          mul = 2.0041683600089728E-292;
          cfromc = cfrom1;
        } else if (fabs(cto1) > fabs(cfromc)) {
          mul = 4.9896007738368E+291;
          ctoc = cto1;
        } else {
          mul = ctoc / cfromc;
          notdone = FALSE;
        }

        for (ii = 0; ii < 9; ii++) {
          A[ii].re *= mul;
          A[ii].im *= mul;
        }
      }
    }

    memcpy(&b_A[0], &A[0], 9U * sizeof(creal_T));
    for (i = 0; i < 3; i++) {
      rscale[i] = 0;
    }

    ilo = 1;
    ihi = 3;
    do {
      exitg2 = 0;
      i = 0;
      j = 0;
      notdone = FALSE;
      ii = ihi;
      exitg5 = FALSE;
      while ((exitg5 == 0U) && (ii > 0)) {
        nzcount = 0;
        i = ii;
        j = ihi;
        jj = 1;
        exitg6 = FALSE;
        while ((exitg6 == 0U) && (jj <= ihi)) {
          guard2 = FALSE;
          if ((b_A[(ii + 3 * (jj - 1)) - 1].re != 0.0) || (b_A[(ii + 3 * (jj - 1))
               - 1].im != 0.0) || (ii == jj)) {
            if (nzcount == 0) {
              j = jj;
              nzcount = 1;
              guard2 = TRUE;
            } else {
              nzcount = 2;
              exitg6 = TRUE;
            }
          } else {
            guard2 = TRUE;
          }

          if (guard2 == TRUE) {
            jj++;
          }
        }

        if (nzcount < 2) {
          notdone = TRUE;
          exitg5 = TRUE;
        } else {
          ii--;
        }
      }

      if (!notdone) {
        exitg2 = 2;
      } else {
        if (i != ihi) {
          for (nzcount = 0; nzcount < 3; nzcount++) {
            atmp = b_A[(i + 3 * nzcount) - 1];
            b_A[(i + 3 * nzcount) - 1] = b_A[(ihi + 3 * nzcount) - 1];
            b_A[(ihi + 3 * nzcount) - 1] = atmp;
          }
        }

        if (j != ihi) {
          for (nzcount = 0; nzcount + 1 <= ihi; nzcount++) {
            atmp = b_A[nzcount + 3 * (j - 1)];
            b_A[nzcount + 3 * (j - 1)] = b_A[nzcount + 3 * (ihi - 1)];
            b_A[nzcount + 3 * (ihi - 1)] = atmp;
          }
        }

        rscale[ihi - 1] = j;
        ihi--;
        if (ihi == 1) {
          rscale[0] = 1;
          exitg2 = 1;
        }
      }
    } while (exitg2 == 0U);

    if (exitg2 == 1U) {
    } else {
      do {
        exitg1 = 0;
        i = 0;
        j = 0;
        notdone = FALSE;
        jj = ilo;
        exitg3 = FALSE;
        while ((exitg3 == 0U) && (jj <= ihi)) {
          nzcount = 0;
          i = ihi;
          j = jj;
          ii = ilo;
          exitg4 = FALSE;
          while ((exitg4 == 0U) && (ii <= ihi)) {
            guard1 = FALSE;
            if ((b_A[(ii + 3 * (jj - 1)) - 1].re != 0.0) || (b_A[(ii + 3 * (jj -
                   1)) - 1].im != 0.0) || (ii == jj)) {
              if (nzcount == 0) {
                i = ii;
                nzcount = 1;
                guard1 = TRUE;
              } else {
                nzcount = 2;
                exitg4 = TRUE;
              }
            } else {
              guard1 = TRUE;
            }

            if (guard1 == TRUE) {
              ii++;
            }
          }

          if (nzcount < 2) {
            notdone = TRUE;
            exitg3 = TRUE;
          } else {
            jj++;
          }
        }

        if (!notdone) {
          exitg1 = 1;
        } else {
          if (i != ilo) {
            for (nzcount = ilo - 1; nzcount + 1 < 4; nzcount++) {
              atmp = b_A[(i + 3 * nzcount) - 1];
              b_A[(i + 3 * nzcount) - 1] = b_A[(ilo + 3 * nzcount) - 1];
              b_A[(ilo + 3 * nzcount) - 1] = atmp;
            }
          }

          if (j != ilo) {
            for (nzcount = 0; nzcount + 1 <= ihi; nzcount++) {
              atmp = b_A[nzcount + 3 * (j - 1)];
              b_A[nzcount + 3 * (j - 1)] = b_A[nzcount + 3 * (ilo - 1)];
              b_A[nzcount + 3 * (ilo - 1)] = atmp;
            }
          }

          rscale[ilo - 1] = j;
          ilo++;
          if (ilo == ihi) {
            rscale[ilo - 1] = ilo;
            exitg1 = 1;
          }
        }
      } while (exitg1 == 0U);
    }

    for (ii = 0; ii < 9; ii++) {
      V[ii].re = (real_T)iv1[ii];
      V[ii].im = 0.0;
    }

    if (ihi < ilo + 2) {
    } else {
      ii = ilo;
      while (ii < 2) {
        eml_matlab_zlartg(b_A[1], b_A[2], &cfrom1, &s, &atmp);
        b_A[1] = atmp;
        b_A[2].re = 0.0;
        b_A[2].im = 0.0;
        for (j = 1; j + 1 < 4; j++) {
          atmp.re = cfrom1 * b_A[1 + 3 * j].re;
          atmp.im = cfrom1 * b_A[1 + 3 * j].im;
          cto1 = s.re * b_A[2 + 3 * j].re - s.im * b_A[2 + 3 * j].im;
          mul = s.re * b_A[2 + 3 * j].im + s.im * b_A[2 + 3 * j].re;
          cfromc = b_A[1 + 3 * j].im;
          ctoc = b_A[1 + 3 * j].re;
          b_A[2 + 3 * j].re = cfrom1 * b_A[2 + 3 * j].re - (s.re * b_A[1 + 3 * j]
            .re + s.im * b_A[1 + 3 * j].im);
          b_A[2 + 3 * j].im = cfrom1 * b_A[2 + 3 * j].im - (s.re * cfromc - s.im
            * ctoc);
          b_A[1 + 3 * j].re = atmp.re + cto1;
          b_A[1 + 3 * j].im = atmp.im + mul;
        }

        s.re = -s.re;
        s.im = -s.im;
        for (i = ilo - 1; i + 1 < 4; i++) {
          atmp.re = cfrom1 * b_A[6 + i].re;
          atmp.im = cfrom1 * b_A[6 + i].im;
          cto1 = s.re * b_A[3 + i].re - s.im * b_A[3 + i].im;
          mul = s.re * b_A[3 + i].im + s.im * b_A[3 + i].re;
          cfromc = b_A[6 + i].im;
          ctoc = b_A[6 + i].re;
          b_A[3 + i].re = cfrom1 * b_A[3 + i].re - (s.re * b_A[6 + i].re + s.im *
            b_A[6 + i].im);
          b_A[3 + i].im = cfrom1 * b_A[3 + i].im - (s.re * cfromc - s.im * ctoc);
          b_A[6 + i].re = atmp.re + cto1;
          b_A[6 + i].im = atmp.im + mul;
        }

        for (i = 0; i < 3; i++) {
          atmp.re = cfrom1 * V[6 + i].re;
          atmp.im = cfrom1 * V[6 + i].im;
          cto1 = s.re * V[3 + i].re - s.im * V[3 + i].im;
          mul = s.re * V[3 + i].im + s.im * V[3 + i].re;
          cfromc = V[6 + i].im;
          ctoc = V[6 + i].re;
          V[3 + i].re = cfrom1 * V[3 + i].re - (s.re * V[6 + i].re + s.im * V[6
            + i].im);
          V[3 + i].im = cfrom1 * V[3 + i].im - (s.re * cfromc - s.im * ctoc);
          V[6 + i].re = atmp.re + cto1;
          V[6 + i].im = atmp.im + mul;
        }

        ii = 2;
      }
    }

    b_eml_matlab_zhgeqz(b_A, ilo, ihi, V, info, alpha1, beta1);
    if (*info != 0.0) {
    } else {
      b_eml_matlab_ztgevc(b_A, V);
      if (ilo > 1) {
        for (i = ilo - 2; i + 1 >= 1; i--) {
          if (rscale[i] != i + 1) {
            for (j = 0; j < 3; j++) {
              atmp = V[i + 3 * j];
              V[i + 3 * j] = V[(rscale[i] + 3 * j) - 1];
              V[(rscale[i] + 3 * j) - 1] = atmp;
            }
          }
        }
      }

      if (ihi < 3) {
        while (ihi + 1 < 4) {
          if (rscale[ihi] != ihi + 1) {
            for (j = 0; j < 3; j++) {
              atmp = V[ihi + 3 * j];
              V[ihi + 3 * j] = V[(rscale[ihi] + 3 * j) - 1];
              V[(rscale[ihi] + 3 * j) - 1] = atmp;
            }
          }

          ihi++;
        }
      }

      for (ii = 0; ii < 3; ii++) {
        cfromc = fabs(V[3 * ii].re) + fabs(V[3 * ii].im);
        for (nzcount = 0; nzcount < 2; nzcount++) {
          ctoc = fabs(V[(nzcount + 3 * ii) + 1].re) + fabs(V[(nzcount + 3 * ii)
            + 1].im);
          if (ctoc > cfromc) {
            cfromc = ctoc;
          }
        }

        if (cfromc >= 6.7178761075670888E-139) {
          cfromc = 1.0 / cfromc;
          for (nzcount = 0; nzcount < 3; nzcount++) {
            V[nzcount + 3 * ii].re *= cfromc;
            V[nzcount + 3 * ii].im *= cfromc;
          }
        }
      }

      if (ilascl) {
        notdone = TRUE;
        while (notdone) {
          cfrom1 = anrmto * 2.0041683600089728E-292;
          cto1 = anrm / 4.9896007738368E+291;
          if ((fabs(cfrom1) > fabs(anrm)) && (anrm != 0.0)) {
            mul = 2.0041683600089728E-292;
            anrmto = cfrom1;
          } else if (fabs(cto1) > fabs(anrmto)) {
            mul = 4.9896007738368E+291;
            anrm = cto1;
          } else {
            mul = anrm / anrmto;
            notdone = FALSE;
          }

          for (ii = 0; ii < 3; ii++) {
            alpha1[ii].re *= mul;
            alpha1[ii].im *= mul;
          }
        }
      }
    }
  }
}

static void b_eml_matlab_zhgeqz(creal_T A[9], int32_T ilo, int32_T ihi, creal_T
  Z[9], real_T *info, creal_T alpha1[3], creal_T beta1[3])
{
  int32_T i;
  real_T eshift_re;
  real_T eshift_im;
  creal_T ctemp;
  real_T rho_re;
  real_T rho_im;
  real_T anorm;
  real_T temp;
  real_T b_atol;
  real_T ascale;
  boolean_T failed;
  int32_T j;
  boolean_T guard1 = FALSE;
  boolean_T guard2 = FALSE;
  int32_T ifirst;
  int32_T istart;
  int32_T ilast;
  int32_T ilastm1;
  int32_T iiter;
  boolean_T goto60;
  boolean_T goto70;
  boolean_T goto90;
  int32_T jiter;
  int32_T exitg1;
  boolean_T exitg3;
  boolean_T ilazro;
  int32_T b_j;
  boolean_T b_guard1 = FALSE;
  creal_T t1;
  creal_T d;
  creal_T sigma1;
  real_T sigma2_re;
  real_T sigma2_im;
  int32_T jp1;
  boolean_T exitg2;
  real_T tempr;
  creal_T dc2;
  for (i = 0; i < 3; i++) {
    alpha1[i].re = 0.0;
    alpha1[i].im = 0.0;
    beta1[i].re = 1.0;
    beta1[i].im = 0.0;
  }

  eshift_re = 0.0;
  eshift_im = 0.0;
  ctemp.re = 0.0;
  ctemp.im = 0.0;
  rho_re = 0.0;
  rho_im = 0.0;
  anorm = b_eml_matlab_zlanhs(A, ilo, ihi);
  temp = 2.2204460492503131E-16 * anorm;
  b_atol = 2.2250738585072014E-308;
  if (temp > 2.2250738585072014E-308) {
    b_atol = temp;
  }

  temp = 2.2250738585072014E-308;
  if (anorm > 2.2250738585072014E-308) {
    temp = anorm;
  }

  ascale = 1.0 / temp;
  failed = TRUE;
  for (j = ihi; j + 1 < 4; j++) {
    alpha1[j] = A[j + 3 * j];
  }

  guard1 = FALSE;
  guard2 = FALSE;
  if (ihi >= ilo) {
    ifirst = ilo;
    istart = ilo;
    ilast = ihi - 1;
    ilastm1 = ihi - 2;
    iiter = 0;
    goto60 = FALSE;
    goto70 = FALSE;
    goto90 = FALSE;
    jiter = 1;
    do {
      exitg1 = 0;
      if (jiter <= 30 * ((ihi - ilo) + 1)) {
        if (ilast + 1 == ilo) {
          goto60 = TRUE;
        } else if (fabs(A[ilast + 3 * ilastm1].re) + fabs(A[ilast + 3 * ilastm1]
                    .im) <= b_atol) {
          A[ilast + 3 * ilastm1].re = 0.0;
          A[ilast + 3 * ilastm1].im = 0.0;
          goto60 = TRUE;
        } else {
          j = ilastm1;
          exitg3 = FALSE;
          while ((exitg3 == 0U) && (j + 1 >= ilo)) {
            i = j - 1;
            if (j + 1 == ilo) {
              ilazro = TRUE;
            } else if (fabs(A[j + 3 * i].re) + fabs(A[j + 3 * i].im) <= b_atol)
            {
              A[j + 3 * i].re = 0.0;
              A[j + 3 * i].im = 0.0;
              ilazro = TRUE;
            } else {
              ilazro = FALSE;
            }

            if (ilazro) {
              ifirst = j + 1;
              goto70 = TRUE;
              exitg3 = TRUE;
            } else {
              j = i;
            }
          }
        }

        if (goto60 || goto70) {
          ilazro = TRUE;
        } else {
          ilazro = FALSE;
        }

        if (!ilazro) {
          for (i = 0; i < 3; i++) {
            alpha1[i].re = rtNaN;
            alpha1[i].im = 0.0;
            beta1[i].re = rtNaN;
            beta1[i].im = 0.0;
          }

          for (i = 0; i < 3; i++) {
            for (b_j = 0; b_j < 3; b_j++) {
              Z[b_j + 3 * i].re = rtNaN;
              Z[b_j + 3 * i].im = 0.0;
            }
          }

          *info = -1.0;
          exitg1 = 1;
        } else {
          b_guard1 = FALSE;
          if (goto60) {
            goto60 = FALSE;
            alpha1[ilast] = A[ilast + 3 * ilast];
            ilast = ilastm1;
            ilastm1--;
            if (ilast + 1 < ilo) {
              failed = FALSE;
              guard2 = TRUE;
              exitg1 = 1;
            } else {
              iiter = 0;
              eshift_re = 0.0;
              eshift_im = 0.0;
              b_guard1 = TRUE;
            }
          } else {
            if (goto70) {
              goto70 = FALSE;
              iiter++;
              if (iiter - div_s32_floor(iiter, 10) * 10 != 0) {
                t1.re = -(A[ilast + 3 * ilast].re - A[ilastm1 + 3 * ilastm1].re);
                t1.im = -(A[ilast + 3 * ilast].im - A[ilastm1 + 3 * ilastm1].im);
                t1 = eml_div(t1, 2.0);
                temp = A[ilastm1 + 3 * ilast].re * A[ilast + 3 * ilastm1].re -
                  A[ilastm1 + 3 * ilast].im * A[ilast + 3 * ilastm1].im;
                anorm = A[ilastm1 + 3 * ilast].re * A[ilast + 3 * ilastm1].im +
                  A[ilastm1 + 3 * ilast].im * A[ilast + 3 * ilastm1].re;
                d.re = (t1.re * t1.re - t1.im * t1.im) + temp;
                d.im = (t1.re * t1.im + t1.im * t1.re) + anorm;
                b_sqrt(&d);
                sigma1.re = A[ilastm1 + 3 * ilastm1].re - (t1.re - d.re);
                sigma1.im = A[ilastm1 + 3 * ilastm1].im - (t1.im - d.im);
                sigma2_re = A[ilastm1 + 3 * ilastm1].re - (t1.re + d.re);
                sigma2_im = A[ilastm1 + 3 * ilastm1].im - (t1.im + d.im);
                rho_re = sigma1.re - A[ilast + 3 * ilast].re;
                rho_im = sigma1.im - A[ilast + 3 * ilast].im;
                temp = sigma2_re - A[ilast + 3 * ilast].re;
                anorm = sigma2_im - A[ilast + 3 * ilast].im;
                if (rt_hypotd_snf(fabs(rho_re), fabs(rho_im)) <= rt_hypotd_snf
                    (fabs(temp), fabs(anorm))) {
                  sigma2_re = sigma1.re;
                  sigma2_im = sigma1.im;
                  rho_re = t1.re - d.re;
                  rho_im = t1.im - d.im;
                } else {
                  rho_re = t1.re + d.re;
                  rho_im = t1.im + d.im;
                }
              } else {
                eshift_re += A[ilast + 3 * ilastm1].re;
                eshift_im += A[ilast + 3 * ilastm1].im;
                sigma2_re = eshift_re;
                sigma2_im = eshift_im;
              }

              j = ilastm1;
              jp1 = ilastm1 + 1;
              exitg2 = FALSE;
              while ((exitg2 == 0U) && (j + 1 > ifirst)) {
                i = j - 1;
                istart = j + 1;
                ctemp.re = A[j + 3 * j].re - sigma2_re;
                ctemp.im = A[j + 3 * j].im - sigma2_im;
                temp = ascale * (fabs(ctemp.re) + fabs(ctemp.im));
                anorm = ascale * (fabs(A[jp1 + 3 * j].re) + fabs(A[jp1 + 3 * j].
                  im));
                tempr = temp;
                if (anorm > temp) {
                  tempr = anorm;
                }

                if ((tempr < 1.0) && (tempr != 0.0)) {
                  temp /= tempr;
                  anorm /= tempr;
                }

                if ((fabs(A[j + 3 * i].re) + fabs(A[j + 3 * i].im)) * anorm <=
                    temp * b_atol) {
                  goto90 = TRUE;
                  exitg2 = TRUE;
                } else {
                  jp1 = j;
                  j = i;
                }
              }

              if (!goto90) {
                istart = ifirst;
                if (ifirst == ilastm1 + 1) {
                  ctemp.re = rho_re;
                  ctemp.im = rho_im;
                } else {
                  ctemp.re = A[(ifirst + 3 * (ifirst - 1)) - 1].re - sigma2_re;
                  ctemp.im = A[(ifirst + 3 * (ifirst - 1)) - 1].im - sigma2_im;
                }

                goto90 = TRUE;
              }
            }

            if (goto90) {
              goto90 = FALSE;
              t1 = A[istart + 3 * (istart - 1)];
              b_eml_matlab_zlartg(ctemp, t1, &sigma2_im, &sigma1);
              j = istart - 1;
              i = istart - 2;
              while (j + 1 < ilast + 1) {
                jp1 = j + 1;
                if (j + 1 > istart) {
                  t1 = A[j + 3 * i];
                  d = A[jp1 + 3 * i];
                  eml_matlab_zlartg(t1, d, &sigma2_im, &sigma1, &dc2);
                  A[j + 3 * i] = dc2;
                  A[jp1 + 3 * i].re = 0.0;
                  A[jp1 + 3 * i].im = 0.0;
                }

                for (b_j = j; b_j + 1 < 4; b_j++) {
                  t1.re = sigma2_im * A[j + 3 * b_j].re;
                  t1.im = sigma2_im * A[j + 3 * b_j].im;
                  d.re = sigma1.re * A[jp1 + 3 * b_j].re - sigma1.im * A[jp1 + 3
                    * b_j].im;
                  d.im = sigma1.re * A[jp1 + 3 * b_j].im + sigma1.im * A[jp1 + 3
                    * b_j].re;
                  temp = A[j + 3 * b_j].re;
                  anorm = A[j + 3 * b_j].im;
                  tempr = A[j + 3 * b_j].im;
                  sigma2_re = A[j + 3 * b_j].re;
                  A[jp1 + 3 * b_j].re = sigma2_im * A[jp1 + 3 * b_j].re -
                    (sigma1.re * temp + sigma1.im * anorm);
                  A[jp1 + 3 * b_j].im = sigma2_im * A[jp1 + 3 * b_j].im -
                    (sigma1.re * tempr - sigma1.im * sigma2_re);
                  A[j + 3 * b_j].re = t1.re + d.re;
                  A[j + 3 * b_j].im = t1.im + d.im;
                }

                sigma1.re = -sigma1.re;
                sigma1.im = -sigma1.im;
                b_j = jp1 + 2;
                if (ilast + 1 < b_j) {
                  b_j = ilast + 1;
                }

                for (i = 0; i + 1 <= b_j; i++) {
                  t1.re = sigma2_im * A[i + 3 * jp1].re;
                  t1.im = sigma2_im * A[i + 3 * jp1].im;
                  d.re = sigma1.re * A[i + 3 * j].re - sigma1.im * A[i + 3 * j].
                    im;
                  d.im = sigma1.re * A[i + 3 * j].im + sigma1.im * A[i + 3 * j].
                    re;
                  temp = A[i + 3 * jp1].re;
                  anorm = A[i + 3 * jp1].im;
                  tempr = A[i + 3 * jp1].im;
                  sigma2_re = A[i + 3 * jp1].re;
                  A[i + 3 * j].re = sigma2_im * A[i + 3 * j].re - (sigma1.re *
                    temp + sigma1.im * anorm);
                  A[i + 3 * j].im = sigma2_im * A[i + 3 * j].im - (sigma1.re *
                    tempr - sigma1.im * sigma2_re);
                  A[i + 3 * jp1].re = t1.re + d.re;
                  A[i + 3 * jp1].im = t1.im + d.im;
                }

                for (i = 0; i < 3; i++) {
                  t1.re = sigma2_im * Z[i + 3 * jp1].re;
                  t1.im = sigma2_im * Z[i + 3 * jp1].im;
                  d.re = sigma1.re * Z[i + 3 * j].re - sigma1.im * Z[i + 3 * j].
                    im;
                  d.im = sigma1.re * Z[i + 3 * j].im + sigma1.im * Z[i + 3 * j].
                    re;
                  anorm = Z[i + 3 * jp1].re;
                  temp = Z[i + 3 * jp1].im;
                  tempr = Z[i + 3 * jp1].im;
                  sigma2_re = Z[i + 3 * jp1].re;
                  Z[i + 3 * j].re = sigma2_im * Z[i + 3 * j].re - (sigma1.re *
                    anorm + sigma1.im * temp);
                  Z[i + 3 * j].im = sigma2_im * Z[i + 3 * j].im - (sigma1.re *
                    tempr - sigma1.im * sigma2_re);
                  Z[i + 3 * jp1].re = t1.re + d.re;
                  Z[i + 3 * jp1].im = t1.im + d.im;
                }

                i = j;
                j = jp1;
              }
            }

            b_guard1 = TRUE;
          }

          if (b_guard1 == TRUE) {
            jiter++;
          }
        }
      } else {
        guard2 = TRUE;
        exitg1 = 1;
      }
    } while (exitg1 == 0U);
  } else {
    guard1 = TRUE;
  }

  if (guard2 == TRUE) {
    if (failed) {
      *info = (real_T)(ilast + 1);
      for (i = 0; i + 1 <= ilast + 1; i++) {
        alpha1[i].re = rtNaN;
        alpha1[i].im = 0.0;
        beta1[i].re = rtNaN;
        beta1[i].im = 0.0;
      }

      for (i = 0; i < 3; i++) {
        for (b_j = 0; b_j < 3; b_j++) {
          Z[b_j + 3 * i].re = rtNaN;
          Z[b_j + 3 * i].im = 0.0;
        }
      }
    } else {
      guard1 = TRUE;
    }
  }

  if (guard1 == TRUE) {
    for (j = 0; j + 1 <= ilo - 1; j++) {
      alpha1[j] = A[j + 3 * j];
    }

    *info = 0.0;
  }
}

static real_T b_eml_matlab_zlanhs(const creal_T A[9], int32_T ilo, int32_T ihi)
{
  real_T f;
  real_T scale;
  real_T sumsq;
  boolean_T firstNonZero;
  int32_T j;
  int32_T c;
  int32_T i;
  real_T temp1;
  real_T temp2;
  f = 0.0;
  if (ilo > ihi) {
  } else {
    scale = 0.0;
    sumsq = 0.0;
    firstNonZero = TRUE;
    for (j = ilo; j <= ihi; j++) {
      c = j + 1;
      if (ihi < c) {
        c = ihi;
      }

      for (i = ilo; i <= c; i++) {
        if (A[(i + 3 * (j - 1)) - 1].re != 0.0) {
          temp1 = fabs(A[(i + 3 * (j - 1)) - 1].re);
          if (firstNonZero) {
            sumsq = 1.0;
            scale = temp1;
            firstNonZero = FALSE;
          } else if (scale < temp1) {
            temp2 = scale / temp1;
            sumsq = 1.0 + sumsq * temp2 * temp2;
            scale = temp1;
          } else {
            temp2 = temp1 / scale;
            sumsq += temp2 * temp2;
          }
        }

        if (A[(i + 3 * (j - 1)) - 1].im != 0.0) {
          temp1 = fabs(A[(i + 3 * (j - 1)) - 1].im);
          if (firstNonZero) {
            sumsq = 1.0;
            scale = temp1;
            firstNonZero = FALSE;
          } else if (scale < temp1) {
            temp2 = scale / temp1;
            sumsq = 1.0 + sumsq * temp2 * temp2;
            scale = temp1;
          } else {
            temp2 = temp1 / scale;
            sumsq += temp2 * temp2;
          }
        }
      }
    }

    f = scale * sqrt(sumsq);
  }

  return f;
}

static void b_eml_matlab_zlartg(const creal_T f, const creal_T g, real_T *cs,
  creal_T *sn)
{
  real_T scale;
  real_T f2s;
  real_T g2;
  real_T fs_re;
  real_T fs_im;
  real_T gs_re;
  real_T gs_im;
  boolean_T guard1 = FALSE;
  real_T g2s;
  scale = fabs(f.re);
  f2s = fabs(f.im);
  if (f2s > scale) {
    scale = f2s;
  }

  f2s = fabs(g.re);
  g2 = fabs(g.im);
  if (g2 > f2s) {
    f2s = g2;
  }

  if (f2s > scale) {
    scale = f2s;
  }

  fs_re = f.re;
  fs_im = f.im;
  gs_re = g.re;
  gs_im = g.im;
  guard1 = FALSE;
  if (scale >= 7.4428285367870146E+137) {
    do {
      fs_re *= 1.3435752215134178E-138;
      fs_im *= 1.3435752215134178E-138;
      gs_re *= 1.3435752215134178E-138;
      gs_im *= 1.3435752215134178E-138;
      scale *= 1.3435752215134178E-138;
    } while (!(scale < 7.4428285367870146E+137));

    guard1 = TRUE;
  } else if (scale <= 1.3435752215134178E-138) {
    if ((g.re == 0.0) && (g.im == 0.0)) {
      *cs = 1.0;
      sn->re = 0.0;
      sn->im = 0.0;
    } else {
      do {
        fs_re *= 7.4428285367870146E+137;
        fs_im *= 7.4428285367870146E+137;
        gs_re *= 7.4428285367870146E+137;
        gs_im *= 7.4428285367870146E+137;
        scale *= 7.4428285367870146E+137;
      } while (!(scale > 1.3435752215134178E-138));

      guard1 = TRUE;
    }
  } else {
    guard1 = TRUE;
  }

  if (guard1 == TRUE) {
    scale = fs_re * fs_re + fs_im * fs_im;
    g2 = gs_re * gs_re + gs_im * gs_im;
    f2s = g2;
    if (1.0 > g2) {
      f2s = 1.0;
    }

    if (scale <= f2s * 2.0041683600089728E-292) {
      if ((f.re == 0.0) && (f.im == 0.0)) {
        *cs = 0.0;
        scale = rt_hypotd_snf(fabs(gs_re), fabs(gs_im));
        sn->re = gs_re / scale;
        sn->im = -gs_im / scale;
      } else {
        g2s = sqrt(g2);
        *cs = rt_hypotd_snf(fabs(fs_re), fabs(fs_im)) / g2s;
        f2s = fabs(f.re);
        g2 = fabs(f.im);
        if (g2 > f2s) {
          f2s = g2;
        }

        if (f2s > 1.0) {
          scale = rt_hypotd_snf(fabs(f.re), fabs(f.im));
          fs_re = f.re / scale;
          fs_im = f.im / scale;
        } else {
          f2s = 7.4428285367870146E+137 * f.re;
          g2 = 7.4428285367870146E+137 * f.im;
          scale = rt_hypotd_snf(fabs(f2s), fabs(g2));
          fs_re = f2s / scale;
          fs_im = g2 / scale;
        }

        gs_re /= g2s;
        gs_im = -gs_im / g2s;
        sn->re = fs_re * gs_re - fs_im * gs_im;
        sn->im = fs_re * gs_im + fs_im * gs_re;
      }
    } else {
      f2s = sqrt(1.0 + g2 / scale);
      *cs = 1.0 / f2s;
      scale += g2;
      fs_re = f2s * fs_re / scale;
      fs_im = f2s * fs_im / scale;
      sn->re = fs_re * gs_re - fs_im * -gs_im;
      sn->im = fs_re * -gs_im + fs_im * gs_re;
    }
  }
}

static void b_eml_matlab_ztgevc(const creal_T A[9], creal_T V[9])
{
  creal_T work1[3];
  creal_T work2[3];
  real_T rworka[3];
  int32_T i;
  real_T anorm;
  int32_T j;
  real_T y;
  real_T xmx;
  real_T ascale;
  int32_T je;
  real_T temp;
  real_T salpha_re;
  real_T salpha_im;
  real_T acoeff;
  boolean_T b2;
  boolean_T b3;
  real_T scale;
  real_T acoefa;
  int32_T jr;
  real_T dmin;
  creal_T d;
  creal_T b_work1;
  for (i = 0; i < 3; i++) {
    work1[i].re = 0.0;
    work1[i].im = 0.0;
    work2[i].re = 0.0;
    work2[i].im = 0.0;
    rworka[i] = 0.0;
  }

  anorm = fabs(A[0].re) + fabs(A[0].im);
  for (j = 0; j < 2; j++) {
    for (i = 0; i <= j; i++) {
      rworka[1 + j] += fabs(A[i + 3 * (1 + j)].re) + fabs(A[i + 3 * (1 + j)].im);
    }

    y = rworka[1 + j] + (fabs(A[(j + 3 * (1 + j)) + 1].re) + fabs(A[(j + 3 * (1
      + j)) + 1].im));
    if (y > anorm) {
      anorm = y;
    }
  }

  xmx = anorm;
  if (2.2250738585072014E-308 > anorm) {
    xmx = 2.2250738585072014E-308;
  }

  ascale = 1.0 / xmx;
  for (je = 0; je < 3; je++) {
    y = (fabs(A[(3 * (2 - je) - je) + 2].re) + fabs(A[(3 * (2 - je) - je) + 2].
          im)) * ascale;
    if (1.0 > y) {
      y = 1.0;
    }

    temp = 1.0 / y;
    salpha_re = ascale * (temp * A[(3 * (2 - je) - je) + 2].re);
    salpha_im = ascale * (temp * A[(3 * (2 - je) - je) + 2].im);
    acoeff = temp * ascale;
    if ((fabs(temp) >= 2.2250738585072014E-308) && (fabs(acoeff) <
         3.0062525400134592E-292)) {
      b2 = TRUE;
    } else {
      b2 = FALSE;
    }

    if ((fabs(salpha_re) + fabs(salpha_im) >= 2.2250738585072014E-308) && (fabs
         (salpha_re) + fabs(salpha_im) < 3.0062525400134592E-292)) {
      b3 = TRUE;
    } else {
      b3 = FALSE;
    }

    scale = 1.0;
    if (b2) {
      xmx = anorm;
      if (3.3264005158911995E+291 < anorm) {
        xmx = 3.3264005158911995E+291;
      }

      scale = 3.0062525400134592E-292 / fabs(temp) * xmx;
    }

    if (b3) {
      xmx = 3.0062525400134592E-292 / (fabs(salpha_re) + fabs(salpha_im));
      if (xmx > scale) {
        scale = xmx;
      }
    }

    if (b2 || b3) {
      y = fabs(acoeff);
      xmx = fabs(salpha_re) + fabs(salpha_im);
      if (1.0 > y) {
        y = 1.0;
      }

      if (xmx > y) {
        y = xmx;
      }

      y = 1.0 / (2.2250738585072014E-308 * y);
      if (y < scale) {
        scale = y;
      }

      if (b2) {
        acoeff = ascale * (scale * temp);
      } else {
        acoeff *= scale;
      }

      if (b3) {
        salpha_re *= scale;
        salpha_im *= scale;
      } else {
        salpha_re *= scale;
        salpha_im *= scale;
      }
    }

    acoefa = fabs(acoeff);
    for (jr = 0; jr < 3; jr++) {
      work1[jr].re = 0.0;
      work1[jr].im = 0.0;
    }

    work1[2 - je].re = 1.0;
    work1[2 - je].im = 0.0;
    dmin = 2.2204460492503131E-16 * acoefa * anorm;
    y = 2.2204460492503131E-16 * (fabs(salpha_re) + fabs(salpha_im));
    if (y > dmin) {
      dmin = y;
    }

    if (2.2250738585072014E-308 > dmin) {
      dmin = 2.2250738585072014E-308;
    }

    for (jr = 0; jr <= 1 - je; jr++) {
      work1[jr].re = acoeff * A[jr + 3 * (2 - je)].re;
      work1[jr].im = acoeff * A[jr + 3 * (2 - je)].im;
    }

    work1[2 - je].re = 1.0;
    work1[2 - je].im = 0.0;
    for (j = 0; j <= 1 - je; j++) {
      i = (-je - j) + 1;
      d.re = acoeff * A[i + 3 * i].re - salpha_re;
      d.im = acoeff * A[i + 3 * i].im - salpha_im;
      if (fabs(d.re) + fabs(d.im) <= dmin) {
        d.re = dmin;
        d.im = 0.0;
      }

      if ((fabs(d.re) + fabs(d.im) < 1.0) && (fabs(work1[i].re) + fabs(work1[i].
            im) >= 1.4980776123852632E+307 * (fabs(d.re) + fabs(d.im)))) {
        temp = 1.0 / (fabs(work1[i].re) + fabs(work1[i].im));
        for (jr = 0; jr <= 2 - je; jr++) {
          work1[jr].re *= temp;
          work1[jr].im *= temp;
        }
      }

      b_work1.re = -work1[i].re;
      b_work1.im = -work1[i].im;
      work1[i] = b_eml_div(b_work1, d);
      if (i + 1 > 1) {
        if (fabs(work1[1].re) + fabs(work1[1].im) > 1.0) {
          temp = 1.0 / (fabs(work1[1].re) + fabs(work1[1].im));
          if (acoefa * rworka[1] >= 1.4980776123852632E+307 * temp) {
            for (jr = 0; jr <= 2 - je; jr++) {
              work1[jr].re *= temp;
              work1[jr].im *= temp;
            }
          }
        }

        xmx = acoeff * work1[1].re;
        scale = acoeff * work1[1].im;
        work1[0].re += xmx * A[3].re - scale * A[3].im;
        work1[0].im += xmx * A[3].im + scale * A[3].re;
      }
    }

    for (jr = 0; jr < 3; jr++) {
      work2[jr].re = 0.0;
      work2[jr].im = 0.0;
    }

    for (i = 0; i <= 2 - je; i++) {
      for (jr = 0; jr < 3; jr++) {
        xmx = V[jr + 3 * i].re * work1[i].re - V[jr + 3 * i].im * work1[i].im;
        scale = V[jr + 3 * i].re * work1[i].im + V[jr + 3 * i].im * work1[i].re;
        work2[jr].re += xmx;
        work2[jr].im += scale;
      }
    }

    xmx = fabs(work2[0].re) + fabs(work2[0].im);
    for (jr = 0; jr < 2; jr++) {
      y = fabs(work2[1 + jr].re) + fabs(work2[1 + jr].im);
      if (y > xmx) {
        xmx = y;
      }
    }

    if (xmx > 2.2250738585072014E-308) {
      temp = 1.0 / xmx;
      for (jr = 0; jr < 3; jr++) {
        V[jr + 3 * (2 - je)].re = temp * work2[jr].re;
        V[jr + 3 * (2 - je)].im = temp * work2[jr].im;
      }
    } else {
      for (jr = 0; jr < 3; jr++) {
        V[jr + 3 * (2 - je)].re = 0.0;
        V[jr + 3 * (2 - je)].im = 0.0;
      }
    }
  }
}

static int32_T div_s32_floor(int32_T numerator, int32_T denominator)
{
  int32_T quotient;
  uint32_T absNumerator;
  uint32_T absDenominator;
  int32_T quotientNeedsNegation;
  uint32_T tempAbsQuotient;
  if (denominator == 0) {
    quotient = numerator >= 0 ? MAX_int32_T : MIN_int32_T;
  } else {
    absNumerator = (uint32_T)(numerator >= 0 ? numerator : -numerator);
    absDenominator = (uint32_T)(denominator >= 0 ? denominator : -denominator);
    quotientNeedsNegation = ((numerator < 0) != (denominator < 0));
    tempAbsQuotient = absNumerator / absDenominator;
    if ((uint32_T)quotientNeedsNegation) {
      absNumerator %= absDenominator;
      if (absNumerator > (uint32_T)0) {
        tempAbsQuotient++;
      }
    }

    quotient = (uint32_T)quotientNeedsNegation ? -(int32_T)tempAbsQuotient :
      (int32_T)tempAbsQuotient;
  }

  return quotient;
}

static void eml_matlab_zggbal(creal_T A[36], int32_T *ilo, int32_T *ihi, int32_T
  rscale[6])
{
  int32_T i;
  int32_T exitg2;
  int32_T j;
  boolean_T found;
  int32_T ii;
  boolean_T exitg5;
  int32_T nzcount;
  int32_T jj;
  boolean_T exitg6;
  boolean_T b_A;
  boolean_T guard2 = FALSE;
  real_T atmp_re;
  real_T atmp_im;
  int32_T exitg1;
  boolean_T exitg3;
  boolean_T exitg4;
  boolean_T guard1 = FALSE;
  for (i = 0; i < 6; i++) {
    rscale[i] = 0;
  }

  *ilo = 1;
  *ihi = 6;
  do {
    exitg2 = 0;
    i = 0;
    j = 0;
    found = FALSE;
    ii = *ihi;
    exitg5 = FALSE;
    while ((exitg5 == 0U) && (ii > 0)) {
      nzcount = 0;
      i = ii;
      j = *ihi;
      jj = 1;
      exitg6 = FALSE;
      while ((exitg6 == 0U) && (jj <= *ihi)) {
        b_A = ((A[(ii + 6 * (jj - 1)) - 1].re != 0.0) || (A[(ii + 6 * (jj - 1))
                - 1].im != 0.0));
        guard2 = FALSE;
        if (b_A || (ii == jj)) {
          if (nzcount == 0) {
            j = jj;
            nzcount = 1;
            guard2 = TRUE;
          } else {
            nzcount = 2;
            exitg6 = TRUE;
          }
        } else {
          guard2 = TRUE;
        }

        if (guard2 == TRUE) {
          jj++;
        }
      }

      if (nzcount < 2) {
        found = TRUE;
        exitg5 = TRUE;
      } else {
        ii--;
      }
    }

    if (!found) {
      exitg2 = 2;
    } else {
      if (i != *ihi) {
        for (ii = 0; ii < 6; ii++) {
          atmp_re = A[(i + 6 * ii) - 1].re;
          atmp_im = A[(i + 6 * ii) - 1].im;
          A[(i + 6 * ii) - 1] = A[(*ihi + 6 * ii) - 1];
          A[(*ihi + 6 * ii) - 1].re = atmp_re;
          A[(*ihi + 6 * ii) - 1].im = atmp_im;
        }
      }

      if (j != *ihi) {
        for (ii = 0; ii + 1 <= *ihi; ii++) {
          atmp_re = A[ii + 6 * (j - 1)].re;
          atmp_im = A[ii + 6 * (j - 1)].im;
          A[ii + 6 * (j - 1)] = A[ii + 6 * (*ihi - 1)];
          A[ii + 6 * (*ihi - 1)].re = atmp_re;
          A[ii + 6 * (*ihi - 1)].im = atmp_im;
        }
      }

      rscale[*ihi - 1] = j;
      (*ihi)--;
      if (*ihi == 1) {
        rscale[0] = 1;
        exitg2 = 1;
      }
    }
  } while (exitg2 == 0U);

  if (exitg2 == 1U) {
  } else {
    do {
      exitg1 = 0;
      i = 0;
      j = 0;
      found = FALSE;
      jj = *ilo;
      exitg3 = FALSE;
      while ((exitg3 == 0U) && (jj <= *ihi)) {
        nzcount = 0;
        i = *ihi;
        j = jj;
        ii = *ilo;
        exitg4 = FALSE;
        while ((exitg4 == 0U) && (ii <= *ihi)) {
          b_A = ((A[(ii + 6 * (jj - 1)) - 1].re != 0.0) || (A[(ii + 6 * (jj - 1))
                  - 1].im != 0.0));
          guard1 = FALSE;
          if (b_A || (ii == jj)) {
            if (nzcount == 0) {
              i = ii;
              nzcount = 1;
              guard1 = TRUE;
            } else {
              nzcount = 2;
              exitg4 = TRUE;
            }
          } else {
            guard1 = TRUE;
          }

          if (guard1 == TRUE) {
            ii++;
          }
        }

        if (nzcount < 2) {
          found = TRUE;
          exitg3 = TRUE;
        } else {
          jj++;
        }
      }

      if (!found) {
        exitg1 = 1;
      } else {
        if (i != *ilo) {
          for (ii = *ilo - 1; ii + 1 < 7; ii++) {
            atmp_re = A[(i + 6 * ii) - 1].re;
            atmp_im = A[(i + 6 * ii) - 1].im;
            A[(i + 6 * ii) - 1] = A[(*ilo + 6 * ii) - 1];
            A[(*ilo + 6 * ii) - 1].re = atmp_re;
            A[(*ilo + 6 * ii) - 1].im = atmp_im;
          }
        }

        if (j != *ilo) {
          for (ii = 0; ii + 1 <= *ihi; ii++) {
            atmp_re = A[ii + 6 * (j - 1)].re;
            atmp_im = A[ii + 6 * (j - 1)].im;
            A[ii + 6 * (j - 1)] = A[ii + 6 * (*ilo - 1)];
            A[ii + 6 * (*ilo - 1)].re = atmp_re;
            A[ii + 6 * (*ilo - 1)].im = atmp_im;
          }
        }

        rscale[*ilo - 1] = j;
        (*ilo)++;
        if (*ilo == *ihi) {
          rscale[*ilo - 1] = *ilo;
          exitg1 = 1;
        }
      }
    } while (exitg1 == 0U);
  }
}

static void eml_matlab_zggev(creal_T A[36], real_T *info, creal_T alpha1[6],
  creal_T beta1[6], creal_T V[36])
{
  real_T anrm;
  int32_T jcol;
  boolean_T exitg1;
  real_T absxk;
  int32_T i;
  boolean_T ilascl;
  real_T anrmto;
  real_T ctoc;
  boolean_T notdone;
  real_T cfrom1;
  real_T cto1;
  real_T mul;
  creal_T b_A[36];
  int32_T rscale[6];
  int32_T ihi;
  int32_T ilo;
  static const int8_T iv0[36] = { 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0,
    0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1 };

  int32_T jcolp1;
  int32_T jrow;
  int32_T jrowm1;
  creal_T tmp;
  creal_T s;
  int32_T j;
  *info = 0.0;
  anrm = 0.0;
  jcol = 0;
  exitg1 = FALSE;
  while ((exitg1 == 0U) && (jcol < 36)) {
    absxk = rt_hypotd_snf(fabs(A[jcol].re), fabs(A[jcol].im));
    if (rtIsNaN(absxk)) {
      anrm = rtNaN;
      exitg1 = TRUE;
    } else {
      if (absxk > anrm) {
        anrm = absxk;
      }

      jcol++;
    }
  }

  if (!((!rtIsInf(anrm)) && (!rtIsNaN(anrm)))) {
    for (i = 0; i < 6; i++) {
      alpha1[i].re = rtNaN;
      alpha1[i].im = 0.0;
      beta1[i].re = rtNaN;
      beta1[i].im = 0.0;
    }

    for (jcol = 0; jcol < 36; jcol++) {
      V[jcol].re = rtNaN;
      V[jcol].im = 0.0;
    }
  } else {
    ilascl = FALSE;
    anrmto = anrm;
    if ((anrm > 0.0) && (anrm < 6.7178761075670888E-139)) {
      anrmto = 6.7178761075670888E-139;
      ilascl = TRUE;
    } else {
      if (anrm > 1.4885657073574029E+138) {
        anrmto = 1.4885657073574029E+138;
        ilascl = TRUE;
      }
    }

    if (ilascl) {
      absxk = anrm;
      ctoc = anrmto;
      notdone = TRUE;
      while (notdone) {
        cfrom1 = absxk * 2.0041683600089728E-292;
        cto1 = ctoc / 4.9896007738368E+291;
        if ((cfrom1 > ctoc) && (ctoc != 0.0)) {
          mul = 2.0041683600089728E-292;
          absxk = cfrom1;
        } else if (cto1 > absxk) {
          mul = 4.9896007738368E+291;
          ctoc = cto1;
        } else {
          mul = ctoc / absxk;
          notdone = FALSE;
        }

        for (jcol = 0; jcol < 36; jcol++) {
          A[jcol].re *= mul;
          A[jcol].im *= mul;
        }
      }
    }

    memcpy(&b_A[0], &A[0], 36U * sizeof(creal_T));
    eml_matlab_zggbal(b_A, &ilo, &ihi, rscale);
    for (jcol = 0; jcol < 36; jcol++) {
      V[jcol].re = (real_T)iv0[jcol];
      V[jcol].im = 0.0;
    }

    if (ihi < ilo + 2) {
    } else {
      jcol = ilo - 1;
      while (jcol + 1 < ihi - 1) {
        jcolp1 = jcol + 1;
        jrow = ihi - 1;
        while (jrow + 1 > jcolp1 + 1) {
          jrowm1 = jrow - 1;
          eml_matlab_zlartg(b_A[jrowm1 + 6 * jcol], b_A[jrow + 6 * jcol],
                            &cfrom1, &s, &tmp);
          b_A[jrowm1 + 6 * jcol] = tmp;
          b_A[jrow + 6 * jcol].re = 0.0;
          b_A[jrow + 6 * jcol].im = 0.0;
          for (j = jcolp1; j + 1 <= ihi; j++) {
            tmp.re = cfrom1 * b_A[jrowm1 + 6 * j].re;
            tmp.im = cfrom1 * b_A[jrowm1 + 6 * j].im;
            cto1 = s.re * b_A[jrow + 6 * j].re - s.im * b_A[jrow + 6 * j].im;
            mul = s.re * b_A[jrow + 6 * j].im + s.im * b_A[jrow + 6 * j].re;
            absxk = b_A[jrowm1 + 6 * j].im;
            ctoc = b_A[jrowm1 + 6 * j].re;
            b_A[jrow + 6 * j].re = cfrom1 * b_A[jrow + 6 * j].re - (s.re *
              b_A[jrowm1 + 6 * j].re + s.im * b_A[jrowm1 + 6 * j].im);
            b_A[jrow + 6 * j].im = cfrom1 * b_A[jrow + 6 * j].im - (s.re * absxk
              - s.im * ctoc);
            b_A[jrowm1 + 6 * j].re = tmp.re + cto1;
            b_A[jrowm1 + 6 * j].im = tmp.im + mul;
          }

          s.re = -s.re;
          s.im = -s.im;
          for (i = ilo - 1; i + 1 <= ihi; i++) {
            tmp.re = cfrom1 * b_A[i + 6 * jrow].re;
            tmp.im = cfrom1 * b_A[i + 6 * jrow].im;
            cto1 = s.re * b_A[i + 6 * jrowm1].re - s.im * b_A[i + 6 * jrowm1].im;
            mul = s.re * b_A[i + 6 * jrowm1].im + s.im * b_A[i + 6 * jrowm1].re;
            absxk = b_A[i + 6 * jrow].im;
            ctoc = b_A[i + 6 * jrow].re;
            b_A[i + 6 * jrowm1].re = cfrom1 * b_A[i + 6 * jrowm1].re - (s.re *
              b_A[i + 6 * jrow].re + s.im * b_A[i + 6 * jrow].im);
            b_A[i + 6 * jrowm1].im = cfrom1 * b_A[i + 6 * jrowm1].im - (s.re *
              absxk - s.im * ctoc);
            b_A[i + 6 * jrow].re = tmp.re + cto1;
            b_A[i + 6 * jrow].im = tmp.im + mul;
          }

          for (i = 0; i < 6; i++) {
            tmp.re = cfrom1 * V[i + 6 * jrow].re;
            tmp.im = cfrom1 * V[i + 6 * jrow].im;
            cto1 = s.re * V[i + 6 * jrowm1].re - s.im * V[i + 6 * jrowm1].im;
            mul = s.re * V[i + 6 * jrowm1].im + s.im * V[i + 6 * jrowm1].re;
            absxk = V[i + 6 * jrow].im;
            ctoc = V[i + 6 * jrow].re;
            V[i + 6 * jrowm1].re = cfrom1 * V[i + 6 * jrowm1].re - (s.re * V[i +
              6 * jrow].re + s.im * V[i + 6 * jrow].im);
            V[i + 6 * jrowm1].im = cfrom1 * V[i + 6 * jrowm1].im - (s.re * absxk
              - s.im * ctoc);
            V[i + 6 * jrow].re = tmp.re + cto1;
            V[i + 6 * jrow].im = tmp.im + mul;
          }

          jrow = jrowm1;
        }

        jcol = jcolp1;
      }
    }

    eml_matlab_zhgeqz(b_A, ilo, ihi, V, info, alpha1, beta1);
    if (*info != 0.0) {
    } else {
      eml_matlab_ztgevc(b_A, V);
      if (ilo > 1) {
        for (i = ilo - 2; i + 1 >= 1; i--) {
          if (rscale[i] != i + 1) {
            for (j = 0; j < 6; j++) {
              tmp = V[i + 6 * j];
              V[i + 6 * j] = V[(rscale[i] + 6 * j) - 1];
              V[(rscale[i] + 6 * j) - 1] = tmp;
            }
          }
        }
      }

      if (ihi < 6) {
        while (ihi + 1 < 7) {
          if (rscale[ihi] != ihi + 1) {
            for (j = 0; j < 6; j++) {
              tmp = V[ihi + 6 * j];
              V[ihi + 6 * j] = V[(rscale[ihi] + 6 * j) - 1];
              V[(rscale[ihi] + 6 * j) - 1] = tmp;
            }
          }

          ihi++;
        }
      }

      for (jcol = 0; jcol < 6; jcol++) {
        absxk = fabs(V[6 * jcol].re) + fabs(V[6 * jcol].im);
        for (jcolp1 = 0; jcolp1 < 5; jcolp1++) {
          ctoc = fabs(V[(jcolp1 + 6 * jcol) + 1].re) + fabs(V[(jcolp1 + 6 * jcol)
            + 1].im);
          if (ctoc > absxk) {
            absxk = ctoc;
          }
        }

        if (absxk >= 6.7178761075670888E-139) {
          absxk = 1.0 / absxk;
          for (jcolp1 = 0; jcolp1 < 6; jcolp1++) {
            V[jcolp1 + 6 * jcol].re *= absxk;
            V[jcolp1 + 6 * jcol].im *= absxk;
          }
        }
      }

      if (ilascl) {
        notdone = TRUE;
        while (notdone) {
          cfrom1 = anrmto * 2.0041683600089728E-292;
          cto1 = anrm / 4.9896007738368E+291;
          if ((cfrom1 > anrm) && (anrm != 0.0)) {
            mul = 2.0041683600089728E-292;
            anrmto = cfrom1;
          } else if (cto1 > anrmto) {
            mul = 4.9896007738368E+291;
            anrm = cto1;
          } else {
            mul = anrm / anrmto;
            notdone = FALSE;
          }

          for (jcol = 0; jcol < 6; jcol++) {
            alpha1[jcol].re *= mul;
            alpha1[jcol].im *= mul;
          }
        }
      }
    }
  }
}

static void eml_matlab_zhgeqz(creal_T A[36], int32_T ilo, int32_T ihi, creal_T
  Z[36], real_T *info, creal_T alpha1[6], creal_T beta1[6])
{
  int32_T i;
  real_T eshift_re;
  real_T eshift_im;
  creal_T ctemp;
  real_T rho_re;
  real_T rho_im;
  real_T anorm;
  real_T temp;
  real_T b_atol;
  real_T ascale;
  boolean_T failed;
  int32_T j;
  boolean_T guard1 = FALSE;
  boolean_T guard2 = FALSE;
  int32_T ifirst;
  int32_T istart;
  int32_T ilast;
  int32_T ilastm1;
  int32_T iiter;
  boolean_T goto60;
  boolean_T goto70;
  boolean_T goto90;
  int32_T jiter;
  int32_T exitg1;
  boolean_T exitg3;
  boolean_T ilazro;
  int32_T b_j;
  boolean_T b_guard1 = FALSE;
  creal_T t1;
  creal_T d;
  creal_T sigma1;
  real_T sigma2_re;
  real_T sigma2_im;
  int32_T jp1;
  boolean_T exitg2;
  real_T tempr;
  creal_T dc0;
  for (i = 0; i < 6; i++) {
    alpha1[i].re = 0.0;
    alpha1[i].im = 0.0;
    beta1[i].re = 1.0;
    beta1[i].im = 0.0;
  }

  eshift_re = 0.0;
  eshift_im = 0.0;
  ctemp.re = 0.0;
  ctemp.im = 0.0;
  rho_re = 0.0;
  rho_im = 0.0;
  anorm = eml_matlab_zlanhs(A, ilo, ihi);
  temp = 2.2204460492503131E-16 * anorm;
  b_atol = 2.2250738585072014E-308;
  if (temp > 2.2250738585072014E-308) {
    b_atol = temp;
  }

  temp = 2.2250738585072014E-308;
  if (anorm > 2.2250738585072014E-308) {
    temp = anorm;
  }

  ascale = 1.0 / temp;
  failed = TRUE;
  for (j = ihi; j + 1 < 7; j++) {
    alpha1[j] = A[j + 6 * j];
  }

  guard1 = FALSE;
  guard2 = FALSE;
  if (ihi >= ilo) {
    ifirst = ilo;
    istart = ilo;
    ilast = ihi - 1;
    ilastm1 = ihi - 2;
    iiter = 0;
    goto60 = FALSE;
    goto70 = FALSE;
    goto90 = FALSE;
    jiter = 1;
    do {
      exitg1 = 0;
      if (jiter <= 30 * ((ihi - ilo) + 1)) {
        if (ilast + 1 == ilo) {
          goto60 = TRUE;
        } else if (fabs(A[ilast + 6 * ilastm1].re) + fabs(A[ilast + 6 * ilastm1]
                    .im) <= b_atol) {
          A[ilast + 6 * ilastm1].re = 0.0;
          A[ilast + 6 * ilastm1].im = 0.0;
          goto60 = TRUE;
        } else {
          j = ilastm1;
          exitg3 = FALSE;
          while ((exitg3 == 0U) && (j + 1 >= ilo)) {
            i = j - 1;
            if (j + 1 == ilo) {
              ilazro = TRUE;
            } else if (fabs(A[j + 6 * i].re) + fabs(A[j + 6 * i].im) <= b_atol)
            {
              A[j + 6 * i].re = 0.0;
              A[j + 6 * i].im = 0.0;
              ilazro = TRUE;
            } else {
              ilazro = FALSE;
            }

            if (ilazro) {
              ifirst = j + 1;
              goto70 = TRUE;
              exitg3 = TRUE;
            } else {
              j = i;
            }
          }
        }

        if (goto60 || goto70) {
          ilazro = TRUE;
        } else {
          ilazro = FALSE;
        }

        if (!ilazro) {
          for (i = 0; i < 6; i++) {
            alpha1[i].re = rtNaN;
            alpha1[i].im = 0.0;
            beta1[i].re = rtNaN;
            beta1[i].im = 0.0;
          }

          for (i = 0; i < 6; i++) {
            for (b_j = 0; b_j < 6; b_j++) {
              Z[b_j + 6 * i].re = rtNaN;
              Z[b_j + 6 * i].im = 0.0;
            }
          }

          *info = -1.0;
          exitg1 = 1;
        } else {
          b_guard1 = FALSE;
          if (goto60) {
            goto60 = FALSE;
            alpha1[ilast] = A[ilast + 6 * ilast];
            ilast = ilastm1;
            ilastm1--;
            if (ilast + 1 < ilo) {
              failed = FALSE;
              guard2 = TRUE;
              exitg1 = 1;
            } else {
              iiter = 0;
              eshift_re = 0.0;
              eshift_im = 0.0;
              b_guard1 = TRUE;
            }
          } else {
            if (goto70) {
              goto70 = FALSE;
              iiter++;
              if (iiter - div_s32_floor(iiter, 10) * 10 != 0) {
                t1.re = -(A[ilast + 6 * ilast].re - A[ilastm1 + 6 * ilastm1].re);
                t1.im = -(A[ilast + 6 * ilast].im - A[ilastm1 + 6 * ilastm1].im);
                t1 = eml_div(t1, 2.0);
                temp = A[ilastm1 + 6 * ilast].re * A[ilast + 6 * ilastm1].re -
                  A[ilastm1 + 6 * ilast].im * A[ilast + 6 * ilastm1].im;
                anorm = A[ilastm1 + 6 * ilast].re * A[ilast + 6 * ilastm1].im +
                  A[ilastm1 + 6 * ilast].im * A[ilast + 6 * ilastm1].re;
                d.re = (t1.re * t1.re - t1.im * t1.im) + temp;
                d.im = (t1.re * t1.im + t1.im * t1.re) + anorm;
                b_sqrt(&d);
                sigma1.re = A[ilastm1 + 6 * ilastm1].re - (t1.re - d.re);
                sigma1.im = A[ilastm1 + 6 * ilastm1].im - (t1.im - d.im);
                sigma2_re = A[ilastm1 + 6 * ilastm1].re - (t1.re + d.re);
                sigma2_im = A[ilastm1 + 6 * ilastm1].im - (t1.im + d.im);
                rho_re = sigma1.re - A[ilast + 6 * ilast].re;
                rho_im = sigma1.im - A[ilast + 6 * ilast].im;
                temp = sigma2_re - A[ilast + 6 * ilast].re;
                anorm = sigma2_im - A[ilast + 6 * ilast].im;
                if (rt_hypotd_snf(fabs(rho_re), fabs(rho_im)) <= rt_hypotd_snf
                    (fabs(temp), fabs(anorm))) {
                  sigma2_re = sigma1.re;
                  sigma2_im = sigma1.im;
                  rho_re = t1.re - d.re;
                  rho_im = t1.im - d.im;
                } else {
                  rho_re = t1.re + d.re;
                  rho_im = t1.im + d.im;
                }
              } else {
                eshift_re += A[ilast + 6 * ilastm1].re;
                eshift_im += A[ilast + 6 * ilastm1].im;
                sigma2_re = eshift_re;
                sigma2_im = eshift_im;
              }

              j = ilastm1;
              jp1 = ilastm1 + 1;
              exitg2 = FALSE;
              while ((exitg2 == 0U) && (j + 1 > ifirst)) {
                i = j - 1;
                istart = j + 1;
                ctemp.re = A[j + 6 * j].re - sigma2_re;
                ctemp.im = A[j + 6 * j].im - sigma2_im;
                temp = ascale * (fabs(ctemp.re) + fabs(ctemp.im));
                anorm = ascale * (fabs(A[jp1 + 6 * j].re) + fabs(A[jp1 + 6 * j].
                  im));
                tempr = temp;
                if (anorm > temp) {
                  tempr = anorm;
                }

                if ((tempr < 1.0) && (tempr != 0.0)) {
                  temp /= tempr;
                  anorm /= tempr;
                }

                if ((fabs(A[j + 6 * i].re) + fabs(A[j + 6 * i].im)) * anorm <=
                    temp * b_atol) {
                  goto90 = TRUE;
                  exitg2 = TRUE;
                } else {
                  jp1 = j;
                  j = i;
                }
              }

              if (!goto90) {
                istart = ifirst;
                if (ifirst == ilastm1 + 1) {
                  ctemp.re = rho_re;
                  ctemp.im = rho_im;
                } else {
                  ctemp.re = A[(ifirst + 6 * (ifirst - 1)) - 1].re - sigma2_re;
                  ctemp.im = A[(ifirst + 6 * (ifirst - 1)) - 1].im - sigma2_im;
                }

                goto90 = TRUE;
              }
            }

            if (goto90) {
              goto90 = FALSE;
              t1 = A[istart + 6 * (istart - 1)];
              b_eml_matlab_zlartg(ctemp, t1, &sigma2_im, &sigma1);
              j = istart - 1;
              i = istart - 2;
              while (j + 1 < ilast + 1) {
                jp1 = j + 1;
                if (j + 1 > istart) {
                  t1 = A[j + 6 * i];
                  d = A[jp1 + 6 * i];
                  eml_matlab_zlartg(t1, d, &sigma2_im, &sigma1, &dc0);
                  A[j + 6 * i] = dc0;
                  A[jp1 + 6 * i].re = 0.0;
                  A[jp1 + 6 * i].im = 0.0;
                }

                for (b_j = j; b_j + 1 < 7; b_j++) {
                  t1.re = sigma2_im * A[j + 6 * b_j].re;
                  t1.im = sigma2_im * A[j + 6 * b_j].im;
                  d.re = sigma1.re * A[jp1 + 6 * b_j].re - sigma1.im * A[jp1 + 6
                    * b_j].im;
                  d.im = sigma1.re * A[jp1 + 6 * b_j].im + sigma1.im * A[jp1 + 6
                    * b_j].re;
                  temp = A[j + 6 * b_j].re;
                  anorm = A[j + 6 * b_j].im;
                  tempr = A[j + 6 * b_j].im;
                  sigma2_re = A[j + 6 * b_j].re;
                  A[jp1 + 6 * b_j].re = sigma2_im * A[jp1 + 6 * b_j].re -
                    (sigma1.re * temp + sigma1.im * anorm);
                  A[jp1 + 6 * b_j].im = sigma2_im * A[jp1 + 6 * b_j].im -
                    (sigma1.re * tempr - sigma1.im * sigma2_re);
                  A[j + 6 * b_j].re = t1.re + d.re;
                  A[j + 6 * b_j].im = t1.im + d.im;
                }

                sigma1.re = -sigma1.re;
                sigma1.im = -sigma1.im;
                b_j = jp1 + 2;
                if (ilast + 1 < b_j) {
                  b_j = ilast + 1;
                }

                for (i = 0; i + 1 <= b_j; i++) {
                  t1.re = sigma2_im * A[i + 6 * jp1].re;
                  t1.im = sigma2_im * A[i + 6 * jp1].im;
                  d.re = sigma1.re * A[i + 6 * j].re - sigma1.im * A[i + 6 * j].
                    im;
                  d.im = sigma1.re * A[i + 6 * j].im + sigma1.im * A[i + 6 * j].
                    re;
                  temp = A[i + 6 * jp1].re;
                  anorm = A[i + 6 * jp1].im;
                  tempr = A[i + 6 * jp1].im;
                  sigma2_re = A[i + 6 * jp1].re;
                  A[i + 6 * j].re = sigma2_im * A[i + 6 * j].re - (sigma1.re *
                    temp + sigma1.im * anorm);
                  A[i + 6 * j].im = sigma2_im * A[i + 6 * j].im - (sigma1.re *
                    tempr - sigma1.im * sigma2_re);
                  A[i + 6 * jp1].re = t1.re + d.re;
                  A[i + 6 * jp1].im = t1.im + d.im;
                }

                for (i = 0; i < 6; i++) {
                  t1.re = sigma2_im * Z[i + 6 * jp1].re;
                  t1.im = sigma2_im * Z[i + 6 * jp1].im;
                  d.re = sigma1.re * Z[i + 6 * j].re - sigma1.im * Z[i + 6 * j].
                    im;
                  d.im = sigma1.re * Z[i + 6 * j].im + sigma1.im * Z[i + 6 * j].
                    re;
                  anorm = Z[i + 6 * jp1].re;
                  temp = Z[i + 6 * jp1].im;
                  tempr = Z[i + 6 * jp1].im;
                  sigma2_re = Z[i + 6 * jp1].re;
                  Z[i + 6 * j].re = sigma2_im * Z[i + 6 * j].re - (sigma1.re *
                    anorm + sigma1.im * temp);
                  Z[i + 6 * j].im = sigma2_im * Z[i + 6 * j].im - (sigma1.re *
                    tempr - sigma1.im * sigma2_re);
                  Z[i + 6 * jp1].re = t1.re + d.re;
                  Z[i + 6 * jp1].im = t1.im + d.im;
                }

                i = j;
                j = jp1;
              }
            }

            b_guard1 = TRUE;
          }

          if (b_guard1 == TRUE) {
            jiter++;
          }
        }
      } else {
        guard2 = TRUE;
        exitg1 = 1;
      }
    } while (exitg1 == 0U);
  } else {
    guard1 = TRUE;
  }

  if (guard2 == TRUE) {
    if (failed) {
      *info = (real_T)(ilast + 1);
      for (i = 0; i + 1 <= ilast + 1; i++) {
        alpha1[i].re = rtNaN;
        alpha1[i].im = 0.0;
        beta1[i].re = rtNaN;
        beta1[i].im = 0.0;
      }

      for (i = 0; i < 6; i++) {
        for (b_j = 0; b_j < 6; b_j++) {
          Z[b_j + 6 * i].re = rtNaN;
          Z[b_j + 6 * i].im = 0.0;
        }
      }
    } else {
      guard1 = TRUE;
    }
  }

  if (guard1 == TRUE) {
    for (j = 0; j + 1 <= ilo - 1; j++) {
      alpha1[j] = A[j + 6 * j];
    }

    *info = 0.0;
  }
}

static real_T eml_matlab_zlangeM(const creal_T x[9])
{
  real_T y;
  int32_T k;
  boolean_T exitg1;
  real_T absxk;
  y = 0.0;
  k = 0;
  exitg1 = FALSE;
  while ((exitg1 == 0U) && (k < 9)) {
    absxk = rt_hypotd_snf(fabs(x[k].re), fabs(x[k].im));
    if (rtIsNaN(absxk)) {
      y = rtNaN;
      exitg1 = TRUE;
    } else {
      if (absxk > y) {
        y = absxk;
      }

      k++;
    }
  }

  return y;
}

static real_T eml_matlab_zlanhs(const creal_T A[36], int32_T ilo, int32_T ihi)
{
  real_T f;
  real_T scale;
  real_T sumsq;
  boolean_T firstNonZero;
  int32_T j;
  int32_T c;
  int32_T i;
  real_T temp1;
  real_T temp2;
  f = 0.0;
  if (ilo > ihi) {
  } else {
    scale = 0.0;
    sumsq = 0.0;
    firstNonZero = TRUE;
    for (j = ilo; j <= ihi; j++) {
      c = j + 1;
      if (ihi < c) {
        c = ihi;
      }

      for (i = ilo; i <= c; i++) {
        if (A[(i + 6 * (j - 1)) - 1].re != 0.0) {
          temp1 = fabs(A[(i + 6 * (j - 1)) - 1].re);
          if (firstNonZero) {
            sumsq = 1.0;
            scale = temp1;
            firstNonZero = FALSE;
          } else if (scale < temp1) {
            temp2 = scale / temp1;
            sumsq = 1.0 + sumsq * temp2 * temp2;
            scale = temp1;
          } else {
            temp2 = temp1 / scale;
            sumsq += temp2 * temp2;
          }
        }

        if (A[(i + 6 * (j - 1)) - 1].im != 0.0) {
          temp1 = fabs(A[(i + 6 * (j - 1)) - 1].im);
          if (firstNonZero) {
            sumsq = 1.0;
            scale = temp1;
            firstNonZero = FALSE;
          } else if (scale < temp1) {
            temp2 = scale / temp1;
            sumsq = 1.0 + sumsq * temp2 * temp2;
            scale = temp1;
          } else {
            temp2 = temp1 / scale;
            sumsq += temp2 * temp2;
          }
        }
      }
    }

    f = scale * sqrt(sumsq);
  }

  return f;
}

static void eml_matlab_zlartg(const creal_T f, const creal_T g, real_T *cs,
  creal_T *sn, creal_T *r)
{
  real_T scale;
  real_T f2s;
  real_T g2;
  real_T fs_re;
  real_T fs_im;
  real_T gs_re;
  real_T gs_im;
  int32_T count;
  int32_T rescaledir;
  boolean_T guard1 = FALSE;
  real_T g2s;
  scale = fabs(f.re);
  f2s = fabs(f.im);
  if (f2s > scale) {
    scale = f2s;
  }

  f2s = fabs(g.re);
  g2 = fabs(g.im);
  if (g2 > f2s) {
    f2s = g2;
  }

  if (f2s > scale) {
    scale = f2s;
  }

  fs_re = f.re;
  fs_im = f.im;
  gs_re = g.re;
  gs_im = g.im;
  count = 0;
  rescaledir = 0;
  guard1 = FALSE;
  if (scale >= 7.4428285367870146E+137) {
    do {
      count++;
      fs_re *= 1.3435752215134178E-138;
      fs_im *= 1.3435752215134178E-138;
      gs_re *= 1.3435752215134178E-138;
      gs_im *= 1.3435752215134178E-138;
      scale *= 1.3435752215134178E-138;
    } while (!(scale < 7.4428285367870146E+137));

    rescaledir = 1;
    guard1 = TRUE;
  } else if (scale <= 1.3435752215134178E-138) {
    if ((g.re == 0.0) && (g.im == 0.0)) {
      *cs = 1.0;
      sn->re = 0.0;
      sn->im = 0.0;
      *r = f;
    } else {
      do {
        count++;
        fs_re *= 7.4428285367870146E+137;
        fs_im *= 7.4428285367870146E+137;
        gs_re *= 7.4428285367870146E+137;
        gs_im *= 7.4428285367870146E+137;
        scale *= 7.4428285367870146E+137;
      } while (!(scale > 1.3435752215134178E-138));

      rescaledir = -1;
      guard1 = TRUE;
    }
  } else {
    guard1 = TRUE;
  }

  if (guard1 == TRUE) {
    scale = fs_re * fs_re + fs_im * fs_im;
    g2 = gs_re * gs_re + gs_im * gs_im;
    f2s = g2;
    if (1.0 > g2) {
      f2s = 1.0;
    }

    if (scale <= f2s * 2.0041683600089728E-292) {
      if ((f.re == 0.0) && (f.im == 0.0)) {
        *cs = 0.0;
        r->re = rt_hypotd_snf(fabs(g.re), fabs(g.im));
        r->im = 0.0;
        f2s = rt_hypotd_snf(fabs(gs_re), fabs(gs_im));
        sn->re = gs_re / f2s;
        sn->im = -gs_im / f2s;
      } else {
        g2s = sqrt(g2);
        *cs = rt_hypotd_snf(fabs(fs_re), fabs(fs_im)) / g2s;
        f2s = fabs(f.re);
        g2 = fabs(f.im);
        if (g2 > f2s) {
          f2s = g2;
        }

        if (f2s > 1.0) {
          f2s = rt_hypotd_snf(fabs(f.re), fabs(f.im));
          fs_re = f.re / f2s;
          fs_im = f.im / f2s;
        } else {
          g2 = 7.4428285367870146E+137 * f.re;
          scale = 7.4428285367870146E+137 * f.im;
          f2s = rt_hypotd_snf(fabs(g2), fabs(scale));
          fs_re = g2 / f2s;
          fs_im = scale / f2s;
        }

        gs_re /= g2s;
        gs_im = -gs_im / g2s;
        sn->re = fs_re * gs_re - fs_im * gs_im;
        sn->im = fs_re * gs_im + fs_im * gs_re;
        r->re = *cs * f.re + (sn->re * g.re - sn->im * g.im);
        r->im = *cs * f.im + (sn->re * g.im + sn->im * g.re);
      }
    } else {
      f2s = sqrt(1.0 + g2 / scale);
      r->re = f2s * fs_re;
      r->im = f2s * fs_im;
      *cs = 1.0 / f2s;
      f2s = scale + g2;
      g2 = r->re / f2s;
      f2s = r->im / f2s;
      sn->re = g2 * gs_re - f2s * -gs_im;
      sn->im = g2 * -gs_im + f2s * gs_re;
      if (rescaledir > 0) {
        for (rescaledir = 1; rescaledir <= count; rescaledir++) {
          r->re *= 7.4428285367870146E+137;
          r->im *= 7.4428285367870146E+137;
        }
      } else {
        if (rescaledir < 0) {
          for (rescaledir = 1; rescaledir <= count; rescaledir++) {
            r->re *= 1.3435752215134178E-138;
            r->im *= 1.3435752215134178E-138;
          }
        }
      }
    }
  }
}

static void eml_matlab_ztgevc(const creal_T A[36], creal_T V[36])
{
  creal_T work1[6];
  creal_T work2[6];
  real_T rworka[6];
  int32_T i;
  real_T anorm;
  int32_T j;
  real_T y;
  real_T acoefa;
  real_T ascale;
  int32_T je;
  real_T temp;
  real_T salpha_re;
  real_T salpha_im;
  real_T acoeff;
  boolean_T b0;
  boolean_T b1;
  real_T scale;
  int32_T jr;
  creal_T d;
  creal_T b_work1;
  for (i = 0; i < 6; i++) {
    work1[i].re = 0.0;
    work1[i].im = 0.0;
    work2[i].re = 0.0;
    work2[i].im = 0.0;
    rworka[i] = 0.0;
  }

  anorm = fabs(A[0].re) + fabs(A[0].im);
  for (j = 0; j < 5; j++) {
    for (i = 0; i <= j; i++) {
      rworka[1 + j] += fabs(A[i + 6 * (1 + j)].re) + fabs(A[i + 6 * (1 + j)].im);
    }

    y = rworka[1 + j] + (fabs(A[(j + 6 * (1 + j)) + 1].re) + fabs(A[(j + 6 * (1
      + j)) + 1].im));
    if (y > anorm) {
      anorm = y;
    }
  }

  acoefa = anorm;
  if (2.2250738585072014E-308 > anorm) {
    acoefa = 2.2250738585072014E-308;
  }

  ascale = 1.0 / acoefa;
  for (je = 0; je < 6; je++) {
    y = (fabs(A[(6 * (5 - je) - je) + 5].re) + fabs(A[(6 * (5 - je) - je) + 5].
          im)) * ascale;
    if (1.0 > y) {
      y = 1.0;
    }

    temp = 1.0 / y;
    salpha_re = ascale * (temp * A[(6 * (5 - je) - je) + 5].re);
    salpha_im = ascale * (temp * A[(6 * (5 - je) - je) + 5].im);
    acoeff = temp * ascale;
    if ((fabs(temp) >= 2.2250738585072014E-308) && (fabs(acoeff) <
         6.0125050800269183E-292)) {
      b0 = TRUE;
    } else {
      b0 = FALSE;
    }

    if ((fabs(salpha_re) + fabs(salpha_im) >= 2.2250738585072014E-308) && (fabs
         (salpha_re) + fabs(salpha_im) < 6.0125050800269183E-292)) {
      b1 = TRUE;
    } else {
      b1 = FALSE;
    }

    scale = 1.0;
    if (b0) {
      acoefa = anorm;
      if (1.6632002579455998E+291 < anorm) {
        acoefa = 1.6632002579455998E+291;
      }

      scale = 6.0125050800269183E-292 / fabs(temp) * acoefa;
    }

    if (b1) {
      acoefa = 6.0125050800269183E-292 / (fabs(salpha_re) + fabs(salpha_im));
      if (acoefa > scale) {
        scale = acoefa;
      }
    }

    if (b0 || b1) {
      y = fabs(acoeff);
      acoefa = fabs(salpha_re) + fabs(salpha_im);
      if (1.0 > y) {
        y = 1.0;
      }

      if (acoefa > y) {
        y = acoefa;
      }

      y = 1.0 / (2.2250738585072014E-308 * y);
      if (y < scale) {
        scale = y;
      }

      if (b0) {
        acoeff = ascale * (scale * temp);
      } else {
        acoeff *= scale;
      }

      if (b1) {
        salpha_re *= scale;
        salpha_im *= scale;
      } else {
        salpha_re *= scale;
        salpha_im *= scale;
      }
    }

    acoefa = fabs(acoeff);
    for (jr = 0; jr < 6; jr++) {
      work1[jr].re = 0.0;
      work1[jr].im = 0.0;
    }

    work1[5 - je].re = 1.0;
    work1[5 - je].im = 0.0;
    scale = 2.2204460492503131E-16 * acoefa * anorm;
    y = 2.2204460492503131E-16 * (fabs(salpha_re) + fabs(salpha_im));
    if (y > scale) {
      scale = y;
    }

    if (2.2250738585072014E-308 > scale) {
      scale = 2.2250738585072014E-308;
    }

    for (jr = 0; jr <= 4 - je; jr++) {
      work1[jr].re = acoeff * A[jr + 6 * (5 - je)].re;
      work1[jr].im = acoeff * A[jr + 6 * (5 - je)].im;
    }

    work1[5 - je].re = 1.0;
    work1[5 - je].im = 0.0;
    for (j = 0; j <= 4 - je; j++) {
      i = (-je - j) + 4;
      d.re = acoeff * A[i + 6 * i].re - salpha_re;
      d.im = acoeff * A[i + 6 * i].im - salpha_im;
      if (fabs(d.re) + fabs(d.im) <= scale) {
        d.re = scale;
        d.im = 0.0;
      }

      if ((fabs(d.re) + fabs(d.im) < 1.0) && (fabs(work1[i].re) + fabs(work1[i].
            im) >= 7.4903880619263159E+306 * (fabs(d.re) + fabs(d.im)))) {
        temp = 1.0 / (fabs(work1[i].re) + fabs(work1[i].im));
        for (jr = 0; jr <= 5 - je; jr++) {
          work1[jr].re *= temp;
          work1[jr].im *= temp;
        }
      }

      b_work1.re = -work1[i].re;
      b_work1.im = -work1[i].im;
      work1[i] = b_eml_div(b_work1, d);
      if (i + 1 > 1) {
        if (fabs(work1[i].re) + fabs(work1[i].im) > 1.0) {
          temp = 1.0 / (fabs(work1[i].re) + fabs(work1[i].im));
          if (acoefa * rworka[i] >= 7.4903880619263159E+306 * temp) {
            for (jr = 0; jr <= 5 - je; jr++) {
              work1[jr].re *= temp;
              work1[jr].im *= temp;
            }
          }
        }

        d.re = acoeff * work1[i].re;
        d.im = acoeff * work1[i].im;
        for (jr = 0; jr <= i - 1; jr++) {
          work1[jr].re += d.re * A[jr + 6 * i].re - d.im * A[jr + 6 * i].im;
          work1[jr].im += d.re * A[jr + 6 * i].im + d.im * A[jr + 6 * i].re;
        }
      }
    }

    for (jr = 0; jr < 6; jr++) {
      work2[jr].re = 0.0;
      work2[jr].im = 0.0;
    }

    for (i = 0; i <= 5 - je; i++) {
      for (jr = 0; jr < 6; jr++) {
        acoefa = V[jr + 6 * i].re * work1[i].re - V[jr + 6 * i].im * work1[i].im;
        scale = V[jr + 6 * i].re * work1[i].im + V[jr + 6 * i].im * work1[i].re;
        work2[jr].re += acoefa;
        work2[jr].im += scale;
      }
    }

    acoefa = fabs(work2[0].re) + fabs(work2[0].im);
    for (jr = 0; jr < 5; jr++) {
      y = fabs(work2[1 + jr].re) + fabs(work2[1 + jr].im);
      if (y > acoefa) {
        acoefa = y;
      }
    }

    if (acoefa > 2.2250738585072014E-308) {
      temp = 1.0 / acoefa;
      for (jr = 0; jr < 6; jr++) {
        V[jr + 6 * (5 - je)].re = temp * work2[jr].re;
        V[jr + 6 * (5 - je)].im = temp * work2[jr].im;
      }
    } else {
      for (jr = 0; jr < 6; jr++) {
        V[jr + 6 * (5 - je)].re = 0.0;
        V[jr + 6 * (5 - je)].im = 0.0;
      }
    }
  }
}

void b_eig(const creal_T A[9], creal_T V[9], creal_T D[9])
{
  creal_T b_A[9];
  int32_T j;
  creal_T beta1[3];
  creal_T alpha1[3];
  real_T colnorm;
  int32_T coltop;
  real_T scale;
  real_T absxk;
  real_T t;
  real_T alpha1_im;
  for (j = 0; j < 9; j++) {
    b_A[j].re = A[j].re;
    b_A[j].im = A[j].im;
  }

  b_eml_matlab_zggev(b_A, &colnorm, alpha1, beta1, V);
  for (coltop = 0; coltop < 8; coltop += 3) {
    colnorm = 0.0;
    scale = 2.2250738585072014E-308;
    for (j = coltop; j + 1 <= coltop + 3; j++) {
      absxk = fabs(V[j].re);
      if (absxk > scale) {
        t = scale / absxk;
        colnorm = 1.0 + colnorm * t * t;
        scale = absxk;
      } else {
        t = absxk / scale;
        colnorm += t * t;
      }

      absxk = fabs(V[j].im);
      if (absxk > scale) {
        t = scale / absxk;
        colnorm = 1.0 + colnorm * t * t;
        scale = absxk;
      } else {
        t = absxk / scale;
        colnorm += t * t;
      }
    }

    colnorm = scale * sqrt(colnorm);
    for (j = coltop; j + 1 <= coltop + 3; j++) {
      V[j] = eml_div(V[j], colnorm);
    }
  }

  for (j = 0; j < 3; j++) {
    t = alpha1[j].re;
    alpha1_im = alpha1[j].im;
    if (beta1[j].im == 0.0) {
      if (alpha1[j].im == 0.0) {
        alpha1[j].re /= beta1[j].re;
        alpha1[j].im = 0.0;
      } else if (alpha1[j].re == 0.0) {
        alpha1[j].re = 0.0;
        alpha1[j].im = alpha1_im / beta1[j].re;
      } else {
        alpha1[j].re /= beta1[j].re;
        alpha1[j].im = alpha1_im / beta1[j].re;
      }
    } else if (beta1[j].re == 0.0) {
      if (alpha1[j].re == 0.0) {
        alpha1[j].re = alpha1[j].im / beta1[j].im;
        alpha1[j].im = 0.0;
      } else if (alpha1[j].im == 0.0) {
        alpha1[j].re = 0.0;
        alpha1[j].im = -(t / beta1[j].im);
      } else {
        alpha1[j].re = alpha1[j].im / beta1[j].im;
        alpha1[j].im = -(t / beta1[j].im);
      }
    } else {
      absxk = fabs(beta1[j].re);
      colnorm = fabs(beta1[j].im);
      if (absxk > colnorm) {
        colnorm = beta1[j].im / beta1[j].re;
        scale = beta1[j].re + colnorm * beta1[j].im;
        alpha1[j].re = (alpha1[j].re + colnorm * alpha1[j].im) / scale;
        alpha1[j].im = (alpha1_im - colnorm * t) / scale;
      } else if (colnorm == absxk) {
        colnorm = beta1[j].re > 0.0 ? 0.5 : -0.5;
        scale = beta1[j].im > 0.0 ? 0.5 : -0.5;
        alpha1[j].re = (alpha1[j].re * colnorm + alpha1[j].im * scale) / absxk;
        alpha1[j].im = (alpha1_im * colnorm - t * scale) / absxk;
      } else {
        colnorm = beta1[j].re / beta1[j].im;
        scale = beta1[j].im + colnorm * beta1[j].re;
        alpha1[j].re = (colnorm * alpha1[j].re + alpha1[j].im) / scale;
        alpha1[j].im = (colnorm * alpha1_im - t) / scale;
      }
    }
  }

  for (j = 0; j < 9; j++) {
    D[j].re = 0.0;
    D[j].im = 0.0;
  }

  for (j = 0; j < 3; j++) {
    D[j + 3 * j] = alpha1[j];
  }
}

creal_T b_eml_div(const creal_T x, const creal_T y)
{
  creal_T z;
  real_T brm;
  real_T bim;
  real_T d;
  if (y.im == 0.0) {
    if (x.im == 0.0) {
      z.re = x.re / y.re;
      z.im = 0.0;
    } else if (x.re == 0.0) {
      z.re = 0.0;
      z.im = x.im / y.re;
    } else {
      z.re = x.re / y.re;
      z.im = x.im / y.re;
    }
  } else if (y.re == 0.0) {
    if (x.re == 0.0) {
      z.re = x.im / y.im;
      z.im = 0.0;
    } else if (x.im == 0.0) {
      z.re = 0.0;
      z.im = -(x.re / y.im);
    } else {
      z.re = x.im / y.im;
      z.im = -(x.re / y.im);
    }
  } else {
    brm = fabs(y.re);
    bim = fabs(y.im);
    if (brm > bim) {
      bim = y.im / y.re;
      d = y.re + bim * y.im;
      z.re = (x.re + bim * x.im) / d;
      z.im = (x.im - bim * x.re) / d;
    } else if (bim == brm) {
      bim = y.re > 0.0 ? 0.5 : -0.5;
      d = y.im > 0.0 ? 0.5 : -0.5;
      z.re = (x.re * bim + x.im * d) / brm;
      z.im = (x.im * bim - x.re * d) / brm;
    } else {
      bim = y.re / y.im;
      d = y.im + bim * y.re;
      z.re = (bim * x.re + x.im) / d;
      z.im = (bim * x.im - x.re) / d;
    }
  }

  return z;
}

void eig(const real_T A[36], creal_T V[36], creal_T D[36])
{
  creal_T b_A[36];
  int32_T j;
  creal_T beta1[6];
  creal_T alpha1[6];
  real_T colnorm;
  int32_T coltop;
  real_T scale;
  real_T absxk;
  real_T t;
  real_T alpha1_im;
  for (j = 0; j < 36; j++) {
    b_A[j].re = A[j];
    b_A[j].im = 0.0;
  }

  eml_matlab_zggev(b_A, &colnorm, alpha1, beta1, V);
  for (coltop = 0; coltop < 32; coltop += 6) {
    colnorm = 0.0;
    scale = 2.2250738585072014E-308;
    for (j = coltop; j + 1 <= coltop + 6; j++) {
      absxk = fabs(V[j].re);
      if (absxk > scale) {
        t = scale / absxk;
        colnorm = 1.0 + colnorm * t * t;
        scale = absxk;
      } else {
        t = absxk / scale;
        colnorm += t * t;
      }

      absxk = fabs(V[j].im);
      if (absxk > scale) {
        t = scale / absxk;
        colnorm = 1.0 + colnorm * t * t;
        scale = absxk;
      } else {
        t = absxk / scale;
        colnorm += t * t;
      }
    }

    colnorm = scale * sqrt(colnorm);
    for (j = coltop; j + 1 <= coltop + 6; j++) {
      V[j] = eml_div(V[j], colnorm);
    }
  }

  for (j = 0; j < 6; j++) {
    t = alpha1[j].re;
    alpha1_im = alpha1[j].im;
    if (beta1[j].im == 0.0) {
      if (alpha1[j].im == 0.0) {
        alpha1[j].re /= beta1[j].re;
        alpha1[j].im = 0.0;
      } else if (alpha1[j].re == 0.0) {
        alpha1[j].re = 0.0;
        alpha1[j].im = alpha1_im / beta1[j].re;
      } else {
        alpha1[j].re /= beta1[j].re;
        alpha1[j].im = alpha1_im / beta1[j].re;
      }
    } else if (beta1[j].re == 0.0) {
      if (alpha1[j].re == 0.0) {
        alpha1[j].re = alpha1[j].im / beta1[j].im;
        alpha1[j].im = 0.0;
      } else if (alpha1[j].im == 0.0) {
        alpha1[j].re = 0.0;
        alpha1[j].im = -(t / beta1[j].im);
      } else {
        alpha1[j].re = alpha1[j].im / beta1[j].im;
        alpha1[j].im = -(t / beta1[j].im);
      }
    } else {
      absxk = fabs(beta1[j].re);
      colnorm = fabs(beta1[j].im);
      if (absxk > colnorm) {
        colnorm = beta1[j].im / beta1[j].re;
        scale = beta1[j].re + colnorm * beta1[j].im;
        alpha1[j].re = (alpha1[j].re + colnorm * alpha1[j].im) / scale;
        alpha1[j].im = (alpha1_im - colnorm * t) / scale;
      } else if (colnorm == absxk) {
        colnorm = beta1[j].re > 0.0 ? 0.5 : -0.5;
        scale = beta1[j].im > 0.0 ? 0.5 : -0.5;
        alpha1[j].re = (alpha1[j].re * colnorm + alpha1[j].im * scale) / absxk;
        alpha1[j].im = (alpha1_im * colnorm - t * scale) / absxk;
      } else {
        colnorm = beta1[j].re / beta1[j].im;
        scale = beta1[j].im + colnorm * beta1[j].re;
        alpha1[j].re = (colnorm * alpha1[j].re + alpha1[j].im) / scale;
        alpha1[j].im = (colnorm * alpha1_im - t) / scale;
      }
    }
  }

  for (j = 0; j < 36; j++) {
    D[j].re = 0.0;
    D[j].im = 0.0;
  }

  for (j = 0; j < 6; j++) {
    D[j + 6 * j] = alpha1[j];
  }
}

creal_T eml_div(const creal_T x, real_T y)
{
  creal_T z;
  if (x.im == 0.0) {
    z.re = x.re / y;
    z.im = 0.0;
  } else if (x.re == 0.0) {
    z.re = 0.0;
    z.im = x.im / y;
  } else {
    z.re = x.re / y;
    z.im = x.im / y;
  }

  return z;
}

/* End of code generation (eig.c) */
