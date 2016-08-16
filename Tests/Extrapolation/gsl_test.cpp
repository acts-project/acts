#include <cmath>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <stdio.h>

void
bfield(const double pos[], double* B)
{
  double y = pos[1];

  if (fabs(y) < 10 or fabs(y) > 20)
    B[2] = 1;
  else
    B[2] = 0;
}

int
func(double t, const double y[], double f[], void* params)
{
  (void)(t); /* avoid unused parameter warning */
  double B[3] = {0, 0, 0};
  bfield(y, B);
  double Bx  = B[0];
  double By  = B[1];
  double Bz  = B[2];
  double m   = ((double*)params)[4];
  double qOm = ((double*)params)[3] / ((double*)params)[4];
  double _x  = y[0];
  double _y  = y[1];
  double _z  = y[2];
  double px  = y[3];
  double py  = y[4];
  double pz  = y[5];
  f[0]       = px / m;
  f[1]       = py / m;
  f[2]       = pz / m;
  f[3]       = qOm * (py * Bz - pz * By);
  f[4]       = qOm * (pz * Bx - px * Bz);
  f[5]       = qOm * (px * By - py * Bx);
  return GSL_SUCCESS;
}

int
jac(double t, const double y[], double* dfdy, double dfdt[], void* params)
{
  (void)(t); /* avoid unused parameter warning */
  double B[3] = {0, 0, 0};
  bfield(y, B);
  double Bx  = B[0];
  double By  = B[1];
  double Bz  = B[2];
  double m   = ((double*)params)[4];
  double qOm = ((double*)params)[3] / ((double*)params)[4];
  double _x  = y[0];
  double _y  = y[1];
  double _z  = y[2];
  double px  = y[3];
  double py  = y[4];
  double pz  = y[5];

  gsl_matrix_view dfdy_mat = gsl_matrix_view_array(dfdy, 6, 6);
  gsl_matrix*     mat      = &dfdy_mat.matrix;
  gsl_matrix_set(mat, 0, 0, 0.0);
  gsl_matrix_set(mat, 0, 1, 0.0);
  gsl_matrix_set(mat, 0, 2, 0.0);
  gsl_matrix_set(mat, 0, 3, 1 / m);
  gsl_matrix_set(mat, 0, 4, 0.0);
  gsl_matrix_set(mat, 0, 5, 0.0);

  gsl_matrix_set(mat, 1, 0, 0.0);
  gsl_matrix_set(mat, 1, 1, 0.0);
  gsl_matrix_set(mat, 1, 2, 0.0);
  gsl_matrix_set(mat, 1, 3, 0.0);
  gsl_matrix_set(mat, 1, 4, 1 / m);
  gsl_matrix_set(mat, 1, 5, 0.0);

  gsl_matrix_set(mat, 2, 0, 0.0);
  gsl_matrix_set(mat, 2, 1, 0.0);
  gsl_matrix_set(mat, 2, 2, 0.0);
  gsl_matrix_set(mat, 2, 3, 0.0);
  gsl_matrix_set(mat, 2, 4, 0.0);
  gsl_matrix_set(mat, 2, 5, 1 / m);

  gsl_matrix_set(mat, 3, 0, 0.0);
  gsl_matrix_set(mat, 3, 1, 0.0);
  gsl_matrix_set(mat, 3, 2, 0.0);
  gsl_matrix_set(mat, 3, 3, 0.0);
  gsl_matrix_set(mat, 3, 4, qOm * Bz);
  gsl_matrix_set(mat, 3, 5, -qOm * By);

  gsl_matrix_set(mat, 4, 0, 0.0);
  gsl_matrix_set(mat, 4, 1, 0.0);
  gsl_matrix_set(mat, 4, 2, 0.0);
  gsl_matrix_set(mat, 4, 3, -qOm * Bz);
  gsl_matrix_set(mat, 4, 4, 0.0);
  gsl_matrix_set(mat, 4, 5, qOm * Bx);

  gsl_matrix_set(mat, 5, 0, 0.0);
  gsl_matrix_set(mat, 5, 1, 0.0);
  gsl_matrix_set(mat, 5, 2, 0.0);
  gsl_matrix_set(mat, 5, 3, qOm * By);
  gsl_matrix_set(mat, 5, 4, -qOm * Bx);
  gsl_matrix_set(mat, 5, 5, 0.0);

  dfdt[0] = 0;
  dfdt[1] = 0;
  dfdt[2] = 0;
  dfdt[3] = 0;
  dfdt[4] = 0;
  dfdt[5] = 0;

  return GSL_SUCCESS;
}

int
main(void)
{
  double            Bx = 0, By = 0, Bz = 1;
  double            q         = 1;
  double            m         = 0.146;
  double            params[5] = {Bx, By, Bz, q, m};
  gsl_odeiv2_system sys       = {func, jac, 6, params};

  gsl_odeiv2_driver* d = gsl_odeiv2_driver_alloc_y_new(
      &sys, gsl_odeiv2_step_rkck, 1e-9, 1e-6, 0.0);
  int    i;
  double t = 0.0, t1 = 1.0;
  double y[6] = {0, 0, 0, 10.0, 0.0, 0.0};

  for (i = 1; i <= 100; i++) {
    double ti     = i * t1 / 100.0;
    int    status = gsl_odeiv2_driver_apply(d, &t, ti, y);

    if (status != GSL_SUCCESS) {
      printf("error, return value=%d\n", status);
      break;
    }
    //    printf("B = (%.2f %.2f %.2f), q = %.1f, m = %.3f\n",
    //           params[0],
    //           params[1],
    //           params[2],
    //           params[3],
    //           params[4]);
    printf("%.3e %.3e %.3e\n", y[0], y[1], y[2]);
    //    printf("p = (%.3e %.3e %.3e)\n", y[3], y[4], y[5]);
    //    printf("pT = %.3e, phi = %.3e \n",
    //           sqrt(y[3] * y[3] + y[4] * y[4]),
    //           atan2(y[4], y[3]));
  }

  gsl_odeiv2_driver_free(d);
  return 0;
}
