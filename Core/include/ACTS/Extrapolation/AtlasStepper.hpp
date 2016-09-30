#ifndef ACTS_ATLAS_STEPPER_HPP
#define ACTS_ATLAS_STEPPER_HPP 1

#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Extrapolation/Direction.hpp"
#include "ACTS/Utilities/Units.hpp"

namespace Acts {

template <typename BField>
class AtlasStepper
{
  struct Cache
  {
    // configuration
    double direction;
    bool   useJacobian;
    double step;
    double maxPathLength;
    bool   mcondition;
    bool   needgradient;
    bool   newfield;
    // internal parameters to be used
    double field[3] = {0., 0., 0.};
    double pVector[64];
    // result
    double             parameters[5] = {0., 0., 0., 0., 0.};
    ActsSymMatrixD<5>* covariance;
    double             jacobian[25];

    Cache(const CurvilinearParameters& pars, Direction dir)
      : direction(alongMomentum)
      , useJacobian(false)
      , step(0.)
      , maxPathLength(0.)
      , mcondition(false)
      , needgradient(false)
      , newfield(true)
      , covariance(nullptr)
    {
      if (pars.covariance())
        covariance = new ActsSymMatrixD<5>(*pars.covariance());

      const ActsVectorD<3> pos = pars.position();
      ActsVectorD<5>       Vp  = pars.parameters();

      // invert direction for backward propagation
      if (dir == backward) {
        // phi
        Vp(2) += M_PI;
        if (Vp(2) > M_PI) Vp(2) -= 2 * M_PI;
        Vp(3) = M_PI - Vp(3);
      }
      double Sf, Cf, Ce, Se;
      sincos(Vp(2), &Sf, &Cf);
      sincos(Vp(3), &Se, &Ce);

      double Ax[3] = {-Sf, Cf, 0.};
      double Ay[3] = {-Cf * Ce, -Sf * Ce, Se};

      pVector[0] = pos(0);
      pVector[1] = pos(1);
      pVector[2] = pos(2);
      pVector[3] = Cf * Se;      // Ax
      pVector[4] = Sf * Se;      // Ay
      pVector[5] = Ce;           // Az
      pVector[6] = dir * Vp[4];  // CM
      if (fabs(pVector[6]) < .000000000000001) {
        pVector[6] < 0. ? pVector[6] = -.000000000000001 : pVector[6]
            = .000000000000001;
      }

      //   /dL1     |   /dL2       |    /dPhi     |    /dThe     |    /dCM     |
      //
      pVector[7]  = Ax[0];
      pVector[14] = Ay[0];
      pVector[21] = 0.;
      pVector[28] = 0.;
      pVector[35] = 0.;  // dX /
      pVector[8]  = Ax[1];
      pVector[15] = Ay[1];
      pVector[22] = 0.;
      pVector[29] = 0.;
      pVector[36] = 0.;  // dY /
      pVector[9]  = Ax[2];
      pVector[16] = Ay[2];
      pVector[23] = 0.;
      pVector[30] = 0.;
      pVector[37] = 0.;  // dZ /
      pVector[10] = 0.;
      pVector[17] = 0.;
      pVector[24] = -Sf * Se;
      pVector[31] = -Ay[0];
      pVector[38] = 0.;  // dAx/
      pVector[11] = 0.;
      pVector[18] = 0.;
      pVector[25] = Cf * Se;
      pVector[32] = -Ay[1];
      pVector[39] = 0.;  // dAy/
      pVector[12] = 0.;
      pVector[19] = 0.;
      pVector[26] = 0.;
      pVector[33] = -Ay[2];
      pVector[40] = 0.;  // dAz/
      pVector[13] = 0.;
      pVector[20] = 0.;
      pVector[27] = 0.;
      pVector[34] = 0.;
      pVector[41] = 1.;  // dCM/
      pVector[42] = 0.;
      pVector[43] = 0.;
      pVector[44] = 0.;
    }
  };

public:
  template <typename T>
  using cache_type = Cache;

  template <typename T>
  using return_parameter_type = CurvilinearParameters;

  static CurvilinearParameters
  convert(const Cache& cache, Direction dir)
  {
    double         charge = cache.pVector[6] > 0. ? dir : -dir;
    Acts::Vector3D gp(cache.pVector[0], cache.pVector[1], cache.pVector[2]);
    Acts::Vector3D mom(cache.pVector[3], cache.pVector[4], cache.pVector[5]);
    mom /= fabs(cache.pVector[6]) * dir;

    return CurvilinearParameters(
        std::unique_ptr<ActsSymMatrixD<5>>(cache.covariance), gp, mom, charge);
  }

  AtlasStepper(BField&& bField = BField()) : m_bField(std::move(bField)){};

  double
  step_backward(Cache& cache, double& stepMax) const
  {
    return step_forward(cache, stepMax);
  }

  double
  step_forward(Cache& cache, double& stepMax) const
  {
    bool Jac = cache.useJacobian;

    double* R  = &(cache.pVector[0]);  // Coordinates
    double* A  = &(cache.pVector[3]);  // Directions
    double* sA = &(cache.pVector[42]);
    // Invert mometum/2.
    double Pi   = 0.5 / units::Nat2SI<units::MOMENTUM>(1. / cache.pVector[6]);
    double dltm = 0.0002 * .03;

    double f0[3], f[3];

    // if new field is required get it
    if (cache.newfield)
      m_bField(R, f0);
    else {
      f0[0] = cache.field[0];
      f0[1] = cache.field[1];
      f0[2] = cache.field[2];
    }

    bool Helix = false;
    // if (fabs(S) < m_cfg.helixStep) Helix = true;

    while (stepMax != 0.) {
      double S3 = (1. / 3.) * stepMax, S4 = .25 * stepMax, PS2 = Pi * stepMax;

      // First point
      //
      double H0[3] = {f0[0] * PS2, f0[1] * PS2, f0[2] * PS2};
      double A0    = A[1] * H0[2] - A[2] * H0[1];
      double B0    = A[2] * H0[0] - A[0] * H0[2];
      double C0    = A[0] * H0[1] - A[1] * H0[0];
      double A2    = A0 + A[0];
      double B2    = B0 + A[1];
      double C2    = C0 + A[2];
      double A1    = A2 + A[0];
      double B1    = B2 + A[1];
      double C1    = C2 + A[2];

      // Second point
      //
      if (!Helix) {
        double gP[3] = {R[0] + A1 * S4, R[1] + B1 * S4, R[2] + C1 * S4};
        m_bField(gP, f);
      } else {
        f[0] = f0[0];
        f[1] = f0[1];
        f[2] = f0[2];
      }

      double H1[3] = {f[0] * PS2, f[1] * PS2, f[2] * PS2};
      double A3    = (A[0] + B2 * H1[2]) - C2 * H1[1];
      double B3    = (A[1] + C2 * H1[0]) - A2 * H1[2];
      double C3    = (A[2] + A2 * H1[1]) - B2 * H1[0];
      double A4    = (A[0] + B3 * H1[2]) - C3 * H1[1];
      double B4    = (A[1] + C3 * H1[0]) - A3 * H1[2];
      double C4    = (A[2] + A3 * H1[1]) - B3 * H1[0];
      double A5    = 2. * A4 - A[0];
      double B5    = 2. * B4 - A[1];
      double C5    = 2. * C4 - A[2];

      // Last point
      //
      if (!Helix) {
        double gP[3]
            = {R[0] + stepMax * A4, R[1] + stepMax * B4, R[2] + stepMax * C4};
        m_bField(gP, f);
      } else {
        f[0] = f0[0];
        f[1] = f0[1];
        f[2] = f0[2];
      }

      double H2[3] = {f[0] * PS2, f[1] * PS2, f[2] * PS2};
      double A6    = B5 * H2[2] - C5 * H2[1];
      double B6    = C5 * H2[0] - A5 * H2[2];
      double C6    = A5 * H2[1] - B5 * H2[0];

      // Test approximation quality on give step and possible step reduction
      //
      double EST = fabs((A1 + A6) - (A3 + A4)) + fabs((B1 + B6) - (B3 + B4))
          + fabs((C1 + C6) - (C3 + C4));
      if (EST > 0.0002) {
        stepMax *= .5;
        dltm = 0.;
        continue;
      }

      // Parameters calculation
      //
      double A00 = A[0], A11 = A[1], A22 = A[2];

      A[0] = 2. * A3 + (A0 + A5 + A6);
      A[1] = 2. * B3 + (B0 + B5 + B6);
      A[2] = 2. * C3 + (C0 + C5 + C6);

      double D  = (A[0] * A[0] + A[1] * A[1]) + (A[2] * A[2] - 9.);
      double Sl = 2. / stepMax;
      D         = (1. / 3.) - ((1. / 648.) * D) * (12. - D);

      R[0] += (A2 + A3 + A4) * S3;
      R[1] += (B2 + B3 + B4) * S3;
      R[2] += (C2 + C3 + C4) * S3;
      A[0] *= D;
      A[1] *= D;
      A[2] *= D;
      sA[0] = A6 * Sl;
      sA[1] = B6 * Sl;
      sA[2] = C6 * Sl;

      cache.field[0] = f[0];
      cache.field[1] = f[1];
      cache.field[2] = f[2];
      cache.newfield = false;

      // stepMax *= 2;
      if (!Jac) return stepMax;

      // Jacobian calculation
      //
      double* d2A  = &cache.pVector[24];
      double* d3A  = &cache.pVector[31];
      double* d4A  = &cache.pVector[38];
      double  d2A0 = H0[2] * d2A[1] - H0[1] * d2A[2];
      double  d2B0 = H0[0] * d2A[2] - H0[2] * d2A[0];
      double  d2C0 = H0[1] * d2A[0] - H0[0] * d2A[1];
      double  d3A0 = H0[2] * d3A[1] - H0[1] * d3A[2];
      double  d3B0 = H0[0] * d3A[2] - H0[2] * d3A[0];
      double  d3C0 = H0[1] * d3A[0] - H0[0] * d3A[1];
      double  d4A0 = (A0 + H0[2] * d4A[1]) - H0[1] * d4A[2];
      double  d4B0 = (B0 + H0[0] * d4A[2]) - H0[2] * d4A[0];
      double  d4C0 = (C0 + H0[1] * d4A[0]) - H0[0] * d4A[1];
      double  d2A2 = d2A0 + d2A[0];
      double  d2B2 = d2B0 + d2A[1];
      double  d2C2 = d2C0 + d2A[2];
      double  d3A2 = d3A0 + d3A[0];
      double  d3B2 = d3B0 + d3A[1];
      double  d3C2 = d3C0 + d3A[2];
      double  d4A2 = d4A0 + d4A[0];
      double  d4B2 = d4B0 + d4A[1];
      double  d4C2 = d4C0 + d4A[2];
      double  d0   = d4A[0] - A00;
      double  d1   = d4A[1] - A11;
      double  d2   = d4A[2] - A22;
      double  d2A3 = (d2A[0] + d2B2 * H1[2]) - d2C2 * H1[1];
      double  d2B3 = (d2A[1] + d2C2 * H1[0]) - d2A2 * H1[2];
      double  d2C3 = (d2A[2] + d2A2 * H1[1]) - d2B2 * H1[0];
      double  d3A3 = (d3A[0] + d3B2 * H1[2]) - d3C2 * H1[1];
      double  d3B3 = (d3A[1] + d3C2 * H1[0]) - d3A2 * H1[2];
      double  d3C3 = (d3A[2] + d3A2 * H1[1]) - d3B2 * H1[0];
      double  d4A3 = ((A3 + d0) + d4B2 * H1[2]) - d4C2 * H1[1];
      double  d4B3 = ((B3 + d1) + d4C2 * H1[0]) - d4A2 * H1[2];
      double  d4C3 = ((C3 + d2) + d4A2 * H1[1]) - d4B2 * H1[0];
      double  d2A4 = (d2A[0] + d2B3 * H1[2]) - d2C3 * H1[1];
      double  d2B4 = (d2A[1] + d2C3 * H1[0]) - d2A3 * H1[2];
      double  d2C4 = (d2A[2] + d2A3 * H1[1]) - d2B3 * H1[0];
      double  d3A4 = (d3A[0] + d3B3 * H1[2]) - d3C3 * H1[1];
      double  d3B4 = (d3A[1] + d3C3 * H1[0]) - d3A3 * H1[2];
      double  d3C4 = (d3A[2] + d3A3 * H1[1]) - d3B3 * H1[0];
      double  d4A4 = ((A4 + d0) + d4B3 * H1[2]) - d4C3 * H1[1];
      double  d4B4 = ((B4 + d1) + d4C3 * H1[0]) - d4A3 * H1[2];
      double  d4C4 = ((C4 + d2) + d4A3 * H1[1]) - d4B3 * H1[0];
      double  d2A5 = 2. * d2A4 - d2A[0];
      double  d2B5 = 2. * d2B4 - d2A[1];
      double  d2C5 = 2. * d2C4 - d2A[2];
      double  d3A5 = 2. * d3A4 - d3A[0];
      double  d3B5 = 2. * d3B4 - d3A[1];
      double  d3C5 = 2. * d3C4 - d3A[2];
      double  d4A5 = 2. * d4A4 - d4A[0];
      double  d4B5 = 2. * d4B4 - d4A[1];
      double  d4C5 = 2. * d4C4 - d4A[2];
      double  d2A6 = d2B5 * H2[2] - d2C5 * H2[1];
      double  d2B6 = d2C5 * H2[0] - d2A5 * H2[2];
      double  d2C6 = d2A5 * H2[1] - d2B5 * H2[0];
      double  d3A6 = d3B5 * H2[2] - d3C5 * H2[1];
      double  d3B6 = d3C5 * H2[0] - d3A5 * H2[2];
      double  d3C6 = d3A5 * H2[1] - d3B5 * H2[0];
      double  d4A6 = d4B5 * H2[2] - d4C5 * H2[1];
      double  d4B6 = d4C5 * H2[0] - d4A5 * H2[2];
      double  d4C6 = d4A5 * H2[1] - d4B5 * H2[0];

      double* dR = &cache.pVector[21];
      dR[0] += (d2A2 + d2A3 + d2A4) * S3;
      dR[1] += (d2B2 + d2B3 + d2B4) * S3;
      dR[2] += (d2C2 + d2C3 + d2C4) * S3;
      d2A[0] = ((d2A0 + 2. * d2A3) + (d2A5 + d2A6)) * (1. / 3.);
      d2A[1] = ((d2B0 + 2. * d2B3) + (d2B5 + d2B6)) * (1. / 3.);
      d2A[2] = ((d2C0 + 2. * d2C3) + (d2C5 + d2C6)) * (1. / 3.);

      dR = &cache.pVector[28];
      dR[0] += (d3A2 + d3A3 + d3A4) * S3;
      dR[1] += (d3B2 + d3B3 + d3B4) * S3;
      dR[2] += (d3C2 + d3C3 + d3C4) * S3;
      d3A[0] = ((d3A0 + 2. * d3A3) + (d3A5 + d3A6)) * (1. / 3.);
      d3A[1] = ((d3B0 + 2. * d3B3) + (d3B5 + d3B6)) * (1. / 3.);
      d3A[2] = ((d3C0 + 2. * d3C3) + (d3C5 + d3C6)) * (1. / 3.);

      dR = &cache.pVector[35];
      dR[0] += (d4A2 + d4A3 + d4A4) * S3;
      dR[1] += (d4B2 + d4B3 + d4B4) * S3;
      dR[2] += (d4C2 + d4C3 + d4C4) * S3;
      d4A[0] = ((d4A0 + 2. * d4A3) + (d4A5 + d4A6 + A6)) * (1. / 3.);
      d4A[1] = ((d4B0 + 2. * d4B3) + (d4B5 + d4B6 + B6)) * (1. / 3.);
      d4A[2] = ((d4C0 + 2. * d4C3) + (d4C5 + d4C6 + C6)) * (1. / 3.);
      return stepMax;
    }

    return stepMax;
  }

private:
  BField m_bField;
};

}  // namespace Acts
#endif  // ACTS_ATLAS_STEPPER_HPP
