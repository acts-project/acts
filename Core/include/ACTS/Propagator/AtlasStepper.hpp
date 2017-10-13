// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_ATLAS_STEPPER_HPP
#define ACTS_ATLAS_STEPPER_HPP 1

#include <cmath>
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/Units.hpp"

namespace Acts {

template <typename BField>
class AtlasStepper
{
  struct Cache
  {
    // configuration
    double dir;
    bool   useJacobian;
    double step;
    double maxPathLength;
    bool   mcondition;
    bool   needgradient;
    bool   newfield;
    // internal parameters to be used
    Vector3D field;
    double   pVector[64];
    // result
    double parameters[NGlobalPars] = {0., 0., 0., 0., 0.};
    const ActsSymMatrixD<NGlobalPars>* covariance;
    double                             jacobian[NGlobalPars * NGlobalPars];

    Vector3D
    position() const
    {
      return Vector3D(pVector[0], pVector[1], pVector[2]);
    }

    Vector3D
    direction() const
    {
      return Vector3D(pVector[3], pVector[4], pVector[5]);
    }

    Cache(const CurvilinearParameters& pars)
      : dir(alongMomentum)
      , useJacobian(false)
      , step(0.)
      , maxPathLength(0.)
      , mcondition(false)
      , needgradient(false)
      , newfield(true)
      , field(0., 0., 0.)
      , covariance(nullptr)
    {
      if (pars.covariance()) {
        covariance  = new ActsSymMatrixD<NGlobalPars>(*pars.covariance());
        useJacobian = true;
      }

      const ActsVectorD<3>     pos = pars.position();
      ActsVectorD<NGlobalPars> Vp  = pars.parameters();

      double Sf, Cf, Ce, Se;
      Sf = sin(Vp(2));
      Cf = cos(Vp(2));
      Se = sin(Vp(3));
      Ce = cos(Vp(3));

      double Ax[3] = {-Sf, Cf, 0.};
      double Ay[3] = {-Cf * Ce, -Sf * Ce, Se};

      pVector[0] = pos(0);
      pVector[1] = pos(1);
      pVector[2] = pos(2);
      pVector[3] = Cf * Se;  // Ax
      pVector[4] = Sf * Se;  // Ay
      pVector[5] = Ce;       // Az
      pVector[6] = Vp[4];    // CM
      if (std::abs(pVector[6]) < .000000000000001) {
        pVector[6] < 0. ? pVector[6] = -.000000000000001
                        : pVector[6] = .000000000000001;
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
  template <typename T, typename S = int>
  using cache_type = Cache;

  template <typename T>
  using step_parameter_type = CurvilinearParameters;

  template <typename T, typename S>
  struct s
  {
    typedef BoundParameters type;
  };

  template <typename T>
  struct s<T, int>
  {
    typedef CurvilinearParameters type;
  };

  template <typename T, typename S = int>
  using return_parameter_type = typename s<T, S>::type;

  static CurvilinearParameters
  convert(Cache& cache)
  {
    double         charge = cache.pVector[6] > 0. ? 1. : -1.;
    Acts::Vector3D gp(cache.pVector[0], cache.pVector[1], cache.pVector[2]);
    Acts::Vector3D mom(cache.pVector[3], cache.pVector[4], cache.pVector[5]);
    mom /= std::abs(cache.pVector[6]);

    double P[45];
    for (unsigned int i = 0; i < 45; ++i) P[i] = cache.pVector[i];

    std::unique_ptr<const ActsSymMatrixD<NGlobalPars>> cov = nullptr;
    if (cache.covariance) {
      double p = 1. / P[6];
      P[35] *= p;
      P[36] *= p;
      P[37] *= p;
      P[38] *= p;
      P[39] *= p;
      P[40] *= p;

      double An = sqrt(P[3] * P[3] + P[4] * P[4]);
      double Ax[3];
      if (An != 0.) {
        Ax[0] = -P[4] / An;
        Ax[1] = P[3] / An;
        Ax[2] = 0.;
      } else {
        Ax[0] = 1.;
        Ax[1] = 0.;
        Ax[2] = 0.;
      }

      double Ay[3] = {-Ax[1] * P[5], Ax[0] * P[5], An};
      double S[3]  = {P[3], P[4], P[5]};

      double A       = P[3] * S[0] + P[4] * S[1] + P[5] * S[2];
      if (A != 0.) A = 1. / A;
      S[0] *= A;
      S[1] *= A;
      S[2] *= A;

      double s0 = P[7] * S[0] + P[8] * S[1] + P[9] * S[2];
      double s1 = P[14] * S[0] + P[15] * S[1] + P[16] * S[2];
      double s2 = P[21] * S[0] + P[22] * S[1] + P[23] * S[2];
      double s3 = P[28] * S[0] + P[29] * S[1] + P[30] * S[2];
      double s4 = P[35] * S[0] + P[36] * S[1] + P[37] * S[2];

      P[7] -= (s0 * P[3]);
      P[8] -= (s0 * P[4]);
      P[9] -= (s0 * P[5]);
      P[10] -= (s0 * P[42]);
      P[11] -= (s0 * P[43]);
      P[12] -= (s0 * P[44]);
      P[14] -= (s1 * P[3]);
      P[15] -= (s1 * P[4]);
      P[16] -= (s1 * P[5]);
      P[17] -= (s1 * P[42]);
      P[18] -= (s1 * P[43]);
      P[19] -= (s1 * P[44]);
      P[21] -= (s2 * P[3]);
      P[22] -= (s2 * P[4]);
      P[23] -= (s2 * P[5]);
      P[24] -= (s2 * P[42]);
      P[25] -= (s2 * P[43]);
      P[26] -= (s2 * P[44]);
      P[28] -= (s3 * P[3]);
      P[29] -= (s3 * P[4]);
      P[30] -= (s3 * P[5]);
      P[31] -= (s3 * P[42]);
      P[32] -= (s3 * P[43]);
      P[33] -= (s3 * P[44]);
      P[35] -= (s4 * P[3]);
      P[36] -= (s4 * P[4]);
      P[37] -= (s4 * P[5]);
      P[38] -= (s4 * P[42]);
      P[39] -= (s4 * P[43]);
      P[40] -= (s4 * P[44]);

      double P3, P4, C = P[3] * P[3] + P[4] * P[4];
      if (C > 1.e-20) {
        C  = 1. / C;
        P3 = P[3] * C;
        P4 = P[4] * C;
        C  = -sqrt(C);
      } else {
        C  = -1.e10;
        P3 = 1.;
        P4 = 0.;
      }

      // Jacobian production
      //
      cache.jacobian[0] = Ax[0] * P[7] + Ax[1] * P[8];    // dL0/dL0
      cache.jacobian[1] = Ax[0] * P[14] + Ax[1] * P[15];  // dL0/dL1
      cache.jacobian[2] = Ax[0] * P[21] + Ax[1] * P[22];  // dL0/dPhi
      cache.jacobian[3] = Ax[0] * P[28] + Ax[1] * P[29];  // dL0/dThe
      cache.jacobian[4] = Ax[0] * P[35] + Ax[1] * P[36];  // dL0/dCM
      cache.jacobian[5]
          = Ay[0] * P[7] + Ay[1] * P[8] + Ay[2] * P[9];  // dL1/dL0
      cache.jacobian[6]
          = Ay[0] * P[14] + Ay[1] * P[15] + Ay[2] * P[16];  // dL1/dL1
      cache.jacobian[7]
          = Ay[0] * P[21] + Ay[1] * P[22] + Ay[2] * P[23];  // dL1/dPhi
      cache.jacobian[8]
          = Ay[0] * P[28] + Ay[1] * P[29] + Ay[2] * P[30];  // dL1/dThe
      cache.jacobian[9]
          = Ay[0] * P[35] + Ay[1] * P[36] + Ay[2] * P[37];  // dL1/dCM
      cache.jacobian[10] = P3 * P[11] - P4 * P[10];         // dPhi/dL0
      cache.jacobian[11] = P3 * P[18] - P4 * P[17];         // dPhi/dL1
      cache.jacobian[12] = P3 * P[25] - P4 * P[24];         // dPhi/dPhi
      cache.jacobian[13] = P3 * P[32] - P4 * P[31];         // dPhi/dThe
      cache.jacobian[14] = P3 * P[39] - P4 * P[38];         // dPhi/dCM
      cache.jacobian[15] = C * P[12];                       // dThe/dL0
      cache.jacobian[16] = C * P[19];                       // dThe/dL1
      cache.jacobian[17] = C * P[26];                       // dThe/dPhi
      cache.jacobian[18] = C * P[33];                       // dThe/dThe
      cache.jacobian[19] = C * P[40];                       // dThe/dCM
      cache.jacobian[20] = 0;                               // dCM /dL0
      cache.jacobian[21] = 0;                               // dCM /dL1
      cache.jacobian[22] = 0;                               // dCM /dPhi
      cache.jacobian[23] = 0;                               // dCM /dTheta
      cache.jacobian[24] = P[41];                           // dCM /dCM
      Eigen::
          Map<Eigen::Matrix<double, NGlobalPars, NGlobalPars, Eigen::RowMajor>>
              J(cache.jacobian);

      cov = std::make_unique<const ActsSymMatrixD<NGlobalPars>>(
          J * (*cache.covariance) * J.transpose());
    }

    return CurvilinearParameters(std::move(cov), gp, mom, charge);
  }

  static BoundParameters
  convert(Cache& cache, const Surface& s)
  {
    double         charge = cache.pVector[6] > 0. ? 1. : -1.;
    Acts::Vector3D gp(cache.pVector[0], cache.pVector[1], cache.pVector[2]);
    Acts::Vector3D mom(cache.pVector[3], cache.pVector[4], cache.pVector[5]);
    mom /= std::abs(cache.pVector[6]);

    std::unique_ptr<const ActsSymMatrixD<5>> cov = nullptr;
    if (cache.covariance) {
      double p = 1. / cache.pVector[6];
      cache.pVector[35] *= p;
      cache.pVector[36] *= p;
      cache.pVector[37] *= p;
      cache.pVector[38] *= p;
      cache.pVector[39] *= p;
      cache.pVector[40] *= p;

      double An = sqrt(cache.pVector[3] * cache.pVector[3]
                       + cache.pVector[4] * cache.pVector[4]);
      double Ax[3];
      if (An != 0.) {
        Ax[0] = -cache.pVector[4] / An;
        Ax[1] = cache.pVector[3] / An;
        Ax[2] = 0.;
      } else {
        Ax[0] = 1.;
        Ax[1] = 0.;
        Ax[2] = 0.;
      }

      double Ay[3] = {-Ax[1] * cache.pVector[5], Ax[0] * cache.pVector[5], An};
      double S[3]  = {cache.pVector[3], cache.pVector[4], cache.pVector[5]};

      double A = cache.pVector[3] * S[0] + cache.pVector[4] * S[1]
          + cache.pVector[5] * S[2];
      if (A != 0.) A = 1. / A;
      S[0] *= A;
      S[1] *= A;
      S[2] *= A;

      double s0 = cache.pVector[7] * S[0] + cache.pVector[8] * S[1]
          + cache.pVector[9] * S[2];
      double s1 = cache.pVector[14] * S[0] + cache.pVector[15] * S[1]
          + cache.pVector[16] * S[2];
      double s2 = cache.pVector[21] * S[0] + cache.pVector[22] * S[1]
          + cache.pVector[23] * S[2];
      double s3 = cache.pVector[28] * S[0] + cache.pVector[29] * S[1]
          + cache.pVector[30] * S[2];
      double s4 = cache.pVector[35] * S[0] + cache.pVector[36] * S[1]
          + cache.pVector[37] * S[2];

      cache.pVector[7] -= (s0 * cache.pVector[3]);
      cache.pVector[8] -= (s0 * cache.pVector[4]);
      cache.pVector[9] -= (s0 * cache.pVector[5]);
      cache.pVector[10] -= (s0 * cache.pVector[42]);
      cache.pVector[11] -= (s0 * cache.pVector[43]);
      cache.pVector[12] -= (s0 * cache.pVector[44]);
      cache.pVector[14] -= (s1 * cache.pVector[3]);
      cache.pVector[15] -= (s1 * cache.pVector[4]);
      cache.pVector[16] -= (s1 * cache.pVector[5]);
      cache.pVector[17] -= (s1 * cache.pVector[42]);
      cache.pVector[18] -= (s1 * cache.pVector[43]);
      cache.pVector[19] -= (s1 * cache.pVector[44]);
      cache.pVector[21] -= (s2 * cache.pVector[3]);
      cache.pVector[22] -= (s2 * cache.pVector[4]);
      cache.pVector[23] -= (s2 * cache.pVector[5]);
      cache.pVector[24] -= (s2 * cache.pVector[42]);
      cache.pVector[25] -= (s2 * cache.pVector[43]);
      cache.pVector[26] -= (s2 * cache.pVector[44]);
      cache.pVector[28] -= (s3 * cache.pVector[3]);
      cache.pVector[29] -= (s3 * cache.pVector[4]);
      cache.pVector[30] -= (s3 * cache.pVector[5]);
      cache.pVector[31] -= (s3 * cache.pVector[42]);
      cache.pVector[32] -= (s3 * cache.pVector[43]);
      cache.pVector[33] -= (s3 * cache.pVector[44]);
      cache.pVector[35] -= (s4 * cache.pVector[3]);
      cache.pVector[36] -= (s4 * cache.pVector[4]);
      cache.pVector[37] -= (s4 * cache.pVector[5]);
      cache.pVector[38] -= (s4 * cache.pVector[42]);
      cache.pVector[39] -= (s4 * cache.pVector[43]);
      cache.pVector[40] -= (s4 * cache.pVector[44]);

      double P3, P4,
          C = cache.pVector[3] * cache.pVector[3]
          + cache.pVector[4] * cache.pVector[4];
      if (C > 1.e-20) {
        C  = 1. / C;
        P3 = cache.pVector[3] * C;
        P4 = cache.pVector[4] * C;
        C  = -sqrt(C);
      } else {
        C  = -1.e10;
        P3 = 1.;
        P4 = 0.;
      }

      // Jacobian production
      //
      cache.jacobian[0]
          = Ax[0] * cache.pVector[7] + Ax[1] * cache.pVector[8];  // dL0/dL0
      cache.jacobian[1]
          = Ax[0] * cache.pVector[14] + Ax[1] * cache.pVector[15];  // dL0/dL1
      cache.jacobian[2]
          = Ax[0] * cache.pVector[21] + Ax[1] * cache.pVector[22];  // dL0/dPhi
      cache.jacobian[3]
          = Ax[0] * cache.pVector[28] + Ax[1] * cache.pVector[29];  // dL0/dThe
      cache.jacobian[4]
          = Ax[0] * cache.pVector[35] + Ax[1] * cache.pVector[36];  // dL0/dCM
      cache.jacobian[5] = Ay[0] * cache.pVector[7] + Ay[1] * cache.pVector[8]
          + Ay[2] * cache.pVector[9];  // dL1/dL0
      cache.jacobian[6] = Ay[0] * cache.pVector[14] + Ay[1] * cache.pVector[15]
          + Ay[2] * cache.pVector[16];  // dL1/dL1
      cache.jacobian[7] = Ay[0] * cache.pVector[21] + Ay[1] * cache.pVector[22]
          + Ay[2] * cache.pVector[23];  // dL1/dPhi
      cache.jacobian[8] = Ay[0] * cache.pVector[28] + Ay[1] * cache.pVector[29]
          + Ay[2] * cache.pVector[30];  // dL1/dThe
      cache.jacobian[9] = Ay[0] * cache.pVector[35] + Ay[1] * cache.pVector[36]
          + Ay[2] * cache.pVector[37];  // dL1/dCM
      cache.jacobian[10]
          = P3 * cache.pVector[11] - P4 * cache.pVector[10];  // dPhi/dL0
      cache.jacobian[11]
          = P3 * cache.pVector[18] - P4 * cache.pVector[17];  // dPhi/dL1
      cache.jacobian[12]
          = P3 * cache.pVector[25] - P4 * cache.pVector[24];  // dPhi/dPhi
      cache.jacobian[13]
          = P3 * cache.pVector[32] - P4 * cache.pVector[31];  // dPhi/dThe
      cache.jacobian[14]
          = P3 * cache.pVector[39] - P4 * cache.pVector[38];  // dPhi/dCM
      cache.jacobian[15] = C * cache.pVector[12];             // dThe/dL0
      cache.jacobian[16] = C * cache.pVector[19];             // dThe/dL1
      cache.jacobian[17] = C * cache.pVector[26];             // dThe/dPhi
      cache.jacobian[18] = C * cache.pVector[33];             // dThe/dThe
      cache.jacobian[19] = C * cache.pVector[40];             // dThe/dCM
      cache.jacobian[20] = 0;                                 // dCM /dL0
      cache.jacobian[21] = 0;                                 // dCM /dL1
      cache.jacobian[22] = 0;                                 // dCM /dPhi
      cache.jacobian[23] = 0;                                 // dCM /dTheta
      cache.jacobian[24] = cache.pVector[41];                 // dCM /dCM
      Eigen::
          Map<Eigen::Matrix<double, NGlobalPars, NGlobalPars, Eigen::RowMajor>>
              J(cache.jacobian);

      cov = std::make_unique<const ActsSymMatrixD<NGlobalPars>>(
          J * (*cache.covariance) * J.transpose());
    }

    return BoundParameters(std::move(cov), gp, mom, charge, s);
  }

  AtlasStepper(BField bField = BField()) : m_bField(std::move(bField)){};

  static double
  distance(const Surface& s, const Vector3D& pos, const Vector3D& dir)
  {
    const Intersection i = s.intersectionEstimate(pos, dir);
    return i.pathLength;
  }

  double
  step(Cache& cache, double& h) const
  {
    bool Jac = cache.useJacobian;

    double* R  = &(cache.pVector[0]);  // Coordinates
    double* A  = &(cache.pVector[3]);  // Directions
    double* sA = &(cache.pVector[42]);
    // Invert mometum/2.
    double Pi = 0.5 / units::Nat2SI<units::MOMENTUM>(1. / cache.pVector[6]);
    //    double dltm = 0.0002 * .03;

    Vector3D f0, f;

    // if new field is required get it
    if (cache.newfield) {
      const Vector3D pos(R[0], R[1], R[2]);
      f0 = m_bField.getField(pos);
    } else {
      f0 = cache.field;
    }

    bool Helix = false;
    // if (std::abs(S) < m_cfg.helixStep) Helix = true;

    while (h != 0.) {
      double S3 = (1. / 3.) * h, S4 = .25 * h, PS2 = Pi * h;

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
        const Vector3D pos(R[0] + A1 * S4, R[1] + B1 * S4, R[2] + C1 * S4);
        f = m_bField.getField(pos);
      } else {
        f = f0;
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
        const Vector3D pos(R[0] + h * A4, R[1] + h * B4, R[2] + h * C4);
        f = m_bField.getField(pos);
      } else {
        f = f0;
      }

      double H2[3] = {f[0] * PS2, f[1] * PS2, f[2] * PS2};
      double A6    = B5 * H2[2] - C5 * H2[1];
      double B6    = C5 * H2[0] - A5 * H2[2];
      double C6    = A5 * H2[1] - B5 * H2[0];

      // Test approximation quality on give step and possible step reduction
      //
      double EST = 2.
          * (std::abs((A1 + A6) - (A3 + A4)) + std::abs((B1 + B6) - (B3 + B4))
             + std::abs((C1 + C6) - (C3 + C4)));
      if (EST > 0.0002) {
        h *= .5;
        //        dltm = 0.;
        continue;
      }

      //      if (EST < dltm) h *= 2.;

      // Parameters calculation
      //
      double A00 = A[0], A11 = A[1], A22 = A[2];

      A[0] = 2. * A3 + (A0 + A5 + A6);
      A[1] = 2. * B3 + (B0 + B5 + B6);
      A[2] = 2. * C3 + (C0 + C5 + C6);

      double D  = (A[0] * A[0] + A[1] * A[1]) + (A[2] * A[2] - 9.);
      double Sl = 2. / h;
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

      cache.field    = f;
      cache.newfield = false;

      // h *= 2;
      if (!Jac) return h;

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
      return h;
    }

    return h;
  }

private:
  BField m_bField;
};

}  // namespace Acts
#endif  // ACTS_ATLAS_STEPPER_HPP
