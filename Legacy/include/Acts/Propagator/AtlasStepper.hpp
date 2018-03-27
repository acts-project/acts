// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
#include <cmath>
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/concept/AnyFieldLookup.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Utilities/Definitions.hpp"

// This is based original stepper code from the ATLAS RungeKuttePropagagor
namespace Acts {

template <typename BField>
class AtlasStepper
{

public:
  struct Cache
  {
    // optimisation that init is not called twice
    bool cache_ready = false;
    // configuration
    NavigationDirection nav_dir;
    bool                useJacobian;
    double              step;
    double              maxPathLength;
    bool                mcondition;
    bool                needgradient;
    bool                newfield;
    // internal parameters to be used
    Vector3D field;
    double   pVector[64];
    // result
    double parameters[NGlobalPars] = {0., 0., 0., 0., 0.};
    const ActsSymMatrixD<NGlobalPars>* covariance;
    double                             jacobian[NGlobalPars * NGlobalPars];

    /// Lazily initialized cache for the magnetic field
    /// It caches the current magnetic field cell and stays interpolates within
    /// as long as this is valid. See step() code for details.
    bool                    field_cache_ready = false;
    concept::AnyFieldCell<> field_cache;

    // accummulated path length cache
    double accumulated_path = 0.;

    // adaptive step size of the runge-kutta integration
    double step_size = std::numeric_limits<double>::max();

    // maximal step size of the runge-kutta integration
    double max_step_size = std::numeric_limits<double>::max();

    /// Navigation cache: the start surface
    const Surface* start_surface = nullptr;

    /// Navigation cache: the current surface
    const Surface* current_surface = nullptr;

    /// Navigation cache: the target surface
    const Surface* target_surface = nullptr;

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

    /// Constructor
    template <typename Parameters>
    Cache(const Parameters&   pars,
          NavigationDirection ndir  = forward,
          double              ssize = std::numeric_limits<double>::max())
      : cache_ready(false)
      , nav_dir(ndir)
      , useJacobian(false)
      , step(0.)
      , maxPathLength(0.)
      , mcondition(false)
      , needgradient(false)
      , newfield(true)
      , field(0., 0., 0.)
      , covariance(nullptr)
      , step_size(ssize)
      , max_step_size(ssize)
    {
      update(pars);
    }

    /// The cache update method
    ///
    /// @param [in] pars The new track parameters at start
    ///
    /// @todo check to identify an reuse of start/cache
    template <typename Parameters>
    void
    update(const Parameters& pars)
    {
      // cache is ready - noting to do
      if (cache_ready) return;

      const ActsVectorD<3> pos = pars.position();
      const auto           Vp  = pars.parameters();

      double Sf, Cf, Ce, Se;
      Sf = sin(Vp(2));
      Cf = cos(Vp(2));
      Se = sin(Vp(3));
      Ce = cos(Vp(3));

      pVector[0] = pos(0);
      pVector[1] = pos(1);
      pVector[2] = pos(2);
      pVector[3] = Cf * Se;
      pVector[4] = Sf * Se;
      pVector[5] = Ce;
      pVector[6] = Vp[4];

      // @todo: remove magic numbers - is that the charge ?
      if (std::abs(pVector[6]) < .000000000000001) {
        pVector[6] < 0. ? pVector[6] = -.000000000000001
                        : pVector[6] = .000000000000001;
      }

      // prepare the jacobian if we have a covariance
      if (pars.covariance()) {
        // copy the covariance matrix
        covariance  = new ActsSymMatrixD<NGlobalPars>(*pars.covariance());
        useJacobian = true;
        const auto transform = pars.referenceFrame();

        pVector[7]  = transform(0, eLOC_0);
        pVector[14] = transform(0, eLOC_1);
        pVector[21] = 0.;
        pVector[28] = 0.;
        pVector[35] = 0.;  // dX /

        pVector[8]  = transform(1, eLOC_0);
        pVector[15] = transform(1, eLOC_1);
        pVector[22] = 0.;
        pVector[29] = 0.;
        pVector[36] = 0.;  // dY /

        pVector[9]  = transform(2, eLOC_0);
        pVector[16] = transform(2, eLOC_1);
        pVector[23] = 0.;
        pVector[30] = 0.;
        pVector[37] = 0.;  // dZ /

        pVector[10] = 0.;
        pVector[17] = 0.;
        pVector[24] = -Sf * Se;  // - sin(phi) * cos(theta)
        pVector[31] = Cf * Ce;   // cos(phi) * cos(theta)
        pVector[38] = 0.;        // dAx/

        pVector[11] = 0.;
        pVector[18] = 0.;
        pVector[25] = Cf * Se;  // cos(phi) * sin(theta)
        pVector[32] = Sf * Ce;  // sin(phi) * cos(theta)
        pVector[39] = 0.;       // dAy/

        pVector[12] = 0.;
        pVector[19] = 0.;
        pVector[26] = 0.;
        pVector[33] = -Se;  // - sin(theta)
        pVector[40] = 0.;   // dAz/

        pVector[13] = 0.;
        pVector[20] = 0.;
        pVector[27] = 0.;
        pVector[34] = 0.;
        pVector[41] = 1.;  // dCM/

        pVector[42] = 0.;
        pVector[43] = 0.;
        pVector[44] = 0.;

        // special treatment for surface types
        const auto& surface = pars.referenceSurface();
        // the disc needs polar coordinate adaptations
        if (surface.type() == Surface::Disc) {
          double lCf   = cos(Vp[1]);
          double lSf   = sin(Vp[1]);
          double Ax[3] = {transform(0, 0), transform(1, 0), transform(2, 0)};
          double Ay[3] = {transform(0, 1), transform(1, 1), transform(2, 1)};
          double d0    = lCf * Ax[0] + lSf * Ay[0];
          double d1    = lCf * Ax[1] + lSf * Ay[1];
          double d2    = lCf * Ax[2] + lSf * Ay[2];
          pVector[7]   = d0;
          pVector[8]   = d1;
          pVector[9]   = d2;
          pVector[14]  = Vp[0] * (lCf * Ay[0] - lSf * Ax[0]);
          pVector[15]  = Vp[0] * (lCf * Ay[1] - lSf * Ax[1]);
          pVector[16]  = Vp[0] * (lCf * Ay[2] - lSf * Ax[2]);
        }
        // the line needs components that relate direction change
        // with global frame change
        if (surface.type() == Surface::Perigee
            || surface.type() == Surface::Straw) {

          // sticking to the nomenclature of the original RkPropagator
          // - axis pointing along the drift/transverse direction
          double B[3] = {transform(0, 0), transform(1, 0), transform(2, 0)};
          // - axis along the straw
          double A[3] = {transform(0, 1), transform(1, 1), transform(2, 1)};
          // - normal vector of the reference frame
          double C[3] = {transform(0, 2), transform(1, 2), transform(2, 2)};

          // projection of direction onto normal vector of reference frame
          double PC = pVector[3] * C[0] + pVector[4] * C[1] + pVector[5] * C[2];
          double Bn = 1. / PC;

          double Bx2 = -A[2] * pVector[25];
          double Bx3 = A[1] * pVector[33] - A[2] * pVector[32];

          double By2 = A[2] * pVector[24];
          double By3 = A[2] * pVector[31] - A[0] * pVector[33];

          double Bz2 = A[0] * pVector[25] - A[1] * pVector[24];
          double Bz3 = A[0] * pVector[32] - A[1] * pVector[31];

          double B2 = B[0] * Bx2 + B[1] * By2 + B[2] * Bz2;
          double B3 = B[0] * Bx3 + B[1] * By3 + B[2] * Bz3;

          Bx2 = (Bx2 - B[0] * B2) * Bn;
          Bx3 = (Bx3 - B[0] * B3) * Bn;
          By2 = (By2 - B[1] * B2) * Bn;
          By3 = (By3 - B[1] * B3) * Bn;
          Bz2 = (Bz2 - B[2] * B2) * Bn;
          Bz3 = (Bz3 - B[2] * B3) * Bn;

          //  /dPhi      |     /dThe       |
          pVector[21] = Bx2 * Vp[0];
          pVector[28] = Bx3 * Vp[0];  // dX/
          pVector[22] = By2 * Vp[0];
          pVector[29] = By3 * Vp[0];  // dY/
          pVector[23] = Bz2 * Vp[0];
          pVector[30] = Bz3 * Vp[0];  // dZ/
        }
      }
      // now declare the cache as ready
      cache_ready = true;
    }
  };

  template <typename T, typename S = int>
  using cache_type = Cache;

  template <typename T>
  using step_parameter_type = CurvilinearParameters;

  // This struct is a meta-function which normally maps to BoundParameters...
  template <typename T, typename S>
  struct s
  {
    typedef BoundParameters type;
  };

  // Unless S is int, then it maps to CurvilinearParameters ...
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
    // the convert method invalidates the cache (in case it's reused)
    cache.cache_ready = false;
    //
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

      cache.jacobian[10] = P3 * P[11] - P4 * P[10];  // dPhi/dL0
      cache.jacobian[11] = P3 * P[18] - P4 * P[17];  // dPhi/dL1
      cache.jacobian[12] = P3 * P[25] - P4 * P[24];  // dPhi/dPhi
      cache.jacobian[13] = P3 * P[32] - P4 * P[31];  // dPhi/dThe
      cache.jacobian[14] = P3 * P[39] - P4 * P[38];  // dPhi/dCM

      cache.jacobian[15] = C * P[12];  // dThe/dL0
      cache.jacobian[16] = C * P[19];  // dThe/dL1
      cache.jacobian[17] = C * P[26];  // dThe/dPhi
      cache.jacobian[18] = C * P[33];  // dThe/dThe
      cache.jacobian[19] = C * P[40];  // dThe/dCM

      cache.jacobian[20] = 0.;     // dCM /dL0
      cache.jacobian[21] = 0.;     // dCM /dL1
      cache.jacobian[22] = 0.;     // dCM /dPhi
      cache.jacobian[23] = 0.;     // dCM /dTheta
      cache.jacobian[24] = P[41];  // dCM /dCM

      Eigen::
          Map<Eigen::Matrix<double, NGlobalPars, NGlobalPars, Eigen::RowMajor>>
              J(cache.jacobian);

      cov = std::make_unique<const ActsSymMatrixD<NGlobalPars>>(
          J * (*cache.covariance) * J.transpose());
    }

    return CurvilinearParameters(std::move(cov), gp, mom, charge);
  }

  /// convert method into bound parameters
  /// @param cache is the propagation cache to be converted
  /// @param s the target surface
  static BoundParameters
  convert(Cache& cache, const Surface& s)
  {

    // the convert method invalidates the cache (in case it's reused)
    cache.cache_ready = false;

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

      const auto fFrame = s.referenceFrame(gp, mom);

      double Ax[3] = {fFrame(0, 0), fFrame(1, 0), fFrame(2, 0)};
      double Ay[3] = {fFrame(0, 1), fFrame(1, 1), fFrame(2, 1)};
      double S[3]  = {fFrame(0, 2), fFrame(1, 2), fFrame(2, 2)};

      // this is the projection of direction onto the local normal vector
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

      // in case of line-type surfaces - we need to take into account that
      // the reference frame changes with variations of all local
      // parameters
      if (s.type() == Surface::Straw || s.type() == Surface::Perigee) {
        // vector from position to center
        double x = cache.pVector[0] - s.center().x();
        double y = cache.pVector[1] - s.center().y();
        double z = cache.pVector[2] - s.center().z();

        // this is the projection of the direction onto the local y axis
        double d = cache.pVector[3] * Ay[0] + cache.pVector[4] * Ay[1]
            + cache.pVector[5] * Ay[2];

        // this is cos(beta)
        double a       = (1. - d) * (1. + d);
        if (a != 0.) a = 1. / a;  // i.e. 1./(1-d^2)

        // that's the modified norm vector
        double X = d * Ay[0] - cache.pVector[3];  //
        double Y = d * Ay[1] - cache.pVector[4];  //
        double Z = d * Ay[2] - cache.pVector[5];  //

        // d0 to d1
        double d0 = cache.pVector[10] * Ay[0] + cache.pVector[11] * Ay[1]
            + cache.pVector[12] * Ay[2];
        double d1 = cache.pVector[17] * Ay[0] + cache.pVector[18] * Ay[1]
            + cache.pVector[19] * Ay[2];
        double d2 = cache.pVector[24] * Ay[0] + cache.pVector[25] * Ay[1]
            + cache.pVector[26] * Ay[2];
        double d3 = cache.pVector[31] * Ay[0] + cache.pVector[32] * Ay[1]
            + cache.pVector[33] * Ay[2];
        double d4 = cache.pVector[38] * Ay[0] + cache.pVector[39] * Ay[1]
            + cache.pVector[40] * Ay[2];

        s0 = (((cache.pVector[7] * X + cache.pVector[8] * Y
                + cache.pVector[9] * Z)
               + x * (d0 * Ay[0] - cache.pVector[10]))
              + (y * (d0 * Ay[1] - cache.pVector[11])
                 + z * (d0 * Ay[2] - cache.pVector[12])))
            * (-a);

        s1 = (((cache.pVector[14] * X + cache.pVector[15] * Y
                + cache.pVector[16] * Z)
               + x * (d1 * Ay[0] - cache.pVector[17]))
              + (y * (d1 * Ay[1] - cache.pVector[18])
                 + z * (d1 * Ay[2] - cache.pVector[19])))
            * (-a);
        s2 = (((cache.pVector[21] * X + cache.pVector[22] * Y
                + cache.pVector[23] * Z)
               + x * (d2 * Ay[0] - cache.pVector[24]))
              + (y * (d2 * Ay[1] - cache.pVector[25])
                 + z * (d2 * Ay[2] - cache.pVector[26])))
            * (-a);
        s3 = (((cache.pVector[28] * X + cache.pVector[29] * Y
                + cache.pVector[30] * Z)
               + x * (d3 * Ay[0] - cache.pVector[31]))
              + (y * (d3 * Ay[1] - cache.pVector[32])
                 + z * (d3 * Ay[2] - cache.pVector[33])))
            * (-a);
        s4 = (((cache.pVector[35] * X + cache.pVector[36] * Y
                + cache.pVector[37] * Z)
               + x * (d4 * Ay[0] - cache.pVector[38]))
              + (y * (d4 * Ay[1] - cache.pVector[39])
                 + z * (d4 * Ay[2] - cache.pVector[40])))
            * (-a);
      }

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

      double MA[3] = {Ax[0], Ax[1], Ax[2]};
      double MB[3] = {Ay[0], Ay[1], Ay[2]};
      // Jacobian production of transport and to_local
      if (s.type() == Surface::Disc) {
        // the vector from the disc surface to the p
        const auto& sfc  = s.center();
        double      d[3] = {cache.pVector[0] - sfc(0),
                       cache.pVector[1] - sfc(1),
                       cache.pVector[2] - sfc(2)};
        // this needs the transformation to polar coordinates
        double RC = d[0] * Ax[0] + d[1] * Ax[1] + d[2] * Ax[2];
        double RS = d[0] * Ay[0] + d[1] * Ay[1] + d[2] * Ay[2];
        double R2 = RC * RC + RS * RS;

        // inverse radius
        double Ri = 1. / sqrt(R2);
        MA[0]     = (RC * Ax[0] + RS * Ay[0]) * Ri;
        MA[1]     = (RC * Ax[1] + RS * Ay[1]) * Ri;
        MA[2]     = (RC * Ax[2] + RS * Ay[2]) * Ri;
        MB[0] = (RC * Ay[0] - RS * Ax[0]) * (Ri = 1. / R2);
        MB[1] = (RC * Ay[1] - RS * Ax[1]) * Ri;
        MB[2] = (RC * Ay[2] - RS * Ax[2]) * Ri;
      }

      cache.jacobian[0] = MA[0] * cache.pVector[7] + MA[1] * cache.pVector[8]
          + MA[2] * cache.pVector[9];  // dL0/dL0
      cache.jacobian[1] = MA[0] * cache.pVector[14] + MA[1] * cache.pVector[15]
          + MA[2] * cache.pVector[16];  // dL0/dL1
      cache.jacobian[2] = MA[0] * cache.pVector[21] + MA[1] * cache.pVector[22]
          + MA[2] * cache.pVector[23];  // dL0/dPhi
      cache.jacobian[3] = MA[0] * cache.pVector[28] + MA[1] * cache.pVector[29]
          + MA[2] * cache.pVector[30];  // dL0/dThe
      cache.jacobian[4] = MA[0] * cache.pVector[35] + MA[1] * cache.pVector[36]
          + MA[2] * cache.pVector[37];  // dL0/dCM

      cache.jacobian[5] = MB[0] * cache.pVector[7] + MB[1] * cache.pVector[8]
          + MB[2] * cache.pVector[9];  // dL1/dL0
      cache.jacobian[6] = MB[0] * cache.pVector[14] + MB[1] * cache.pVector[15]
          + MB[2] * cache.pVector[16];  // dL1/dL1
      cache.jacobian[7] = MB[0] * cache.pVector[21] + MB[1] * cache.pVector[22]
          + MB[2] * cache.pVector[23];  // dL1/dPhi
      cache.jacobian[8] = MB[0] * cache.pVector[28] + MB[1] * cache.pVector[29]
          + MB[2] * cache.pVector[30];  // dL1/dThe
      cache.jacobian[9] = MB[0] * cache.pVector[35] + MB[1] * cache.pVector[36]
          + MB[2] * cache.pVector[37];  // dL1/dCM

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
      cache.jacobian[20] = 0.;                                // dCM /dL0
      cache.jacobian[21] = 0.;                                // dCM /dL1
      cache.jacobian[22] = 0.;                                // dCM /dPhi
      cache.jacobian[23] = 0.;                                // dCM /dTheta
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

  /// Get the field for the stepping
  /// It checks first if the access is still within the Cell,
  /// and updates the cell if necessary, then it takes the field
  /// from the cell
  /// @param [in,out] cache is the propagation cache associated with the track
  ///                 the magnetic field cell is used (and potentially updated)
  /// @param [in] pos is the field position
  Vector3D
  getField(Cache& cache, const Vector3D& pos) const
  {
    if (!cache.field_cache_ready || !cache.field_cache.isInside(pos)) {
      cache.field_cache_ready = true;
      cache.field_cache       = m_bField.getFieldCell(pos);
    }
    // get the field from the cell
    cache.field = cache.field_cache.getField(pos);
    return cache.field;
  }

  double
  step(Cache& cache) const
  {

    // we use h for keeping the nominclature with the original atlas code
    double& h   = cache.step_size;
    bool    Jac = cache.useJacobian;

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
      f0 = getField(cache, pos);
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
        f = getField(cache, pos);
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
        f = getField(cache, pos);
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

      if (Jac) {
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
      }
      cache.accumulated_path += h;
      return h;
    }
    cache.accumulated_path += h;
    return h;
  }

private:
  BField m_bField;
};

}  // namespace Acts