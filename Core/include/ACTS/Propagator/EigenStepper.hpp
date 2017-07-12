// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_EIGEN_STEPPER_HPP
#define ACTS_EIGEN_STEPPER_HPP 1

#include <cmath>
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Surfaces/Surface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
#include "ACTS/Utilities/Units.hpp"

namespace Acts {

ActsMatrixD<3, 3>
cross(const ActsMatrixD<3, 3>& m, const Vector3D& v)
{
  ActsMatrixD<3, 3> r;
  r.col(0) = m.col(0).cross(v);
  r.col(1) = m.col(1).cross(v);
  r.col(2) = m.col(2).cross(v);

  return r;
}

/// Runge-Kutta-Nystroem stepper for the following ODE:
///
/// r = (x,y,z)    ... global position
/// T = (Ax,Ay,Az) ... momentum direction (normalized)
///
/// dr/ds = T
/// dT/ds = q/p * (T x B)
///
/// with s being the arc length of the track, q the charge of the particle,
/// p its momentum and B the magnetic field
///
template <typename BField>
class EigenStepper
{
private:

  /// Internal cache for track parameter propagation
  struct Cache
  {
    /// Constructor from the initial track parameters
    template <typename T>
    explicit Cache(const T& par)
      : pos(par.position())
      , dir(par.momentum().normalized())
      , qop(par.charge() / par.momentum().norm())
      , cov(par.covariance())
    {
      const double phi   = dir.phi();
      const double theta = dir.theta();

      if (cov) {
        const auto transform = par.referenceSurface().transform().matrix();
        jacobian(0, eLOC_0) = transform(0, eLOC_0);
        jacobian(0, eLOC_1) = transform(0, eLOC_1);
        jacobian(1, eLOC_0) = transform(1, eLOC_0);
        jacobian(1, eLOC_1) = transform(1, eLOC_1);
        jacobian(2, eLOC_0) = transform(2, eLOC_0);
        jacobian(2, eLOC_1) = transform(2, eLOC_1);
        jacobian(3, ePHI)   = -sin(theta) * sin(phi);
        jacobian(3, eTHETA) = cos(theta) * cos(phi);
        jacobian(4, ePHI)   = sin(theta) * cos(phi);
        jacobian(4, eTHETA) = cos(theta) * sin(phi);
        jacobian(5, eTHETA) = -sin(theta);
        jacobian(6, eQOP)   = 1;
      }
    }

    /// Global particle position accessor
    Vector3D
    position() const
    {
      return pos;
    }

    /// Momentum direction accessor
    Vector3D
    direction() const
    {
      return dir;
    }

    /// Global particle position
    Vector3D pos = Vector3D(0, 0, 0);

    /// Momentum direction (normalized)
    Vector3D dir = Vector3D(1, 0, 0);

    /// Charge-momentum ratio, in natural units
    double qop = 1;

    /// Jacobian used to transport the covariance matrix
    ActsMatrixD<7, 5> jacobian = ActsMatrixD<7, 5>::Zero();

    /// ???
    ActsVectorD<7> derivative = ActsVectorD<7>::Zero();

    /// Covariance matrix assocated with the initial error on track parameters
    const ActsSymMatrixD<5>* cov = nullptr;
  };


  // This struct is a meta-function which normally maps to BoundParameters...
  template <typename T, typename S>
  struct s
  {
    typedef BoundParameters type;
  };

  // ...unless type S is int, in which case it maps to Curvilinear parameters
  template <typename T>
  struct s<T, int>
  {
    typedef CurvilinearParameters type;
  };


  /// Rotation matrix going from global coordinates to local surface coordinates
  static ActsMatrixD<3, 3>
  dLocaldGlobal(const Surface& p, const Vector3D& gpos)
  {
    // A surface's associated transform is a 4x4 affine transformation matrix,
    // whose top-left corner is the 3x3 linear local-to-global rotation matrix,
    // and whose right column is a translation vector, with a trailing 1.
    //
    // So to get the global-to-local rotation matrix, we only need to take the
    // transpose of the top-left corner of the surface's associated transform.
    //
    ActsMatrixD<3, 3> j = ActsMatrixD<3, 3>::Zero();
    j.block<1, 3>(0, 0) = p.transform().matrix().block<3, 1>(0, 0).transpose();
    j.block<1, 3>(1, 0) = p.transform().matrix().block<3, 1>(0, 1).transpose();
    j.block<1, 3>(2, 0) = p.transform().matrix().block<3, 1>(0, 2).transpose();

    return j;
  }


public:

  /// Always use the same propagation cache type, independently of the initial
  /// track parameter type and of the target surface
  template <typename T, typename S = int>
  using cache_type = Cache;

  /// Intermediate track parameters are always in curvilinear parametrization
  template <typename T>
  using step_parameter_type = CurvilinearParameters;

  /// Return parameter types depend on the propagation mode: when propagating to
  /// a surface we return BoundParameters, otherwise CurvilinearParameters
  template <typename T, typename S = int>
  using return_parameter_type = typename s<T, S>::type;


  /// Constructor requires knowledge of the detector's magnetic field
  EigenStepper(BField bField = BField()) : m_bField(std::move(bField)){};


  /// Convert the propagation cache to curvilinear parameters
  static CurvilinearParameters
  convert(const Cache& c)
  {
    double charge = c.qop > 0. ? 1. : -1.;
    std::unique_ptr<const ActsSymMatrixD<5>> cov = nullptr;

    // Perform error propagation if an initial covariance matrix was provided
    if (c.cov) {
      const double phi   = c.dir.phi();
      const double theta = c.dir.theta();

      ActsMatrixD<5, 7> J = ActsMatrixD<5, 7>::Zero();
      if (std::abs(cos(theta)) < 0.99) {
        // We normally operate in curvilinear coordinates defined as follows
        J(0, 0) = -sin(phi);
        J(0, 1) = cos(phi);
        J(1, 0) = -cos(phi) * cos(theta);
        J(1, 1) = -sin(phi) * cos(theta);
        J(1, 2) = sin(theta);
      } else {
        // Under grazing incidence to z, the above coordinate system definition
        // becomes numerically unstable, and we need to switch to another one
        const double c
            = sqrt(pow(cos(theta), 2) + pow(sin(phi) * sin(theta), 2));
        J(0, 1) = -cos(theta) / c;
        J(0, 2) = sin(phi) * sin(theta) / c;
        J(1, 0) = c;
        J(1, 1) = -cos(phi) * sin(phi) * pow(sin(theta), 2) / c;
        J(1, 2) = -cos(phi) * sin(theta) * cos(theta) / c;
      }
      J(2, 3) = -sin(phi) / sin(theta);
      J(2, 4) = cos(phi) / sin(theta);
      J(3, 5) = -1. / sin(theta);
      J(4, 6) = 1;

      ActsRowVectorD<3> norm_vec(
          cos(phi) * sin(theta), sin(phi) * sin(theta), cos(theta));

      norm_vec /= (norm_vec * c.dir);
      ActsRowVectorD<5> scale_factors
          = (c.jacobian.template block<3, 5>(0, 0).array().colwise()
             * norm_vec.transpose().array())
                .colwise()
                .sum();
      ActsMatrixD<7, 5> tmp;
      tmp.col(0) = c.derivative;
      tmp.col(1) = c.derivative;
      tmp.col(2) = c.derivative;
      tmp.col(3) = c.derivative;
      tmp.col(4) = c.derivative;
      tmp *= scale_factors.asDiagonal();
      auto jacobian = c.jacobian;
      jacobian -= tmp;
      auto jac = J * jacobian;

      cov = std::make_unique<const ActsSymMatrixD<5>>(jac * (*c.cov)
                                                      * jac.transpose());
    }

    return CurvilinearParameters(
        std::move(cov), c.pos, c.dir / std::abs(c.qop), charge);
  }


  /// Convert the propagation cache to track parameters at a certain surface
  template <typename S>
  static BoundParameters
  convert(Cache& c, const S& surface)
  {
    double charge = c.qop > 0. ? 1. : -1.;
    std::unique_ptr<const ActsSymMatrixD<5>> cov = nullptr;

    // Perform error propagation if an initial covariance matrix was provided
    if (c.cov) {
      const double phi   = c.dir.phi();
      const double theta = c.dir.theta();

      ActsMatrixD<5, 7> J = ActsMatrixD<5, 7>::Zero();
      const auto dLdG = dLocaldGlobal(surface, c.pos);
      J.block<2, 3>(0, 0) = dLdG.template block<2, 3>(0, 0);
      J(2, 3) = -sin(phi) / sin(theta);
      J(2, 4) = cos(phi) / sin(theta);
      J(3, 5) = -1. / sin(theta);
      J(4, 6) = 1;

      ActsRowVectorD<3> norm_vec = dLdG.template block<1, 3>(2, 0);
      norm_vec /= (norm_vec * c.dir);
      ActsRowVectorD<5> scale_factors
          = (c.jacobian.template block<3, 5>(0, 0).array().colwise()
             * norm_vec.transpose().array())
                .colwise()
                .sum();
      ActsMatrixD<7, 5> tmp;
      tmp.col(0) = c.derivative;
      tmp.col(1) = c.derivative;
      tmp.col(2) = c.derivative;
      tmp.col(3) = c.derivative;
      tmp.col(4) = c.derivative;
      tmp *= scale_factors.asDiagonal();
      c.jacobian -= tmp;
      auto jac = J * c.jacobian;

      cov = std::make_unique<const ActsSymMatrixD<5>>(jac * (*c.cov)
                                                      * jac.transpose());
    }

    return BoundParameters(
        std::move(cov), c.pos, c.dir / std::abs(c.qop), charge, surface);
  }


  /// Estimate the (signed) distance to a certain surface
  static double
  distance(const Surface& s, const Vector3D& pos, const Vector3D& dir)
  {
    const Intersection i = s.intersectionEstimate(pos, dir);
    return i.pathLength;
  }


  /// Perform a Runge-Kutta track parameter propagation step
  ///
  /// @param[in,out] c is the propagation cache associated with the track
  ///                  parameters that are being propagated.
  ///
  /// @param[in,out] h is the desired step size. It can be negative during
  ///                  backwards track propagation, and since we're using an
  ///                  adaptive integrator, it can be modified by the stepper.
  ///
  double
  step(Cache& c, double& h) const
  {
    // Charge-momentum ratio, in SI units
    const double qop = 1. / units::Nat2SI<units::MOMENTUM>(1. / c.qop);

    // Runge-Kutta integrator state
    double h2, half_h;
    Vector3D B_middle, B_last, k2, k3, k4;

    // First Runge-Kutta point (at current position)
    const Vector3D B_first = m_bField.getField(c.pos);
    const Vector3D k1 = qop * c.dir.cross(B_first);

    // Select the appropriate Runge-Kutta step size
    do {
      // Cache the square and half of the step size
      h2     = h * h;
      half_h = h / 2;

      // Second Runge-Kutta point
      const Vector3D pos1 = c.pos + half_h * c.dir + h2 / 8 * k1;
      B_middle = m_bField.getField(pos1);
      k2 = qop * (c.dir + half_h * k1).cross(B_middle);

      // Third Runge-Kutta point
      k3 = qop * (c.dir + half_h * k2).cross(B_middle);

      // Last Runge-Kutta point
      const Vector3D pos2 = c.pos + h * c.dir + h2 / 2 * k3;
      B_last = m_bField.getField(pos2);
      k4 = qop * (c.dir + h * k3).cross(B_last);

      // Estimate of the local integration error
      const double EST = h * (k1 - k2 - k3 + k4).template lpNorm<1>();

      // Reduce the step size and repeat until the estimated error is low enough
      if (EST > 0.0002) {
        h *= 0.5;
        continue;
      } else {
        break;
      }
    } while(true);

    // When doing error propagation, update the associated Jacobian matrix
    if (c.cov) {
      ActsMatrixD<7, 7> D = ActsMatrixD<7, 7>::Zero();
      D(6, 6)             = 1;

      const double conv = units::SI2Nat<units::MOMENTUM>(1);

      auto dFdx = D.block<3, 3>(0, 0);
      auto dFdT = D.block<3, 3>(0, 3);
      auto dFdL = D.block<3, 1>(0, 6);
      auto dGdx = D.block<3, 3>(3, 0);
      auto dGdT = D.block<3, 3>(3, 3);
      auto dGdL = D.block<3, 1>(3, 6);

      ActsMatrixD<3, 3> dk1dT = ActsMatrixD<3, 3>::Zero();
      ActsMatrixD<3, 3> dk2dT = ActsMatrixD<3, 3>::Identity();
      ActsMatrixD<3, 3> dk3dT = ActsMatrixD<3, 3>::Identity();
      ActsMatrixD<3, 3> dk4dT = ActsMatrixD<3, 3>::Identity();

      ActsVectorD<3> dk1dL = ActsVectorD<3>::Zero();
      ActsVectorD<3> dk2dL = ActsVectorD<3>::Zero();
      ActsVectorD<3> dk3dL = ActsVectorD<3>::Zero();
      ActsVectorD<3> dk4dL = ActsVectorD<3>::Zero();

      dk1dL = c.dir.cross(B_first);
      dk2dL = (c.dir + half_h * k1).cross(B_middle)
                  + qop * half_h * dk1dL.cross(B_middle);
      dk3dL = (c.dir + half_h * k2).cross(B_middle)
                  + qop * half_h * dk2dL.cross(B_middle);
      dk4dL = (c.dir + h * k3).cross(B_last) + qop * h * dk3dL.cross(B_last);

      dk1dT(0, 1) = B_first.z();
      dk1dT(0, 2) = -B_first.y();
      dk1dT(1, 0) = -B_first.z();
      dk1dT(1, 2) = B_first.x();
      dk1dT(2, 0) = B_first.y();
      dk1dT(2, 1) = -B_first.x();
      dk1dT *= qop;

      dk2dT += h / 2 * dk1dT;
      dk2dT = qop * cross(dk2dT, B_middle);

      dk3dT += h / 2 * dk2dT;
      dk3dT = qop * cross(dk3dT, B_middle);

      dk4dT += h * dk3dT;
      dk4dT = qop * cross(dk4dT, B_last);

      dFdx.setIdentity();

      dFdT.setIdentity();
      dFdT += h / 6 * (dk1dT + dk2dT + dk3dT);
      dFdT *= h;

      dFdL = conv * h2 / 6 * (dk1dL + dk2dL + dk3dL);

      dGdx.setZero();

      dGdT.setIdentity();
      dGdT += h / 6 * (dk1dT + 2 * (dk2dT + dk3dT) + dk4dT);

      dGdL = conv * h / 6 * (dk1dL + 2 * (dk2dL + dk3dL) + dk4dL);

      c.jacobian = D * c.jacobian;
    }

    // Update the track parameters according to the equations of motion
    c.pos += h * c.dir + h2 / 6 * (k1 + k2 + k3);
    c.dir += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
    c.dir /= c.dir.norm();
    c.derivative.template head<3>()     = c.dir;
    c.derivative.template segment<3>(3) = k4;

    // Return the updated step size
    return h;
  }


private:

  /// Magnetic field inside of the detector
  BField m_bField;

};

}  // namespace Acts

#endif  // ACTS_EIGEN_STEPPER_HPP
