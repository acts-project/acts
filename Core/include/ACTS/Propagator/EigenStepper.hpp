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
template <typename BField>
class EigenStepper
{
private:
  struct Cache
  {
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
        const auto& transform = par.referenceSurface().transform().matrix();
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

    Vector3D
    position() const
    {
      return pos;
    }
    Vector3D
    direction() const
    {
      return dir;
    }

    Vector3D pos = Vector3D(0, 0, 0);
    Vector3D dir = Vector3D(1, 0, 0);
    double   qop = 1;
    ActsMatrixD<7, 5> jacobian = ActsMatrixD<7, 5>::Zero();
    ActsVectorD<7>           derivative = ActsVectorD<7>::Zero();
    const ActsSymMatrixD<5>* cov        = nullptr;
  };

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

  static ActsMatrixD<3, 3>
  dLocaldGlobal(const Surface& p, const Vector3D& gpos)
  {
    ActsMatrixD<3, 3> j = ActsMatrixD<3, 3>::Zero();
    j.block<1, 3>(0, 0) = p.transform().matrix().block<3, 1>(0, 0).transpose();
    j.block<1, 3>(1, 0) = p.transform().matrix().block<3, 1>(0, 1).transpose();
    j.block<1, 3>(2, 0) = p.transform().matrix().block<3, 1>(0, 2).transpose();

    return j;
  }

public:
  template <typename T, typename S = int>
  using cache_type = Cache;

  template <typename T>
  using step_parameter_type = CurvilinearParameters;

  template <typename T, typename S = int>
  using return_parameter_type = typename s<T, S>::type;

  EigenStepper(BField bField = BField()) : m_bField(std::move(bField)){};

  static CurvilinearParameters
  convert(const Cache& c)
  {
    double                                   charge = c.qop > 0. ? 1. : -1.;
    std::unique_ptr<const ActsSymMatrixD<5>> cov    = nullptr;
    if (c.cov) {
      const double phi   = c.dir.phi();
      const double theta = c.dir.theta();

      ActsMatrixD<5, 7> J = ActsMatrixD<5, 7>::Zero();
      if (std::abs(cos(theta)) < 0.99) {
        J(0, 0) = -sin(phi);
        J(0, 1) = cos(phi);
        J(1, 0) = -cos(phi) * cos(theta);
        J(1, 1) = -sin(phi) * cos(theta);
        J(1, 2) = sin(theta);
      } else {
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

  template <typename S>
  static BoundParameters
  convert(Cache& c, const S& surface)
  {
    double                             charge = c.qop > 0. ? 1. : -1.;
    std::unique_ptr<const ActsSymMatrixD<5>> cov    = nullptr;
    if (c.cov) {
      const double phi   = c.dir.phi();
      const double theta = c.dir.theta();

      ActsMatrixD<5, 7> J = ActsMatrixD<5, 7>::Zero();
      const auto& dLdG = dLocaldGlobal(surface, c.pos);
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

  static double
  distance(const Surface& s, const Vector3D& pos, const Vector3D& dir)
  {
    const Intersection& i = s.intersectionEstimate(pos, dir);
    return i.pathLength;
  }

  double
  step(Cache& c, double& h) const
  {
    const double qop = 1. / units::Nat2SI<units::MOMENTUM>(1. / c.qop);

    // first point
    const Vector3D B_first = m_bField.getField(c.pos);
    const Vector3D& k1 = qop * c.dir.cross(B_first);

    while (h != 0.) {
      const double& h2     = h * h;
      const double& half_h = h / 2;
      // second point
      const Vector3D& pos1 = c.pos + half_h * c.dir + h2 / 8 * k1;
      const Vector3D B_middle = m_bField.getField(pos1);
      const Vector3D& k2 = qop * (c.dir + half_h * k1).cross(B_middle);

      // third point
      const Vector3D& k3 = qop * (c.dir + half_h * k2).cross(B_middle);

      // last point
      const Vector3D& pos2 = c.pos + h * c.dir + h2 / 2 * k3;
      const Vector3D B_last = m_bField.getField(pos2);
      const Vector3D& k4 = qop * (c.dir + h * k3).cross(B_last);

      // local error estimate
      const double EST = h * (k1 - k2 - k3 + k4).template lpNorm<1>();
      if (EST > 0.0002) {
        h *= 0.5;
        continue;
      }

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

      c.pos += h * c.dir + h2 / 6 * (k1 + k2 + k3);
      c.dir += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
      c.dir /= c.dir.norm();
      c.derivative.template head<3>()     = c.dir;
      c.derivative.template segment<3>(3) = k4;

      return h;
    }

    return h;
  }

private:
  BField m_bField;
};

}  // namespace Acts
#endif  // ACTS_EIGEN_STEPPER_HPP
