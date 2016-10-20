#ifndef ACTS_EIGEN_STEPPER_HPP
#define ACTS_EIGEN_STEPPER_HPP 1

#include "ACTS/EventData/TrackParameters.hpp"
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
    explicit Cache(const CurvilinearParameters& par)
      : pos(par.position())
      , dir(par.momentum().normalized())
      , qop(par.charge() / par.momentum().norm())
      , cov(par.covariance())
    {
      const double phi   = dir.phi();
      const double theta = dir.theta();

      if (cov) {
        const auto& transform = par.associatedSurface().transform().matrix();
        jacobian(0, eLOC_1) = transform(0, eLOC_1);
        jacobian(0, eLOC_2) = transform(0, eLOC_2);
        jacobian(1, eLOC_1) = transform(1, eLOC_1);
        jacobian(1, eLOC_2) = transform(1, eLOC_2);
        jacobian(2, eLOC_1) = transform(2, eLOC_1);
        jacobian(2, eLOC_2) = transform(2, eLOC_2);
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
    ActsMatrixD<7, 5> jacobian = ActsMatrixD<7, 5>::Identity();
    const ActsSymMatrixD<5>* cov = nullptr;
  };

public:
  template <typename T>
  using cache_type = Cache;

  template <typename T, typename S = int>
  using return_parameter_type = CurvilinearParameters;

  EigenStepper(BField bField = BField()) : m_bField(std::move(bField)){};

  CurvilinearParameters
  convert(const Cache& c) const
  {
    double                             charge = c.qop > 0. ? 1. : -1.;
    std::unique_ptr<ActsSymMatrixD<5>> cov    = nullptr;
    if (c.cov) {
      const double phi   = c.dir.phi();
      const double theta = c.dir.theta();

      ActsMatrixD<5, 7> J = ActsMatrixD<5, 7>::Zero();
      J(0, 0)             = -sin(phi);
      J(0, 1)             = cos(phi);
      J(1, 0)             = -cos(phi) * cos(theta);
      J(1, 1)             = -sin(phi) * cos(theta);
      J(1, 2)             = sin(theta);
      J(2, 3)             = -sin(phi) / sin(theta);
      J(2, 4)             = cos(phi) / sin(theta);
      J(3, 5)             = -1. / sin(theta);
      J(4, 6)             = 1;

      auto jac = J * c.jacobian;

      cov = std::make_unique<ActsSymMatrixD<5>>(jac * (*c.cov)
                                                * jac.transpose());
    }

    return CurvilinearParameters(
        std::move(cov), c.pos, c.dir / fabs(c.qop), charge);
  }

  double
  step(Cache& c, double& h) const
  {
    const double qop = 1. / units::Nat2SI<units::MOMENTUM>(1. / c.qop);

    Vector3D B_first(0, 0, 0);
    Vector3D B_middle(0, 0, 0);
    Vector3D B_last(0, 0, 0);
    m_bField(c.pos.data(), B_first.data());

    // first point
    const Vector3D& k1 = qop * c.dir.cross(B_first);

    while (h != 0.) {
      const double& h2     = h * h;
      const double& half_h = h / 2;
      // second point
      const Vector3D& pos1 = c.pos + half_h * c.dir + h2 / 8 * k1;
      m_bField(pos1.data(), B_middle.data());
      const Vector3D& k2 = qop * (c.dir + half_h * k1).cross(B_middle);

      // third point
      const Vector3D& k3 = qop * (c.dir + half_h * k2).cross(B_middle);

      // last point
      const Vector3D& pos2 = c.pos + h * c.dir + h2 / 2 * k3;
      m_bField(pos2.data(), B_last.data());
      const Vector3D& k4 = qop * (c.dir + h * k3).cross(B_last);

      // local error estimate
      const double EST = h * (k1 - k2 - k3 + k4).template lpNorm<1>();
      if (EST > 0.0002) {
        h *= 0.5;
        continue;
      }

      c.pos += h * c.dir + h2 / 6 * (k1 + k2 + k3);
      c.dir += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
      c.dir /= c.dir.norm();

      if (c.cov) {
        ActsMatrixD<7, 7> D = ActsMatrixD<7, 7>::Zero();
        D(6, 6)             = 1;

        auto dFdx = D.block<3, 3>(0, 0);
        auto dFdT = D.block<3, 3>(0, 3);
        auto dGdx = D.block<3, 3>(3, 0);
        auto dGdT = D.block<3, 3>(3, 3);

        ActsMatrixD<3, 3> dk1dT = ActsMatrixD<3, 3>::Zero();
        ActsMatrixD<3, 3> dk2dT = ActsMatrixD<3, 3>::Identity();
        ActsMatrixD<3, 3> dk3dT = ActsMatrixD<3, 3>::Identity();
        ActsMatrixD<3, 3> dk4dT = ActsMatrixD<3, 3>::Identity();

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

        dFdT += h / 6 * (dk1dT + dk2dT + dk3dT);
        dFdT *= h;

        dGdx.setZero();

        dGdT.setIdentity();
        dGdT += h / 6 * (dk1dT + 2 * (dk2dT + dk3dT) + dk4dT);

        c.jacobian = D * c.jacobian;
      }

      return h;
    }

    return h;
  }

private:
  BField m_bField;
};

}  // namespace Acts
#endif  // ACTS_EIGEN_STEPPER_HPP
