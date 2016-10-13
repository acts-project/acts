#ifndef ACTS_EIGEN_STEPPER_HPP
#define ACTS_EIGEN_STEPPER_HPP 1

#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Utilities/Definitions.hpp"

namespace Acts {

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
    {
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
    double charge = c.qop > 0. ? 1. : -1.;
    return CurvilinearParameters(nullptr, c.pos, c.dir / fabs(c.qop), charge);
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
      // second point
      const Vector3D& pos1 = c.pos + h / 2 * c.dir + h * h / 8 * k1;
      m_bField(pos1.data(), B_middle.data());
      const Vector3D& k2 = qop * (c.dir + h / 2 * k1).cross(B_middle);

      // third point
      const Vector3D& k3 = qop * (c.dir + h / 2 * k2).cross(B_middle);

      // last point
      const Vector3D& pos2 = c.pos + h * c.dir + h * h / 2 * k3;
      m_bField(pos2.data(), B_last.data());
      const Vector3D& k4 = qop * (c.dir + h * k3).cross(B_last);

      // local error estimate
      const double EST = h * (k1 - k2 - k3 + k4).template lpNorm<1>();
      if (EST > 0.0002) {
        h *= 0.5;
        continue;
      }

      c.pos += h * c.dir + h * h / 6 * (k1 + k2 + k3);
      c.dir += h / 6 * (k1 + 2 * k2 + 2 * k3 + k4);
      c.dir /= c.dir.norm();

      return h;
    }

    return h;
  }

private:
  BField m_bField;
};

}  // namespace Acts
#endif  // ACTS_EIGEN_STEPPER_HPP
