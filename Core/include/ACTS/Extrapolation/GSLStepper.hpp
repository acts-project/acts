#ifndef ACTS_GSL_STEPPER_HPP
#define ACTS_GSL_STEPPER_HPP 1

#include <array>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <memory>
#include "ACTS/Utilities/Units.hpp"

namespace Acts {

struct DefaultParConverter
{
  template <typename T>
  using internal_parameter_type = T;

  template <typename T>
  using return_parameter_type = T;

  template <typename T>
  static T
  convert(const T& input)
  {
    return input;
  }
};

template <typename BField, typename ParConverter = DefaultParConverter>
class GSLStepper
{
public:
  template <typename T>
  using internal_parameter_type =
      typename ParConverter::template internal_parameter_type<T>;

  template <typename T>
  using return_parameter_type =
      typename ParConverter::template return_parameter_type<T>;

  template <typename T>
  static auto
  convert(const T& input)
  {
    return ParConverter::convert(input);
  }

  GSLStepper(BField&& bField = BField())
    : m_bField(std::move(bField))
    , m_params{0, 0, 0, 1}
    , m_ODE_system({GSLStepper<BField>::ODE_system, 0, 6, m_params.data()})
    , m_GSL_driver(nullptr)
  {
    m_GSL_driver.reset(gsl_odeiv2_driver_alloc_y_new(
        &m_ODE_system, gsl_odeiv2_step_rkck, 0.1 * units::_cm, 1e-6, 0.0));
    gsl_odeiv2_driver_set_nmax(m_GSL_driver.get(), 1);
  }

  GSLStepper(GSLStepper<BField>&& rhs)
    : m_bField(std::move(rhs.m_bField))
    , m_params(std::move(rhs.m_params))
    , m_ODE_system(std::move(rhs.m_ODE_system))
    , m_GSL_driver(std::move(rhs.m_GSL_driver))
  {
    m_ODE_system.params = m_params.data();
    m_GSL_driver->sys   = &m_ODE_system;
  }

  template <typename TrackParameters>
  TrackParameters
  doStep(const TrackParameters& in, double stepMax = 1 * units::_cm)
  {
    Vector3D pos    = in.position();
    Vector3D mom    = in.momentum();
    double   p      = mom.norm();
    double   charge = in.charge();
    double   input[6]
        = {pos(0), pos(1), pos(2), mom(0) / p, mom(1) / p, mom(2) / p};

    // set B-Field and q/p
    m_bField(input, m_params.data());
    // conversion to SI units: p --> (p/c)/J / (kg * m^2 / s^2)
    m_params[3] = charge / units::Nat2SI<units::MOMENTUM>(p);

    performStep(input, stepMax);

    return CurvilinearParameters(nullptr,
                                 Vector3D(input[0], input[1], input[2]),
                                 p * Vector3D(input[3], input[4], input[5]),
                                 in.charge());
    ;
  }

private:
  BField m_bField;
  std::array<double, 4> m_params;
  gsl_odeiv2_system                  m_ODE_system;
  std::unique_ptr<gsl_odeiv2_driver> m_GSL_driver;

  /// input = [x,y,z,Tx,Ty,Tz,q/p]
  bool
  performStep(double input[6], double stepMax) const
  {
    double t = 0;
    gsl_odeiv2_driver_apply(m_GSL_driver.get(), &t, stepMax, input);

    return true;
  }

  static int
  ODE_system(double t, const double y[], double f[], void* params)
  {
    (void)(t); /* avoid unused parameter warning */
    double Bx  = ((double*)params)[0];
    double By  = ((double*)params)[1];
    double Bz  = ((double*)params)[2];
    double qOp = ((double*)params)[3];
    double _x  = y[0];
    double _y  = y[1];
    double _z  = y[2];
    double Tx  = y[3];
    double Ty  = y[4];
    double Tz  = y[5];
    f[0]       = Tx;
    f[1]       = Ty;
    f[2]       = Tz;
    f[3]       = qOp * (Ty * Bz - Tz * By);
    f[4]       = qOp * (Tz * Bx - Tx * Bz);
    f[5]       = qOp * (Tx * By - Ty * Bx);
    return GSL_SUCCESS;
  }
};

}  // namespace Acts

#endif  // ACTS_GSL_STEPPER_HPP
