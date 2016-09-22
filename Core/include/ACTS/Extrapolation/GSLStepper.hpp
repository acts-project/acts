#ifndef ACTS_GSL_STEPPER_HPP
#define ACTS_GSL_STEPPER_HPP 1

#include <array>
#include <gsl/gsl_errno.h>
#include <gsl/gsl_odeiv2.h>
#include <memory>
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Utilities/Units.hpp"

namespace Acts {

struct DefaultParConverter
{
  template <typename T>
  using return_parameter_type = T;

  template <typename T>
  struct cache_type
  {
    cache_type(const T& trackPar)
      : params{0, 0, 0, 0, 0, 0}, charge(trackPar.charge()), p(0)
    {
      const auto& pos = trackPar.position();
      const auto& mom = trackPar.momentum();
      p               = mom.norm();
      params[0]       = pos(0);
      params[1]       = pos(1);
      params[2]       = pos(2);
      params[3]       = mom(0) / p;
      params[4]       = mom(1) / p;
      params[5]       = mom(2) / p;
    }

    double params[6];
    double charge;
    double p;
  };

  template <typename T>
  static return_parameter_type<T>
  convert(const cache_type<T>& cache)
  {
    return CurvilinearParameters(
        nullptr,
        Vector3D(cache.params[0], cache.params[1], cache.params[2]),
        cache.p * Vector3D(cache.params[3], cache.params[4], cache.params[5]),
        cache.charge);
  }
};

template <typename BField, typename ParConverter = DefaultParConverter>
class GSLStepper
{
public:
  template <typename T>
  using cache_type = typename ParConverter::template cache_type<T>;

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
  double
  doStep(cache_type<TrackParameters>& cache, double stepMax = 1 * units::_cm)
  {
    // set B-Field and q/p
    m_bField(cache.params, m_params.data());
    // conversion to SI units: p --> (p/c)/J / (kg * m^2 / s^2)
    m_params[3] = cache.charge / units::Nat2SI<units::MOMENTUM>(cache.p);

    return performStep(cache.params, stepMax);
  }

private:
  BField m_bField;
  std::array<double, 4> m_params;
  gsl_odeiv2_system                  m_ODE_system;
  std::unique_ptr<gsl_odeiv2_driver> m_GSL_driver;

  /// input = [x,y,z,Tx,Ty,Tz,q/p]
  double
  performStep(double input[6], double stepMax) const
  {
    double t = 0;
    gsl_odeiv2_driver_apply(m_GSL_driver.get(), &t, stepMax, input);

    return t;
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
