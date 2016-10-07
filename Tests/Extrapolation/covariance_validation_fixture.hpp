#ifndef ACTS_COVARIANCE_VALIDATION_FIXTURE_HPP
#define ACTS_COVARIANCE_VALIDATION_FIXTURE_HPP 1

#include <array>
#include <iostream>
#include "ACTS/EventData/TrackParameters.hpp"
#include "ACTS/Surfaces/Surface.hpp"

namespace Acts {

namespace Test {

  template <typename T>
  struct covariance_validation_fixture
  {
  public:
    covariance_validation_fixture(T propagator)
      : m_propagator(std::move(propagator))
    {
    }

    template <typename U>
    ActsSymMatrixD<5>
    calculateCovariance(const CurvilinearParameters& trackPars,
                        const U&                     options) const
    {
      // steps for estimating derivatives
      const std::array<double, 4> h_steps = {-2e-4, -1e-4, 1e-4, 2e-4};

      // nominal propagation
      const auto&    r_nominal = m_propagator.propagate(trackPars, options);
      const auto&    nominal   = r_nominal.endParameters->parameters();
      const Surface& dest      = r_nominal.endParameters->associatedSurface();

      U var_options = options;
      var_options.max_path_length *= 2;

      // variation in phi
      std::vector<ActsVectorD<5>> phi_derivatives;
      phi_derivatives.reserve(h_steps.size());
      for (double h : h_steps) {
        CurvilinearParameters tp = trackPars;
        tp.set<Acts::ePHI>(tp.get<Acts::ePHI>() + h);
        const auto& r = m_propagator.propagate(tp, dest, var_options);
        phi_derivatives.push_back((r.endParameters->parameters() - nominal)
                                  / h);
      }

      // variation in theta
      std::vector<ActsVectorD<5>> theta_derivatives;
      theta_derivatives.reserve(h_steps.size());
      for (double h : h_steps) {
        CurvilinearParameters tp = trackPars;
        tp.set<Acts::eTHETA>(tp.get<Acts::eTHETA>() + h);
        const auto& r = m_propagator.propagate(tp, dest, var_options);
        theta_derivatives.push_back((r.endParameters->parameters() - nominal)
                                    / h);
      }

      // variation in q/p
      std::vector<ActsVectorD<5>> qop_derivatives;
      qop_derivatives.reserve(h_steps.size());
      for (double h : h_steps) {
        CurvilinearParameters tp = trackPars;
        tp.set<Acts::eQOP>(tp.get<Acts::eQOP>() + h);
        const auto& r = m_propagator.propagate(tp, dest, var_options);
        qop_derivatives.push_back((r.endParameters->parameters() - nominal)
                                  / h);
      }

      ActsSymMatrixD<5> jacobian;
      jacobian.setIdentity();
      jacobian.col(Acts::ePHI)   = fitLinear(phi_derivatives, h_steps);
      jacobian.col(Acts::eTHETA) = fitLinear(theta_derivatives, h_steps);
      jacobian.col(Acts::eQOP)   = fitLinear(qop_derivatives, h_steps);

      std::cout << "J = " << jacobian << std::endl;
      return jacobian * (*trackPars.covariance()) * jacobian.transpose();
    }

  private:
    template <unsigned long int N>
    static ActsVectorD<5>
    fitLinear(const std::vector<ActsVectorD<5>>& values,
              const std::array<double, N>& h)
    {
      ActsVectorD<5> A;
      ActsVectorD<5> C;
      A.setZero();
      C.setZero();
      double B = 0;
      double D = 0;

      for (unsigned int i = 0; i < N; ++i) {
        A += h.at(i) * values.at(i);
        B += h.at(i);
        C += values.at(i);
        D += h.at(i) * h.at(i);
      }

      ActsVectorD<5> b = (N * A - B * C) / (N * D - B * B);
      ActsVectorD<5> a = (C - B * b) / N;

      return a;
    }

    T m_propagator;
  };

}  // namespace Test

}  // namespace Acts
#endif  // ACTS_COVARIANCE_VALIDATION_FIXTURE_HPP
