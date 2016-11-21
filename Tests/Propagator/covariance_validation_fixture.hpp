// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_COVARIANCE_VALIDATION_FIXTURE_HPP
#define ACTS_COVARIANCE_VALIDATION_FIXTURE_HPP 1

#include <array>
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

      // variation in x
      std::vector<ActsVectorD<5>> x_derivatives;
      x_derivatives.reserve(h_steps.size());
      for (double h : h_steps) {
        Vector3D pos;
        Vector2D loc_pos(h, 0);
        trackPars.associatedSurface().localToGlobal(
            loc_pos, trackPars.momentum(), pos);
        CurvilinearParameters tp(
            nullptr, pos, trackPars.momentum(), trackPars.charge());
        const auto& r = m_propagator.propagate(tp, dest, var_options);
        x_derivatives.push_back((r.endParameters->parameters() - nominal) / h);
      }

      // variation in y
      std::vector<ActsVectorD<5>> y_derivatives;
      y_derivatives.reserve(h_steps.size());
      for (double h : h_steps) {
        Vector3D pos;
        Vector2D loc_pos(0, h);
        trackPars.associatedSurface().localToGlobal(
            loc_pos, trackPars.momentum(), pos);
        CurvilinearParameters tp(
            nullptr, pos, trackPars.momentum(), trackPars.charge());
        const auto& r = m_propagator.propagate(tp, dest, var_options);
        y_derivatives.push_back((r.endParameters->parameters() - nominal) / h);
      }

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
      jacobian.col(Acts::eLOC_1) = fitLinear(x_derivatives, h_steps);
      jacobian.col(Acts::eLOC_2) = fitLinear(y_derivatives, h_steps);
      jacobian.col(Acts::ePHI)   = fitLinear(phi_derivatives, h_steps);
      jacobian.col(Acts::eTHETA) = fitLinear(theta_derivatives, h_steps);
      jacobian.col(Acts::eQOP)   = fitLinear(qop_derivatives, h_steps);

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
