// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {

namespace IntegrationTest {

  template <typename T>
  struct covariance_validation_fixture
  {
  public:
    covariance_validation_fixture(T propagator)
      : m_propagator(std::move(propagator))
    {
    }

    /// Numerical transport of covariance using the ridder's algorithm
    /// this is for covariance propagation validation
    /// it can either be used for curvilinear transport
    template <typename StartParameters, typename EndParameters, typename U>
    ActsSymMatrixD<5>
    calculateCovariance(const StartParameters&   startPars,
                        const ActsSymMatrixD<5>& startCov,
                        const EndParameters&     endPars,
                        const U&                 options) const
    {
      // steps for estimating derivatives
      const std::array<double, 4> h_steps = {{-2e-4, -1e-4, 1e-4, 2e-4}};

      // nominal propagation
      const auto&    nominal = endPars.parameters();
      const Surface& dest    = endPars.referenceSurface();

      // - for planar surfaces the dest surface is a perfect destination
      // surface for the numerical propagation, as reference frame
      // aligns with the referenceSurface.transform().rotation() at
      // at any given time
      //
      // - for straw & cylinder, where the error is given
      // in the reference frame that re-aligns with a slightly different
      // intersection solution

      // avoid stopping before the surface because of path length reached
      U var_options = options;
      var_options.pathLimit *= 2;

      // variation in x
      std::vector<ActsVectorD<5>> x_derivatives;
      x_derivatives.reserve(h_steps.size());
      for (double h : h_steps) {
        StartParameters tp = startPars;
        tp.template set<Acts::eLOC_0>(tp.template get<Acts::eLOC_0>() + h);
        const auto& r = m_propagator.propagate(tp, dest, var_options);
        x_derivatives.push_back((r.endParameters->parameters() - nominal) / h);
      }

      // variation in y
      std::vector<ActsVectorD<5>> y_derivatives;
      y_derivatives.reserve(h_steps.size());
      for (double h : h_steps) {
        StartParameters tp = startPars;
        tp.template set<Acts::eLOC_1>(tp.template get<Acts::eLOC_1>() + h);
        const auto& r = m_propagator.propagate(tp, dest, var_options);
        y_derivatives.push_back((r.endParameters->parameters() - nominal) / h);
      }

      // variation in phi
      std::vector<ActsVectorD<5>> phi_derivatives;
      phi_derivatives.reserve(h_steps.size());
      for (double h : h_steps) {
        StartParameters tp = startPars;
        tp.template set<Acts::ePHI>(tp.template get<Acts::ePHI>() + h);
        const auto& r = m_propagator.propagate(tp, dest, var_options);
        phi_derivatives.push_back((r.endParameters->parameters() - nominal)
                                  / h);
      }

      // variation in theta
      std::vector<ActsVectorD<5>> theta_derivatives;
      theta_derivatives.reserve(h_steps.size());
      for (double h : h_steps) {
        StartParameters tp            = startPars;
        const double    current_theta = tp.template get<Acts::eTHETA>();
        if (current_theta + h > M_PI) {
          h = M_PI - current_theta;
        }
        if (current_theta + h < 0) {
          h = -current_theta;
        }
        tp.template set<Acts::eTHETA>(tp.template get<Acts::eTHETA>() + h);
        const auto& r = m_propagator.propagate(tp, dest, var_options);
        theta_derivatives.push_back((r.endParameters->parameters() - nominal)
                                    / h);
      }

      // variation in q/p
      std::vector<ActsVectorD<5>> qop_derivatives;
      qop_derivatives.reserve(h_steps.size());
      for (double h : h_steps) {
        StartParameters tp = startPars;
        tp.template set<Acts::eQOP>(tp.template get<Acts::eQOP>() + h);
        const auto& r = m_propagator.propagate(tp, dest, var_options);
        qop_derivatives.push_back((r.endParameters->parameters() - nominal)
                                  / h);
      }

      ActsSymMatrixD<5> jacobian;
      jacobian.setIdentity();
      jacobian.col(Acts::eLOC_0) = fitLinear(x_derivatives, h_steps);
      jacobian.col(Acts::eLOC_1) = fitLinear(y_derivatives, h_steps);
      jacobian.col(Acts::ePHI)   = fitLinear(phi_derivatives, h_steps);
      jacobian.col(Acts::eTHETA) = fitLinear(theta_derivatives, h_steps);
      jacobian.col(Acts::eQOP)   = fitLinear(qop_derivatives, h_steps);

      return jacobian * startCov * jacobian.transpose();
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