// This file is part of the Acts project.
//
// Copyright (C) 2017-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {

namespace IntegrationTest {

using Jacobian = BoundMatrix;
using Covariance = BoundSymMatrix;

template <typename T>
struct covariance_validation_fixture {
 public:
  covariance_validation_fixture(T propagator)
      : m_propagator(std::move(propagator)) {}

	bool
	inconsistentDerivativesOnDisc(std::vector<BoundVector>& derivatives) const
	{
		for(unsigned int i = 0; i < derivatives.size(); i++)
		{
			bool jumpedAngle = true;
			for(unsigned int j = 0; j < derivatives.size(); j++)
			{
				if(i != j && std::abs(derivatives[i](1) - derivatives[j](1)) < 0.5 * M_PI)
				{
					jumpedAngle = false;
					break;
				}
			}
			if(jumpedAngle)
			{
				return true;
			}
		}
		return false;
	}

  /// Numerical transport of covariance using the ridder's algorithm
  /// this is for covariance propagation validation
  /// it can either be used for curvilinear transport
  template <typename StartParameters, typename EndParameters, typename U>
  Covariance calculateCovariance(const StartParameters& startPars,
                                 const Covariance& startCov,
                                 const EndParameters& endPars,
                                 const U& options) const {
    // steps for estimating derivatives
    const std::array<double, 4> h_steps = {{-4e-4, -2e-4, 2e-4, 4e-4}};

    // nominal propagation
    const auto& nominal = endPars.parameters();
    const Surface& dest = endPars.referenceSurface();
std::cout << "nominal: " << nominal.transpose() << " | " << endPars.momentum().norm() << " | " << dest.center(options.geoContext).transpose() << std::endl;
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
    std::vector<BoundVector> x_derivatives;
    x_derivatives.reserve(h_steps.size());
    for (double h : h_steps) {
      StartParameters tp = startPars;
      tp.template set<Acts::eLOC_0>(options.geoContext,
                                    tp.template get<Acts::eLOC_0>() + h);
      const auto& r = m_propagator.propagate(tp, dest, var_options).value();
      x_derivatives.push_back((r.endParameters->parameters() - nominal) / h);
    }

    // variation in y
    std::vector<BoundVector> y_derivatives;
    y_derivatives.reserve(h_steps.size());
    for (double h : h_steps) {
      StartParameters tp = startPars;
      tp.template set<Acts::eLOC_1>(options.geoContext,
                                    tp.template get<Acts::eLOC_1>() + h);
      const auto& r = m_propagator.propagate(tp, dest, var_options).value();
      y_derivatives.push_back((r.endParameters->parameters() - nominal) / h);
    }

    // variation in phi
    std::vector<BoundVector> phi_derivatives;
    phi_derivatives.reserve(h_steps.size());
    for (double h : h_steps) {
      StartParameters tp = startPars;
      tp.template set<Acts::ePHI>(options.geoContext,
                                  tp.template get<Acts::ePHI>() + h);
      const auto& r = m_propagator.propagate(tp, dest, var_options).value();
      phi_derivatives.push_back((r.endParameters->parameters() - nominal) / h);
      
      double phi0 = nominal(Acts::ePHI);
      double phi1 = r.endParameters->parameters()(Acts::ePHI);
      if(std::abs(phi1 + 2. * M_PI - phi0) < std::abs(phi1 - phi0))
		phi_derivatives.back()[Acts::ePHI] = (phi1 + 2. * M_PI - phi0) / h;
	  else
		if(std::abs(phi1 - 2. * M_PI - phi0) < std::abs(phi1 - phi0))
			phi_derivatives.back()[Acts::ePHI] = (phi1 - 2. * M_PI - phi0) / h;
std::cout << "phi: " << tp.template get<Acts::ePHI>() << " | " << r.endParameters->parameters().transpose() << " | " <<  phi_derivatives.back().transpose() << " | " << std::endl;
    }

    // variation in theta
    std::vector<BoundVector> theta_derivatives;
    theta_derivatives.reserve(h_steps.size());
    for (double h : h_steps) {
      StartParameters tp = startPars;
      const double current_theta = tp.template get<Acts::eTHETA>();
      if (current_theta + h > M_PI) {
        h = M_PI - current_theta;
      }
      if (current_theta + h < 0) {
        h = -current_theta;
      }
      tp.template set<Acts::eTHETA>(options.geoContext,
                                    tp.template get<Acts::eTHETA>() + h);
      const auto& r = m_propagator.propagate(tp, dest, var_options).value();
      theta_derivatives.push_back((r.endParameters->parameters() - nominal) /
                                  h);
    }

    // variation in q/p
    std::vector<BoundVector> qop_derivatives;
    qop_derivatives.reserve(h_steps.size());
    for (double h : h_steps) {
      StartParameters tp = startPars;
      tp.template set<Acts::eQOP>(options.geoContext,
                                  tp.template get<Acts::eQOP>() + h);
      const auto& r = m_propagator.propagate(tp, dest, var_options).value();
      qop_derivatives.push_back((r.endParameters->parameters() - nominal) / h);
    }

    // variation in t
    std::vector<BoundVector> t_derivatives;
    t_derivatives.reserve(h_steps.size());
    for (double h : h_steps) {
      StartParameters tp = startPars;
      tp.template set<Acts::eT>(options.geoContext,
                                tp.template get<Acts::eT>() + h);
      const auto& r = m_propagator.propagate(tp, dest, var_options).value();
      t_derivatives.push_back((r.endParameters->parameters() - nominal) / h);
    }

	if(dest.type() == Surface::Disc &&
		(inconsistentDerivativesOnDisc(x_derivatives) ||
		inconsistentDerivativesOnDisc(y_derivatives) ||
		inconsistentDerivativesOnDisc(phi_derivatives) ||
		inconsistentDerivativesOnDisc(theta_derivatives) ||
		inconsistentDerivativesOnDisc(qop_derivatives) ||
		inconsistentDerivativesOnDisc(t_derivatives)))
		{
		return startCov;
	}
	
    Jacobian jacobian;
    jacobian.setIdentity();
    jacobian.col(Acts::eLOC_0) = fitLinear(x_derivatives, h_steps);
    jacobian.col(Acts::eLOC_1) = fitLinear(y_derivatives, h_steps);
    jacobian.col(Acts::ePHI) = fitLinear(phi_derivatives, h_steps);
    jacobian.col(Acts::eTHETA) = fitLinear(theta_derivatives, h_steps);
    jacobian.col(Acts::eQOP) = fitLinear(qop_derivatives, h_steps);
    jacobian.col(Acts::eT) = fitLinear(t_derivatives, h_steps);
std::cout << "jac:\n" << jacobian << std::endl;
//~ std::cout << "initial cov:\n" << startCov << std::endl;
//~ std::cout << "cov:\n" << jacobian * startCov * jacobian.transpose() << std::endl;
    return jacobian * startCov * jacobian.transpose();
  }

 private:
  template <unsigned long int N>
  static BoundVector fitLinear(const std::vector<BoundVector>& values,
                               const std::array<double, N>& h) {
    BoundVector A;
    BoundVector C;
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

    BoundVector b = (N * A - B * C) / (N * D - B * B);
    BoundVector a = (C - B * b) / N;

    return a;
  }

  T m_propagator;
};

}  // namespace IntegrationTest

}  // namespace Acts