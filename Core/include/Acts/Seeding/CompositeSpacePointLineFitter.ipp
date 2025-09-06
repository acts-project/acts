// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/CompositeSpacePointLineFitter.hpp"

#include <format>
namespace Acts::Experimental {

template <CompositeSpacePointContainer Cont_t,
          CompositeSpacePointCalibrator<Cont_t, Cont_t> Calibrator_t>
CompositeSpacePointLineFitter::FitResult<Cont_t>
CompositeSpacePointLineFitter::fit(
    FitOptions<Cont_t, Calibrator_t>&& fitOpts) const {
  if (!fitOpts.calibrator) {
    throw std::invalid_argument(
        std::format("{}:{} - Please provide a valid pointer to a calibrator.",
                    __FILE__, __LINE__));
  }

  FitResult<Cont_t> result{};
  result.measurements = std::move(fitOpts.measurements);

  // Declare the auxiliaries object to calculate the residuals
  detail::CompSpacePointAuxiliaries::Config resCfg{};
  resCfg.localToGlobal = std::move(fitOpts.localToGlobal);
  resCfg.useHessian = m_cfg.useHessian;
  resCfg.calcAlongStraw = m_cfg.calcAlongStraw;
  resCfg.calcAlongStrip = m_cfg.calcAlongStrip;
  resCfg.includeToF = resCfg.includeToF;

  std::size_t nLoc0{0};
  std::size_t nLoc1{0};
  std::size_t nTime{0};
  std::size_t nStraw{0};

  /// @brief Calculate the number of degrees of freedom & deduce which
  ///        parameters are to be fitted. Returns false if the parameters
  ///        to fit change w.r.t the currently configured parameters
  auto calculateNdoF = [&]() -> bool {
    nLoc0 = nLoc1 = nTime = nStraw = 0;
    std::vector<FitParIndex> parsBefore = std::move(resCfg.parsToUse);
    for (const auto& spacePoint : result.measurements) {
      // Invalid calibration state
      if (!fitOpts.selector(*spacePoint)) {
        continue;
      }
      // Count the number of straw measurements
      if (spacePoint->isStraw()) {
        ++nStraw;
      }
      // Count the number of measurements to fit x0, phi
      if (spacePoint->measuresLoc0()) {
        ++nLoc0;
      }
      // Count the number of measurements to fit y0, theta
      if (spacePoint->measuresLoc1()) {
        ++nLoc1;
      }
      // Count the number of measurements with explicit time
      if (m_cfg.fitT0 && spacePoint->hasTime()) {
        ++nTime;
      }
    }
    // Full fit procedure -> Check which parameters need to be fitted
    if (nLoc0 > 1) {
      resCfg.parsToUse.insert(resCfg.parsToUse.end(),
                              {FitParIndex::x0, FitParIndex::phi});
    }
    // Measurements in the bending direction
    if (nLoc1 > 1) {
      resCfg.parsToUse.insert(resCfg.parsToUse.end(),
                              {FitParIndex::y0, FitParIndex::theta});
    }
    // Time measurements
    if (m_cfg.fitT0 && nTime + nStraw > 1) {
      resCfg.parsToUse.push_back(FitParIndex::t0);
    }
    std::ranges::sort(resCfg.parsToUse);

    ACTS_DEBUG(__func__ << "() " << __LINE__
                        << ": Number of measurements (loc0/loc1/time): "
                        << std::format("{:d}/{:d}/{:d}", nLoc0, nLoc1, nTime)
                        << " --> parameters to fit: ");
    return parsBefore == resCfg.parsToUse;
  };

  calculateNdoF();
  if (resCfg.parsToUse.empty()) {
    ACTS_WARNING(__func__ << "() " << __LINE__
                          << ": No valid degrees of freedom parsed. Please "
                             "check your measurements");
    return result;
  }
  // First check whether all measurements are straw and the
  // fast fitter shall be used.
  if (m_cfg.useFastFitter && nStraw == nLoc1 + nLoc0) {
  }

  /// Proceed with the usual fit
  ChiSqCache cache{};
  Line_t line{result.parameters};
  detail::CompSpacePointAuxiliaries pullCalculator{resCfg, logger().clone()};
  for (; !result.converged && result.nIter < m_cfg.nIterMax; ++result.nIter) {
    cache.reset();
  }
  return result;
}

}  // namespace Acts::Experimental
