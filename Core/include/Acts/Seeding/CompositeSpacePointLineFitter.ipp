// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/CompositeSpacePointLineFitter.hpp"

#include "Acts/Utilities/AlgebraHelpers.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

#include <format>
namespace Acts::Experimental {

inline std::vector<CompositeSpacePointLineFitter::FitParIndex>
CompositeSpacePointLineFitter::extractFitablePars(
    const std::array<std::size_t, 3>& hitCounts) {
  std::vector<FitParIndex> pars{};
  const auto& [nLoc0, nLoc1, nTime] = hitCounts;
  if (nLoc0 > 1) {
    pars.insert(pars.end(), {FitParIndex::x0, FitParIndex::phi});
  }
  // Measurements in the bending direction
  if (nLoc1 > 1) {
    pars.insert(pars.end(), {FitParIndex::y0, FitParIndex::theta});
  }
  // Time measurements
  if (nTime > 1) {
    pars.push_back(FitParIndex::t0);
  }
  std::ranges::sort(pars);
  return pars;
}
template <CompositeSpacePointContainer Cont_t>
std::array<std::size_t, 3> CompositeSpacePointLineFitter::countDoF(
    const Cont_t& measurements) const {
  using Sp_t = SpacePoint_t<Cont_t>;
  Selector_t<Sp_t> selector{};
  selector.template connect<&detail::passThroughSelector<Sp_t>>();
  return countDof(measurements, selector);
}

template <CompositeSpacePointContainer Cont_t>
std::array<std::size_t, 3> CompositeSpacePointLineFitter::countDoF(
    const Cont_t& measurements,
    const Selector_t<SpacePoint_t<Cont_t>>& selector) const {
  using enum detail::CompSpacePointAuxiliaries::ResidualIdx;
  auto counts = filledArray<std::size_t, 3>(0u);
  std::size_t nValid{0};
  for (const auto& sp : measurements) {
    if (!selector(*sp)) {
      continue;
    }
    ++nValid;
    if (sp->measuresLoc0()) {
      ++counts[toUnderlying(nonBending)];
    }
    if (sp->measuresLoc1()) {
      ++counts[toUnderlying(bending)];
    }
    /// Count the straws as implicit time measurements
    if (m_cfg.fitT0 && (sp->hasTime() || sp->isStraw())) {
      ++counts[toUnderlying(time)];
    }
  }
  ACTS_DEBUG(
      __func__ << "() " << __LINE__ << " - " << nValid << "/"
               << measurements.size() << " valid measurements passed. Found "
               << counts[toUnderlying(nonBending)] << "/"
               << counts[toUnderlying(bending)] << "/"
               << counts[toUnderlying(time)] << " measuring loc0/loc1/time.");
  return counts;
}

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
  using namespace Acts::UnitLiterals;

  FitResult<Cont_t> result{};
  result.measurements = std::move(fitOpts.measurements);
  result.parameters = std::move(fitOpts.startParameters);

  // Declare the auxiliaries object to calculate the residuals
  detail::CompSpacePointAuxiliaries::Config resCfg{};
  resCfg.localToGlobal = std::move(fitOpts.localToGlobal);
  resCfg.useHessian = m_cfg.useHessian;
  resCfg.calcAlongStraw = m_cfg.calcAlongStraw;
  resCfg.calcAlongStrip = m_cfg.calcAlongStrip;
  resCfg.includeToF = resCfg.includeToF;
  resCfg.parsToUse =
      extractFitablePars(countDoF(result.measurements, fitOpts.selector));
  if (resCfg.parsToUse.empty()) {
    ACTS_WARNING(__func__ << "() " << __LINE__
                          << ": No valid degrees of freedom parsed. Please "
                             "check your measurements");
    return result;
  }

  Line_t line{};

  // First check whether all measurements are straw and the
  // fast fitter shall be used.
  if (m_cfg.useFastFitter &&
      std::ranges::none_of(resCfg.parsToUse,
                           [](const FitParIndex idx) {
                             return idx == FitParIndex::x0 ||
                                    idx == FitParIndex::phi;
                           }) &&
      std::ranges::none_of(result.measurements, [&fitOpts](const auto& sp) {
        return !sp->isStraw() && fitOpts.selector(*sp);
      })) {
    detail::FastStrawLineFitter::Config fastCfg{};
    fastCfg.maxIter = m_cfg.maxIter;
    fastCfg.precCutOff = m_cfg.precCutOff;
    detail::FastStrawLineFitter fastFitter{fastCfg, logger().clone()};
    line.updateParameters(result.parameters);
    // Fit without t0
    if (resCfg.parsToUse.back() != FitParIndex::t0) {
      auto fastResult = fastFitter.fit(
          result.measurements, detail::CompSpacePointAuxiliaries::strawSigns(
                                   line, result.measurements));
      if (!fastResult) {
        ACTS_DEBUG(__func__ << "() " << __LINE__ << " - Fast fit failed.");
        return result;
      }
      // Copy the parameters & covariance
      result.parameters[toUnderlying(FitParIndex::y0)] = fastResult->y0;
      result.parameters[toUnderlying(FitParIndex::theta)] = fastResult->theta;
      result.covariance(toUnderlying(FitParIndex::y0),
                        toUnderlying(FitParIndex::y0)) =
          Acts::square(fastResult->dY0);
      result.covariance(toUnderlying(FitParIndex::theta),
                        toUnderlying(FitParIndex::theta)) =
          Acts::square(fastResult->dTheta);
      // Store information about fit quality and convergence speed
      result.nDoF = fastResult->nDoF;
      result.nIter = fastResult->nIter;
      result.chi2 = fastResult->chi2;
      result.converged = true;
      return result;
    }
  }
  // Update the drift signs if no re-calibration per iteration is scheduled
  if (!m_cfg.recalibrate) {
    line.updateParameters(result.parameters);
    fitOpts.calibrator->updateSigns(line.position(), line.direction(),
                                    result.measurements);
  }
  /// Proceed with the usual fit
  ChiSqCache cache{};
  detail::CompSpacePointAuxiliaries pullCalculator{resCfg, logger().clone()};
  for (; !result.converged && result.nIter < m_cfg.maxIter; ++result.nIter) {
    cache.reset();
    // Update the parameters from the last iteration
    line.updateParameters(result.parameters);
    const double t0 = result.parameters[toUnderlying(FitParIndex::t0)];

    ACTS_INFO(
        __func__
        << "() " << __LINE__ << ": Start iteration #" << (result.nIter + 1)
        << ", current parameters "
        << std::format(
               "theta: {:.2f}, phi: {:.2f}, y0: {:.1f}, x0: {:.1f}, t0: {:.2f}",
               result.parameters[toUnderlying(FitParIndex::theta)] / 1._degree,
               result.parameters[toUnderlying(FitParIndex::phi)] / 1._degree,
               result.parameters[toUnderlying(FitParIndex::y0)],
               result.parameters[toUnderlying(FitParIndex::x0)], t0 / 1._ns)
        << " --> " << toString(line.position()) << " + "
        << toString(line.direction()) << ", chi2: " << result.chi2 << ".");
    // Update the measurements if calibration loop is switched on
    if (m_cfg.recalibrate) {
      result.measurements = fitOpts.calibrator->calibrate(
          fitOpts.calibContext, line.position(), line.direction(), t0,
          result.measurements);
      // Check whether the measurements are still valid
      auto parsBkp = std::move(resCfg.parsToUse);
      resCfg.parsToUse =
          extractFitablePars(countDoF(result.measurements, fitOpts.selector));

      // No valid measurement is left
      if (resCfg.parsToUse.empty()) {
        ACTS_WARNING(__func__ << "() " << __LINE__ << ":  Line parameters "
                              << toString(line.position()) << " + "
                              << toString(line.direction()) << ", t0: " << t0
                              << " invalidated all measurements");
        return result;
      }

      if (parsBkp != resCfg.parsToUse) {
        // Instantiate an updated pull calculator
        pullCalculator =
            detail::CompSpacePointAuxiliaries{resCfg, logger().clone()};
      }
      // update the drift signs
      fitOpts.calibrator->updateSigns(line.position(), line.direction(),
                                      result.measurements);
    }
    // Calculate the new chi2
    for (const auto& spacePoint : result.measurements) {
      // Skip bad measurements
      if (!fitOpts.selector(*spacePoint)) {
        continue;
      }
      double driftV{0.};
      double driftA{0.};
      // Calculate the residual & derivatives
      if (resCfg.parsToUse.back() == FitParIndex::t0) {
        pullCalculator.updateFullResidual(line, t0, *spacePoint, driftV,
                                          driftA);
      } else {
        pullCalculator.updateSpatialResidual(line, *spacePoint);
      }
      // Propagte to the chi2
      pullCalculator.updateChiSq(cache, spacePoint->covariance());
    }
    pullCalculator.symmetrizeHessian(cache);
    result.chi2 = cache.chi2;
    // Now update the parameters
    UpdateStep update{UpdateStep::goodStep};
    switch (resCfg.parsToUse.size()) {
      // 2D fit (intercept + inclination angle)
      case 2: {
        update = updateParameters<2>(resCfg.parsToUse.front(), cache,
                                     result.parameters);
        break;
      }
      // 2D fit + time
      case 3: {
        update = updateParameters<3>(resCfg.parsToUse.front(), cache,
                                     result.parameters);
        break;
      }
      // 3D spatial fit (x0, y0, theta, phi)
      case 4: {
        update = updateParameters<4>(resCfg.parsToUse.front(), cache,
                                     result.parameters);
        break;
      }
      // full fit
      case 5: {
        update = updateParameters<5>(resCfg.parsToUse.front(), cache,
                                     result.parameters);
        break;
      }
      default:
        ACTS_WARNING(__func__ << "() " << __LINE__
                              << ": Invalid parameter size "
                              << resCfg.parsToUse.size());
        return result;
    }
    switch (update) {
      case UpdateStep::converged: {
        result.converged = true;
        break;
      }
      case UpdateStep::goodStep: {
        break;
      }
      case UpdateStep::outOfBounds: {
        return result;
      }
    }
  }
  // Parameters converged
  if (result.converged) {
    const auto doF = countDoF(result.measurements, fitOpts.selector);
    result.nDoF = (doF[0] + doF[1] + doF[2]) - resCfg.parsToUse.size();
    line.updateParameters(result.parameters);
    const double t0 = result.parameters[toUnderlying(FitParIndex::t0)];
    // Recalibrate the measurements before returning
    result.measurements = fitOpts.calibrator->calibrate(
        fitOpts.calibContext, line.position(), line.direction(), t0,
        result.measurements);

    // Assign the Hessian
    switch (resCfg.parsToUse.size()) {
      // 2D fit (intercept + inclination angle)
      case 2: {
        fillCovariance<2>(resCfg.parsToUse.front(), cache.hessian,
                          result.covariance);
        break;
      }
      // 2D fit + time
      case 3: {
        fillCovariance<3>(resCfg.parsToUse.front(), cache.hessian,
                          result.covariance);
        break;
      }
      // 3D spatial fit (x0, y0, theta, phi)
      case 4: {
        fillCovariance<4>(resCfg.parsToUse.front(), cache.hessian,
                          result.covariance);
        break;
      }
      // full fit
      case 5: {
        fillCovariance<5>(resCfg.parsToUse.front(), cache.hessian,
                          result.covariance);
        break;
      }
      // No need to warn here -> captured by the fit iterations
      default:
        break;
    }
  }
  return result;
}

template <unsigned N>
CompositeSpacePointLineFitter::UpdateStep
CompositeSpacePointLineFitter::updateParameters(const FitParIndex firstPar,
                                                const ChiSqCache& cache,
                                                ParamVec_t& currentPars) const
  requires(N >= 2 && N <= s_nPars)
{
  auto firstIdx = toUnderlying(firstPar);
  assert(firstIdx + N < s_nPars);
  // Current parameters mapped to an Eigen interface
  Eigen::Map<ActsVector<N>> miniPars{currentPars.data() + firstIdx};
  ACTS_INFO(__func__ << "<" << N << ">() - " << __LINE__
                     << ": Current parameters " << toString(miniPars)
                     << " with chi2: " << cache.chi2 << ",  gradient: "
                     << toString(cache.gradient) << ", hessian: \n"
                     << cache.hessian);

  // Take out the filled block from the gradient
  Eigen::Map<const ActsVector<N>> miniGradient{cache.gradient.data() +
                                               firstIdx};
  // The gradient is already small enough
  if (miniGradient.norm() < m_cfg.precCutOff) {
    ACTS_INFO(__func__ << "<" << N << ">() - " << __LINE__
                       << ": Gradient is small enough");
    return UpdateStep::converged;
  }
  // Take out the filled block from the hessian
  Acts::ActsSquareMatrix<N> miniHessian{
      cache.hessian.block<N, N>(firstIdx, firstIdx)};
  ACTS_INFO(__func__ << "<" << N << ">() - " << __LINE__
                     << ": Projected parameters: " << toString(miniPars)
                     << " gradient: " << toString(miniGradient)
                     << ", hessian: \n"
                     << miniHessian << "\n, determinant"
                     << miniHessian.determinant());

  auto inverseH = safeInverse(miniHessian);
  // The Hessian can safely be inverted
  if (inverseH) {
    const ActsVector<N> update{(*inverseH) * miniGradient};

    if (update.norm() < m_cfg.precCutOff) {
      ACTS_INFO(__func__ << "<" << N << ">() - " << __LINE__ << ": Update "
                         << toString(update) << " is negligible small.");
      return UpdateStep::converged;
    }
    ACTS_INFO(__func__ << "<" << N << ">() - " << __LINE__
                       << ": Update parameters by " << toString(update));
    miniPars -= update;

  } else {
    // Fall back to gradient decent with a fixed damping factor
    const ActsVector<N> update{
        std::min(m_cfg.gradientStep, miniGradient.norm()) *
        miniGradient.normalized()};

    ACTS_INFO(__func__ << "<" << N << ">() - " << __LINE__
                       << ": Update parameters by " << toString(update));
    miniPars -= update;
  }
  // Check parameter ranges

  return UpdateStep::goodStep;
}

template <unsigned N>
void CompositeSpacePointLineFitter::fillCovariance(const FitParIndex firstPar,
                                                   const CovMat_t& hessian,
                                                   CovMat_t& covariance) const
  requires(N >= 2 && N <= s_nPars)
{
  auto firstIdx = toUnderlying(firstPar);
  assert(firstIdx + N < s_nPars);

  Acts::ActsSquareMatrix<N> miniHessian{
      hessian.block<N, N>(firstIdx, firstIdx)};

  auto inverseH = safeInverse(miniHessian);
  // The Hessian can safely be inverted
  if (inverseH) {
    covariance.block<N, N>(firstIdx, firstIdx) = (*inverseH);
  }
}

}  // namespace Acts::Experimental
