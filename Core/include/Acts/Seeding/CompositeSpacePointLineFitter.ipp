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

template <CompositeSpacePointContainer Cont_t>
CompositeSpacePointLineFitter::DoFcounts
CompositeSpacePointLineFitter::countDoF(const Cont_t& measurements) const {
  return countDoF(measurements, Selector_t<SpacePoint_t<Cont_t>>{});
}

template <CompositeSpacePointContainer Cont_t>
CompositeSpacePointLineFitter::DoFcounts
CompositeSpacePointLineFitter::countDoF(
    const Cont_t& measurements,
    const Selector_t<SpacePoint_t<Cont_t>>& selector) const {
  DoFcounts counts{};
  std::size_t nValid{0};
  for (const auto& sp : measurements) {
    if (selector.connected() && !selector(*sp)) {
      ACTS_VERBOSE(__func__ << "() " << __LINE__
                            << " -  Skip invalid measurement @"
                            << toString(sp->localPosition()));
      continue;
    }
    ++nValid;
    if (sp->measuresLoc0()) {
      ++counts.nonBending;
    }
    if (sp->measuresLoc1()) {
      ++counts.bending;
    }
    if (sp->isStraw()) {
      ++counts.straw;
    } else if (sp->hasTime()) {
      ++counts.time;
    }
  }
  ACTS_DEBUG(__func__ << "() " << __LINE__ << " - " << nValid << "/"
                      << measurements.size()
                      << " valid measurements passed. Found "
                      << counts.nonBending << "/" << counts.bending << "/"
                      << counts.time << " measuring loc0/loc1/time with "
                      << counts.straw << " straw measurements.");
  return counts;
}

template <CompositeSpacePointContainer Cont_t>
CompositeSpacePointLineFitter::FastFitResult
CompositeSpacePointLineFitter::fastPrecFit(
    const Cont_t& measurements, const Line_t& initialGuess,
    const std::size_t nStraws,
    const std::vector<FitParIndex>& parsToUse) const {
  if (std::ranges::none_of(parsToUse, [](const FitParIndex idx) {
        using enum FitParIndex;
        return idx == theta || idx == y0;
      })) {
    return std::nullopt;
  }
  ACTS_DEBUG(__func__ << "() " << __LINE__ << " - Fetched " << nStraws
                      << " measurements for fast fit without time.");

  if (nStraws > 2) {
    std::vector<int> signs = detail::CompSpacePointAuxiliaries::strawSigns(
        initialGuess, measurements);
    FastFitResult precResult{m_fastFitter.fit(measurements, signs)};

    // Fast fit gave a bad chi2 -> Maybe L/R solution is swapped?
    if (precResult && precResult->chi2 / static_cast<double>(precResult->nDoF) >
                          m_cfg.badFastChi2SignSwap) {
      // Swap signs
      ACTS_DEBUG(__func__ << "() " << __LINE__
                          << " - The fit result is of poor quality "
                          << (*precResult) << " attempt with L<->R swapping");
      for (auto& s : signs) {
        s = -s;
      }
      // Retry & check whether the chi2 is better
      FastFitResult swappedPrecResult = {m_fastFitter.fit(measurements, signs)};
      if (swappedPrecResult && swappedPrecResult->chi2 < precResult->chi2) {
        ACTS_DEBUG(__func__ << "() " << __LINE__
                            << " - Swapped fit is of better quality "
                            << (*swappedPrecResult));
        swappedPrecResult->nIter += precResult->nIter;
        return swappedPrecResult;
      } else {
        precResult->nIter += swappedPrecResult->nIter;
        ACTS_DEBUG(__func__ << "() " << __LINE__ << " - Fit did not improve "
                            << (*swappedPrecResult));
      }
    }
    return precResult;
  }

  // Not enough straw measurements, proceeding with strip fast fit
  using ResidualIdx = detail::CompSpacePointAuxiliaries::ResidualIdx;
  return m_fastFitter.fit(measurements, ResidualIdx::bending);
}

template <CompositeSpacePointContainer Cont_t,
          CompositeSpacePointFastCalibrator<
              CompositeSpacePointLineFitter::SpacePoint_t<Cont_t>>
              Calibrator_t>
CompositeSpacePointLineFitter::FastFitResultT0
CompositeSpacePointLineFitter::fastPrecFit(
    const Acts::CalibrationContext& ctx, const Calibrator_t& calibrator,
    const Cont_t& measurements, const Line_t& initialGuess,
    const double initialT0, const std::vector<FitParIndex>& parsToUse) const {
  if (std::ranges::none_of(parsToUse, [](const FitParIndex idx) {
        using enum FitParIndex;
        return idx == theta || idx == y0;
      })) {
    return std::nullopt;
  }
  ACTS_DEBUG(__func__ << "() " << __LINE__
                      << " Starting precision fast fit with time.");

  std::vector<int> signs =
      detail::CompSpacePointAuxiliaries::strawSigns(initialGuess, measurements);
  FastFitResultT0 precResult{
      m_fastFitter.fit(ctx, calibrator, measurements, signs, initialT0)};

  // Fast fit gave a bad chi2 -> Maybe L/R solution is swapped?
  if (precResult && precResult->chi2 / static_cast<double>(precResult->nDoF) >
                        m_cfg.badFastChi2SignSwap) {
    // Swap signs
    ACTS_DEBUG(__func__ << "() " << __LINE__
                        << " - The fit result is of poor quality "
                        << (*precResult) << " attempt with L<->R swapping");
    for (auto& s : signs) {
      s = -s;
    }
    // Retry & check whether the chi2 is better
    FastFitResultT0 swappedPrecResult{
        m_fastFitter.fit(ctx, calibrator, measurements, signs, initialT0)};
    if (swappedPrecResult && swappedPrecResult->chi2 < precResult->chi2) {
      ACTS_DEBUG(__func__ << "() " << __LINE__
                          << " - Swapped fit is of better quality "
                          << (*swappedPrecResult));
      swappedPrecResult->nIter += precResult->nIter;
      return swappedPrecResult;
    } else {
      precResult->nIter += swappedPrecResult->nIter;
      ACTS_DEBUG(__func__ << "() " << __LINE__ << " - Fit did not improve "
                          << (*swappedPrecResult));
    }
  }
  return precResult;
}

template <CompositeSpacePointContainer Cont_t>
CompositeSpacePointLineFitter::FitParameters
CompositeSpacePointLineFitter::fastFit(
    const Cont_t& measurements, const Line_t& initialGuess,
    const std::size_t nStraws,
    const std::vector<FitParIndex>& parsToUse) const {
  using enum FitParIndex;

  FitParameters result{};
  const FastFitResult precResult{
      fastPrecFit(measurements, initialGuess, nStraws, parsToUse)};
  if (!precResult) {
    ACTS_DEBUG(__func__ << "() " << __LINE__ << " - Fast fit failed.");
    return result;
  }

  // Copy spatial parameters & covariance from precision fit and
  // perform non-bending fit (if required)
  mergePrecAndNonPrec(result, precResult, measurements, parsToUse);

  return result;
}

template <CompositeSpacePointContainer Cont_t,
          CompositeSpacePointFastCalibrator<
              CompositeSpacePointLineFitter::SpacePoint_t<Cont_t>>
              Calibrator_t>
CompositeSpacePointLineFitter::FitParameters
CompositeSpacePointLineFitter::fastFit(
    const Acts::CalibrationContext& ctx, const Calibrator_t& calibrator,
    const Cont_t& measurements, const Line_t& initialGuess,
    const double startT0, const std::vector<FitParIndex>& parsToUse) const {
  using enum FitParIndex;

  FitParameters result{};
  const FastFitResultT0 precResult{fastPrecFit(
      ctx, calibrator, measurements, initialGuess, startT0, parsToUse)};
  if (!precResult) {
    ACTS_DEBUG(__func__ << "() " << __LINE__ << " - Fast fit failed.");
    return result;
  }
  // Copy t0 parameter & covariance from precision fit
  result.parameters[toUnderlying(t0)] = precResult->t0;
  result.covariance(toUnderlying(t0), toUnderlying(t0)) =
      Acts::square(precResult->dT0);

  // Copy spatial parameters & covariance from precision fit and
  // perform non-bending fit (if required)
  const bool fitT0{true};
  mergePrecAndNonPrec(result, precResult, measurements, parsToUse, fitT0);

  return result;
}

template <CompositeSpacePointContainer Cont_t>
void CompositeSpacePointLineFitter::mergePrecAndNonPrec(
    FitParameters& result, const FastFitResult& precResult,
    const Cont_t& measurements, const std::vector<FitParIndex>& parsToUse,
    const bool fitT0) const {
  using enum FitParIndex;
  using namespace Acts::UnitLiterals;

  // Copy the parameters & covariance from precision fit result
  result.parameters[toUnderlying(y0)] = precResult->y0;
  result.parameters[toUnderlying(theta)] = precResult->theta;
  result.covariance(toUnderlying(y0), toUnderlying(y0)) =
      Acts::square(precResult->dY0);
  result.covariance(toUnderlying(theta), toUnderlying(theta)) =
      Acts::square(precResult->dTheta);
  result.nDoF = static_cast<unsigned>(precResult->nDoF);
  result.nIter = static_cast<unsigned>(precResult->nIter);
  result.chi2 = precResult->chi2;
  result.converged = true;

  // Check whether a non-bending fit is required
  if (std::ranges::none_of(parsToUse, [](const FitParIndex idx) {
        return idx == phi || idx == x0;
      })) {
    result.parameters[toUnderlying(phi)] = 90._degree;
    ACTS_DEBUG(
        __func__ << "() " << __LINE__
                 << " - No measurements in non precision direction parsed.");
    return;
  }

  ACTS_DEBUG(__func__ << "() " << __LINE__
                      << " - Start fast non-precision fit.");
  using ResidualIdx = detail::CompSpacePointAuxiliaries::ResidualIdx;
  const FastFitResult nonPrecResult =
      m_fastFitter.fit(measurements, ResidualIdx::nonBending);
  if (!nonPrecResult) {
    result.parameters[toUnderlying(phi)] = 90._degree;
    ACTS_DEBUG(__func__ << "() " << __LINE__
                        << " - Fast non-precision fit failed.");
    return;
  }
  result.parameters[toUnderlying(x0)] = nonPrecResult->y0;
  result.covariance(toUnderlying(x0), toUnderlying(x0)) =
      Acts::square(nonPrecResult->dY0);
  result.nDoF += nonPrecResult->nDoF;
  result.nIter += nonPrecResult->nIter;

  // Combine the two results into a single direction
  double tanTheta = std::tan(result.parameters[toUnderlying(theta)]);
  double tanPhi = std::tan(nonPrecResult->theta);
  auto dir = makeDirectionFromAxisTangents(tanPhi, tanTheta);
  result.parameters[toUnderlying(phi)] = VectorHelpers::phi(dir);
  result.parameters[toUnderlying(theta)] = VectorHelpers::theta(dir);

  // Recompute the chi2
  Line_t recoLine{result.parameters};
  result.chi2 = 0.;
  std::ranges::for_each(
      measurements, [&recoLine, &result, &fitT0](const auto& m) {
        result.chi2 +=
            fitT0 ? detail::CompSpacePointAuxiliaries::chi2Term(
                        recoLine, result.parameters[toUnderlying(t0)], *m)
                  : detail::CompSpacePointAuxiliaries::chi2Term(recoLine, *m);
      });
  ACTS_DEBUG(__func__ << "() " << __LINE__ << ": Fast fit done. Obtained result"
                      << result);
}

template <CompositeSpacePointContainer Cont_t,
          CompositeSpacePointCalibrator<Cont_t, Cont_t> Calibrator_t>
CompositeSpacePointLineFitter::FitResult<Cont_t>
CompositeSpacePointLineFitter::fit(
    FitOptions<Cont_t, Calibrator_t>&& fitOpts) const {
  using namespace Acts::UnitLiterals;

  if (!fitOpts.calibrator) {
    throw std::invalid_argument(
        "CompositeSpacePointLineFitter::fit() - Please provide a valid pointer "
        "to a calibrator.");
  }

  FitResult<Cont_t> result{};
  result.measurements = std::move(fitOpts.measurements);
  result.parameters = std::move(fitOpts.startParameters);

  // Declare the auxiliaries object to calculate the residuals
  detail::CompSpacePointAuxiliaries::Config resCfg{};
  resCfg.localToGlobal = std::move(fitOpts.localToGlobal);
  resCfg.useHessian = m_cfg.useHessian;
  resCfg.calcAlongStraw = m_cfg.calcAlongStraw;
  resCfg.calcAlongStrip = m_cfg.calcAlongStrip;
  resCfg.includeToF = m_cfg.includeToF;

  DoFcounts hitCounts{countDoF(result.measurements, fitOpts.selector)};
  const std::size_t nStraws{hitCounts.straw};
  resCfg.parsToUse = extractFitablePars(hitCounts);
  if (resCfg.parsToUse.empty()) {
    ACTS_WARNING(__func__ << "() " << __LINE__
                          << ": No valid degrees of freedom parsed. Please "
                          << "check your measurements");
    return result;
  }

  Line_t line{};

  // Fast fitter shall be used
  if (m_cfg.useFastFitter) {
    line.updateParameters(result.parameters);
    ACTS_DEBUG(__func__ << "() " << __LINE__
                        << " - Attempt a fast fit, first.");
    const bool fitTime{resCfg.parsToUse.back() == FitParIndex::t0};
    FitParameters fastResult{};
    if (fitTime && nStraws > 2) {
      if constexpr (CompositeSpacePointFastCalibrator<Calibrator_t,
                                                      SpacePoint_t<Cont_t>>) {
        fastResult = fastFit(fitOpts.calibContext, *fitOpts.calibrator,
                             result.measurements, line,
                             result.parameters[toUnderlying(FitParIndex::t0)],
                             resCfg.parsToUse);
        if (fastResult.converged) {
          fastResult.parameters[toUnderlying(FitParIndex::t0)] -=
              (fitOpts.localToGlobal * line.position()).norm() /
              PhysicalConstants::c;
        }
      } else {
        ACTS_WARNING(__func__
                     << "() " << __LINE__
                     << ": Skipping fast fit - calibrator does not satisfy "
                     << "CompositeSpacePointFastCalibrator concept");
      }
    } else {
      fastResult =
          fastFit(result.measurements, line, nStraws, resCfg.parsToUse);
    }
    if (fastResult.converged) {
      static_cast<FitParameters&>(result) = std::move(fastResult);
      ACTS_DEBUG(__func__ << "() " << __LINE__ << " - Fit converged.");
      // Use the result from the fast fitter as final answer
      if (!m_cfg.fastPreFitter) {
        line.updateParameters(result.parameters);
        // Recalibrate the measurements before returning
        result.measurements = fitOpts.calibrator->calibrate(
            fitOpts.calibContext, line.position(), line.direction(),
            fitTime ? result.parameters[toUnderlying(FitParIndex::t0)] : 0.,
            result.measurements);
        return result;
      }
      // Set convergence flag to false because the full fit comes later.
      result.converged = false;
    } else {
      ACTS_DEBUG(__func__ << "() " << __LINE__ << " - Fit failed.");
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

    ACTS_VERBOSE(__func__ << "() " << __LINE__ << ": Start next iteration "
                          << result << " --> " << toString(line.position())
                          << " + " << toString(line.direction()) << ".");
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
                              << toString(line.direction()) << ", t0: "
                              << t0 / 1._ns << " invalidated all measurements");
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
      if (fitOpts.selector.connected() && !fitOpts.selector(*spacePoint)) {
        continue;
      }
      // Calculate the residual & derivatives
      if (pullCalculator.config().parsToUse.back() == FitParIndex::t0) {
        double driftV{fitOpts.calibrator->driftVelocity(fitOpts.calibContext,
                                                        *spacePoint)};
        double driftA{fitOpts.calibrator->driftAcceleration(
            fitOpts.calibContext, *spacePoint)};
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
    switch (pullCalculator.config().parsToUse.size()) {
      // 2D fit (intercept + inclination angle)
      case 2: {
        update = updateParameters<2>(pullCalculator.config().parsToUse.front(),
                                     cache, result.parameters);
        break;
      }
      // 2D fit + time
      case 3: {
        update = updateParameters<3>(pullCalculator.config().parsToUse.front(),
                                     cache, result.parameters);
        break;
      }
      // 3D spatial fit (x0, y0, theta, phi)
      case 4: {
        update = updateParameters<4>(pullCalculator.config().parsToUse.front(),
                                     cache, result.parameters);
        break;
      }
      // full fit
      case 5: {
        update = updateParameters<5>(pullCalculator.config().parsToUse.front(),
                                     cache, result.parameters);
        break;
      }
      default:
        ACTS_WARNING(__func__ << "() " << __LINE__
                              << ": Invalid parameter size "
                              << pullCalculator.config().parsToUse.size());
        return result;
    }
    /// Check whether the fit parameters are within range
    for (const FitParIndex par : pullCalculator.config().parsToUse) {
      const auto p = toUnderlying(par);
      if (m_cfg.ranges[p][0] < m_cfg.ranges[p][1] &&
          (result.parameters[p] < m_cfg.ranges[p][0] ||
           result.parameters[p] > m_cfg.ranges[p][1])) {
        ACTS_VERBOSE(__func__ << "() " << __LINE__ << ": The parameter "
                              << pullCalculator.parName(par) << " "
                              << result.parameters[p] << " is out range ["
                              << m_cfg.ranges[p][0] << ";" << m_cfg.ranges[p][1]
                              << "]");
        update = UpdateStep::outOfBounds;
      }
    }
    switch (update) {
      using enum UpdateStep;
      case converged: {
        result.converged = true;
        break;
      }
      case goodStep: {
        break;
      }
      case outOfBounds: {
        return result;
      }
    }
  }
  // Parameters converged
  if (result.converged) {
    const DoFcounts doF = countDoF(result.measurements, fitOpts.selector);
    result.nDoF = doF.nonBending + doF.bending + doF.time -
                  pullCalculator.config().parsToUse.size();
    line.updateParameters(result.parameters);
    const double t0 = result.parameters[toUnderlying(FitParIndex::t0)];
    // Recalibrate the measurements before returning
    result.measurements = fitOpts.calibrator->calibrate(
        fitOpts.calibContext, line.position(), line.direction(), t0,
        result.measurements);

    // Assign the Hessian
    switch (pullCalculator.config().parsToUse.size()) {
      // 2D fit (intercept + inclination angle)
      case 2: {
        fillCovariance<2>(pullCalculator.config().parsToUse.front(),
                          cache.hessian, result.covariance);
        break;
      }
      // 2D fit + time
      case 3: {
        fillCovariance<3>(pullCalculator.config().parsToUse.front(),
                          cache.hessian, result.covariance);
        break;
      }
      // 3D spatial fit (x0, y0, theta, phi)
      case 4: {
        fillCovariance<4>(pullCalculator.config().parsToUse.front(),
                          cache.hessian, result.covariance);
        break;
      }
      // full fit
      case 5: {
        fillCovariance<5>(pullCalculator.config().parsToUse.front(),
                          cache.hessian, result.covariance);
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
  assert(firstIdx + N <= s_nPars);

  // Current parameters mapped to an Eigen interface
  if constexpr (N == 3) {
    std::swap(currentPars[firstIdx + N - 1],
              currentPars[toUnderlying(FitParIndex::t0)]);
  }
  Eigen::Map<ActsVector<N>> miniPars{currentPars.data() + firstIdx};
  ACTS_VERBOSE(__func__ << "<" << N << ">() - " << __LINE__
                        << ": Current parameters " << toString(miniPars)
                        << " with chi2: " << cache.chi2 << ",  gradient: "
                        << toString(cache.gradient) << ", hessian: \n"
                        << cache.hessian);

  // Define a lambda for the common update logic
  auto doUpdate = [this, &miniPars](
                      const Eigen::Map<const ActsVector<N>>& miniGradient,
                      const ActsSquareMatrix<N>& miniHessian) {
    ACTS_VERBOSE(__func__ << "<" << N << ">() - " << __LINE__
                          << ": Projected parameters: " << toString(miniPars)
                          << " gradient: " << toString(miniGradient)
                          << ", hessian: \n"
                          << miniHessian
                          << ", determinant: " << miniHessian.determinant());
    // The Hessian can safely be inverted
    if (auto inverseH = safeInverse(miniHessian); inverseH != std::nullopt) {
      const ActsVector<N> update{(*inverseH) * miniGradient};
      if (update.norm() < m_cfg.precCutOff) {
        ACTS_DEBUG(__func__ << "<" << N << ">() - " << __LINE__ << ": Update "
                            << toString(update) << " is negligible small.");
        return UpdateStep::converged;
      }
      ACTS_VERBOSE(__func__ << "<" << N << ">() - " << __LINE__ << ": Update: "
                            << toString(update) << " with inverted Hessian \n"
                            << (*inverseH));
      miniPars -= update;
    } else {
      // Fall back to gradient decent with a fixed damping factor
      const ActsVector<N> update{
          std::min(m_cfg.gradientStep, miniGradient.norm()) *
          miniGradient.normalized()};
      ACTS_VERBOSE(__func__ << "<" << N << ">() - " << __LINE__
                            << ": Update parameters by " << toString(update));
      miniPars -= update;
    }
    return UpdateStep::goodStep;
  };

  // Handle the case N=3 separately
  if constexpr (N != 3) {
    // Take out the filled block from the gradient
    Eigen::Map<const ActsVector<N>> miniGradient{cache.gradient.data() +
                                                 firstIdx};
    // The gradient is already small enough
    if (miniGradient.norm() < m_cfg.precCutOff) {
      ACTS_DEBUG(__func__ << "<" << N << ">() - " << __LINE__
                          << ": Gradient is small enough");
      return UpdateStep::converged;
    }
    // Take out the filled block from the hessian
    const ActsSquareMatrix<N> miniHessian{
        cache.hessian.block<N, N>(firstIdx, firstIdx)};

    return doUpdate(miniGradient, miniHessian);

  } else {
    // Take out the filled block from the gradient
    Eigen::Array<unsigned, N, 1> idx{firstIdx, firstIdx + 1u,
                                     toUnderlying(FitParIndex::t0)};
    const ActsVector<N> miniGrad{cache.gradient(idx)};
    Eigen::Map<const ActsVector<N>> miniGradient{miniGrad.data()};

    // The gradient is already small enough
    if (miniGradient.norm() < m_cfg.precCutOff) {
      ACTS_DEBUG(__func__ << "<" << N << ">() - " << __LINE__
                          << ": Gradient is small enough");
      return UpdateStep::converged;
    }

    // Take out the filled block from the hessian
    const Acts::ActsSquareMatrix<N> miniHessian{cache.hessian(idx, idx)};

    auto result{doUpdate(miniGradient, miniHessian)};

    // Swap back the updated parameters
    std::swap(currentPars[firstIdx + N - 1],
              currentPars[toUnderlying(FitParIndex::t0)]);

    return result;
  }
}

template <unsigned N>
void CompositeSpacePointLineFitter::fillCovariance(const FitParIndex firstPar,
                                                   const CovMat_t& hessian,
                                                   CovMat_t& covariance) const
  requires(N >= 2 && N <= s_nPars)
{
  auto firstIdx = toUnderlying(firstPar);
  assert(firstIdx + N <= s_nPars);

  if constexpr (N == 3) {
    Eigen::Array<unsigned, N, 1> idx{firstIdx, firstIdx + 1u,
                                     toUnderlying(FitParIndex::t0)};
    Acts::ActsSquareMatrix<N> miniHessian = hessian(idx, idx);

    auto inverseH = safeInverse(miniHessian);
    // The Hessian can safely be inverted
    if (inverseH) {
      covariance(idx, idx) = *inverseH;
    }

  } else {
    Acts::ActsSquareMatrix<N> miniHessian =
        hessian.block<N, N>(firstIdx, firstIdx).eval();

    auto inverseH = safeInverse(miniHessian);
    // The Hessian can safely be inverted
    if (inverseH) {
      auto covBlock = covariance.block<N, N>(firstIdx, firstIdx);
      covBlock = *inverseH;
    }
  }
  ACTS_DEBUG(__func__ << "<" << N << ">() - " << __LINE__
                      << ": Evaluated covariance: \n"
                      << covariance);
}

}  // namespace Acts::Experimental
