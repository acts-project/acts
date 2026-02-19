// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/CompositeSpacePointLineFitter.hpp"

#include "Acts/Seeding/CompositeSpacePointLineSeeder.hpp"
#include "Acts/Utilities/AlgebraHelpers.hpp"
#include "Acts/Utilities/StringHelpers.hpp"

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
                            << " -  Skip invalid measurement "
                            << toString(*sp));
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
    } else if (m_cfg.fitT0 && sp->hasTime()) {
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

template <bool fitStraws, bool fitTime, CompositeSpacePointContainer Cont_t>
CompositeSpacePointLineFitter::DelegateRet_t<fitTime>
CompositeSpacePointLineFitter::fastPrecFit(
    const Cont_t& measurements, const Line_t& initialGuess,
    const FastFitDelegate_t<Cont_t, fitTime>& delegate) const {
  std::vector<int> strawSigns{};
  if constexpr (fitStraws) {
    strawSigns = detail::CompSpacePointAuxiliaries::strawSigns(initialGuess,
                                                               measurements);
  }
  ACTS_DEBUG(__func__ << " <fit " << (fitStraws ? "straws" : "strips") << ", "
                      << (fitTime ? "with" : "no") << " time>() " << __LINE__
                      << " - Start precision fit.");
  auto result = delegate(measurements, strawSigns);

  if constexpr (!fitStraws) {
    ACTS_DEBUG(__func__ << " <fit " << (fitStraws ? "straws" : "strips") << ", "
                        << (fitTime ? "with" : "no") << " time>() " << __LINE__
                        << " - Precision fit "
                        << (result ? "succeeded" : "failed"));

    return result;
  }
  if (result) {
    ACTS_DEBUG(__func__ << " <fit " << (fitStraws ? "straws" : "strips") << ", "
                        << (fitTime ? "with" : "no") << " time>() " << __LINE__
                        << " - Precision fit converged. " << (*result));

    if (result->chi2 / static_cast<double>(result->nDoF) <
        m_cfg.badFastChi2SignSwap) {
      return result;
    }
  }

  ACTS_DEBUG(
      __func__
      << " <fit " << (fitStraws ? "straws" : "strips") << ", "
      << (fitTime ? "with" : "no") << " time>() " << __LINE__
      << " - The fit result is of poor quality attempt with L<->R swapping");
  // Swap signs

  for (int& s : strawSigns) {
    s = -s;
  }
  // Retry & check whether the chi2 is better
  auto swappedPrecResult = delegate(measurements, strawSigns);
  if (!swappedPrecResult) {
    ACTS_DEBUG(__func__ << " <fit " << (fitStraws ? "straws" : "strips") << ", "
                        << (fitTime ? "with" : "no") << " time>() " << __LINE__
                        << " - The fit did not improve");

    return result;
  }
  if (!result || swappedPrecResult->chi2 < result->chi2) {
    ACTS_DEBUG(__func__ << " <fit " << (fitStraws ? "straws" : "strips") << ", "
                        << (fitTime ? "with" : "no") << " time>() " << __LINE__
                        << " - Swapped fit is of better quality "
                        << (*swappedPrecResult));
    if (result) {
      swappedPrecResult->nIter += result->nIter;
    }
    return swappedPrecResult;
  } else {
    result->nIter += swappedPrecResult->nIter;
    ACTS_DEBUG(__func__ << " <fit " << (fitStraws ? "straws" : "strips") << ", "
                        << (fitTime ? "with" : "no") << " time>() " << __LINE__
                        << " - Fit did not improve " << (*swappedPrecResult));
  }
  return result;
}

template <bool fitStraws, bool fitTime, CompositeSpacePointContainer Cont_t>
CompositeSpacePointLineFitter::FitParameters
CompositeSpacePointLineFitter::fastFit(
    const Cont_t& measurements, const Line_t& initialGuess,
    const std::vector<FitParIndex>& parsToUse,
    const FastFitDelegate_t<Cont_t, fitTime>& precFitDelegate) const {
  FitParameters result{};

  using enum FitParIndex;

  static_assert(fitTime == false || (fitStraws == true && fitTime == true),
                "Initial t0 can only be fitted with straws");

  const bool doPrecFit = std::ranges::any_of(
      parsToUse,
      [](const FitParIndex idx) { return idx == theta || idx == y0; });
  const bool doNonPrecFit = std::ranges::any_of(
      parsToUse, [](const FitParIndex idx) { return idx == phi || idx == x0; });
  const Vector_t& preFitDir{initialGuess.direction()};
  double tanAlpha = preFitDir.x() / preFitDir.z();
  double tanBeta = preFitDir.y() / preFitDir.z();

  using precResult_t = FastFitDelegate_t<Cont_t, fitTime>::result_type;
  precResult_t precResult =
      doPrecFit ? fastPrecFit<fitStraws, fitTime>(measurements, initialGuess,
                                                  precFitDelegate)
                : std::nullopt;

  if (doPrecFit && !precResult) {
    return result;
  } else if (precResult) {
    if constexpr (fitTime) {
      if (!m_cfg.fastPreFitter) {
        result.covariance(toUnderlying(t0), toUnderlying(t0)) =
            Acts::square(precResult->dT0);
      }
      result.parameters[toUnderlying(t0)] = precResult->t0;
      if (!withinRange(precResult->t0, FitParIndex::t0)) {
        ACTS_DEBUG(__func__ << " <fit " << (fitStraws ? "straws" : "strips")
                            << ", " << (fitTime ? "with" : "no") << " time>() "
                            << __LINE__ << " - Time is out of bounds.");
        return result;
      }
    }
    if (!withinRange(precResult->y0, FitParIndex::y0) ||
        !withinRange(precResult->theta, FitParIndex::theta)) {
      ACTS_DEBUG(__func__ << " <fit " << (fitStraws ? "straws" : "strips")
                          << ", " << (fitTime ? "with" : "no") << " time>() "
                          << __LINE__
                          << " - Fit parameters are out of bounds.");
      return result;
    }
    result.parameters[toUnderlying(y0)] = precResult->y0;
    result.nDoF = precResult->nDoF;
    result.nIter = precResult->nIter;
    result.chi2 = precResult->chi2;
    result.converged = true;
    if (!m_cfg.fastPreFitter) {
      result.covariance(toUnderlying(y0), toUnderlying(y0)) =
          Acts::square(precResult->dY0);
      result.covariance(toUnderlying(theta), toUnderlying(theta)) =
          Acts::square(precResult->dTheta);
    }

    auto firstPrecMeas = std::ranges::find_if(measurements, [](const auto& m) {
      return (fitStraws && m->isStraw()) || (!fitStraws && m->measuresLoc1());
    });
    assert(firstPrecMeas != measurements.end());

    const Vector_t postFitDir = CompositeSpacePointLineSeeder::makeDirection(
        **firstPrecMeas, precResult->theta);
    tanBeta = postFitDir.y() / postFitDir.z();
  }
  // Try to perform a fast fit in non-precision direction
  const FastFitResult nonPrecResult =
      doNonPrecFit ? fastNonPrecFit(measurements) : std::nullopt;

  if (nonPrecResult) {
    tanAlpha = nonPrecResult->theta;
    result.parameters[toUnderlying(x0)] = nonPrecResult->y0;
    result.nDoF += nonPrecResult->nDoF;
    result.nIter += nonPrecResult->nIter;
    if (!m_cfg.fastPreFitter) {
      result.covariance(toUnderlying(x0), toUnderlying(x0)) =
          Acts::square(nonPrecResult->dY0);
      result.covariance(toUnderlying(phi), toUnderlying(phi)) =
          Acts::square(nonPrecResult->dTheta);
    }
  }

  const Vector_t postFitDir = makeDirectionFromAxisTangents(tanAlpha, tanBeta);
  result.parameters[toUnderlying(theta)] = VectorHelpers::theta(postFitDir);
  result.parameters[toUnderlying(phi)] = VectorHelpers::phi(postFitDir);

  return result;
}
template <CompositeSpacePointContainer Cont_t>
CompositeSpacePointLineFitter::FastFitResult
CompositeSpacePointLineFitter::fastNonPrecFit(
    const Cont_t& measurements) const {
  using enum FitParIndex;

  auto firstNonPrecMeas = std::ranges::find_if(
      measurements, [](const auto& m) { return m->measuresLoc0(); });
  if (firstNonPrecMeas == measurements.end()) {
    ACTS_DEBUG(__func__ << "() " << __LINE__
                        << " - No measurement in non-bending direction");
    return std::nullopt;
  }
  ACTS_DEBUG(__func__ << "() " << __LINE__ << " Start non precision fit.");

  FastFitResult result = m_fastFitter.fit(
      measurements, detail::CompSpacePointAuxiliaries::ResidualIdx::nonBending);
  if (!result) {
    ACTS_DEBUG(__func__ << "() " << __LINE__ << " Non precision fit failed.");
    return result;
  }
  if (!withinRange(result->y0, FitParIndex::x0)) {
    ACTS_DEBUG(__func__ << "() " << __LINE__ << " Intercept is out of bounds.");
    return std::nullopt;
  }
  const auto& refHit{**firstNonPrecMeas};
  const Vector_t& eY{!refHit.measuresLoc1() ? refHit.toNextSensor()
                                            : refHit.sensorDirection()};
  const Vector_t& eZ{refHit.planeNormal()};
  const double sinPhi = std::sin(result->theta);
  const auto dir = copySign<Vector_t, double>(
      sinPhi * eY + std::cos(result->theta) * eZ, sinPhi);
  result->theta = dir.x() / dir.z();
  return result;
}

template <CompositeSpacePointContainer Cont_t,
          CompositeSpacePointCalibrator<Cont_t, Cont_t> Calibrator_t>
CompositeSpacePointLineFitter::FitResult<Cont_t>
CompositeSpacePointLineFitter::fit(
    FitOptions<Cont_t, Calibrator_t>&& fitOpts) const {
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

    constexpr bool fastCalibrator =
        CompositeSpacePointFastCalibrator<Calibrator_t, SpacePoint_t<Cont_t>>;
    const bool fitTime{resCfg.parsToUse.back() == FitParIndex::t0 &&
                       fastCalibrator};
    if constexpr (!fastCalibrator) {
      if (resCfg.parsToUse.back() == FitParIndex::t0) {
        ACTS_WARNING(__func__ << "() " << __LINE__
                              << ": Time cannot be fast-fitted because the "
                                 "calibrator does not satisfy"
                              << " CompositeSpacePointFastCalibrator concept.");
      }
    }
    FitParameters fastResult{};
    // Straw fit
    if (nStraws >= 2) {
      if (fitTime) {
        if constexpr (fastCalibrator) {
          FastFitDelegate_t<Cont_t, true> fitDelegate{
              [this, &fitOpts, &result](const Cont_t& measurements,
                                        const std::vector<int>& strawSigns) {
                const double initialT0 =
                    result.parameters[toUnderlying(FitParIndex::t0)];
                return m_fastFitter.fit(fitOpts.calibContext,
                                        *fitOpts.calibrator, measurements,
                                        strawSigns, initialT0);
              }};
          fastResult = fastFit<true, true>(result.measurements, line,
                                           resCfg.parsToUse, fitDelegate);
        }
      } else {
        FastFitDelegate_t<Cont_t, false> fitDelegate{
            [this](const Cont_t& measurements,
                   const std::vector<int>& strawSigns) {
              return m_fastFitter.fit(measurements, strawSigns);
            }};
        fastResult = fastFit<true, false>(result.measurements, line,
                                          resCfg.parsToUse, fitDelegate);
      }
    }
    // Pure strip fit
    else {
      FastFitDelegate_t<Cont_t, false> fitDelegate{
          [this](const Cont_t& measurements,
                 const std::vector<int>& /*strawSigns*/) {
            return m_fastFitter.fit(
                measurements,
                detail::CompSpacePointAuxiliaries::ResidualIdx::bending);
          }};
      fastResult = fastFit<false, false>(result.measurements, line,
                                         resCfg.parsToUse, fitDelegate);
    }
    if (fastResult.converged) {
      if (fitTime) {
        fastResult.parameters[toUnderlying(FitParIndex::t0)] -=
            (fitOpts.localToGlobal * line.position()).norm() /
            PhysicalConstants::c;
      }
      static_cast<FitParameters&>(result) = std::move(fastResult);
      ACTS_DEBUG(__func__ << "() " << __LINE__ << " - Fast fit converged.");
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
      if (!m_cfg.ignoreFailedPreFit) {
        return result;
      }
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
        using namespace Acts::UnitLiterals;
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
        update =
            updateParameters<2>(pullCalculator.config().parsToUse.front(),
                                cache, result.parameters, result.covariance);
        break;
      }
      // 2D fit + time
      case 3: {
        update =
            updateParameters<3>(pullCalculator.config().parsToUse.front(),
                                cache, result.parameters, result.covariance);
        break;
      }
      // 3D spatial fit (x0, y0, theta, phi)
      case 4: {
        update =
            updateParameters<4>(pullCalculator.config().parsToUse.front(),
                                cache, result.parameters, result.covariance);
        break;
      }
      // full fit
      case 5: {
        update =
            updateParameters<5>(pullCalculator.config().parsToUse.front(),
                                cache, result.parameters, result.covariance);
        break;
      }
      default:
        ACTS_WARNING(__func__ << "() " << __LINE__
                              << ": Invalid parameter size "
                              << pullCalculator.config().parsToUse.size());
        return result;
    }
    /// Check whether the fit parameters are within range
    if (!withinRange(result.parameters, pullCalculator.config().parsToUse)) {
      update = UpdateStep::outOfBounds;
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
  }
  return result;
}

template <unsigned N>
CompositeSpacePointLineFitter::UpdateStep
CompositeSpacePointLineFitter::updateParameters(const FitParIndex firstPar,
                                                ChiSqCache& cache,
                                                ParamVec_t& currentPars,
                                                CovMat_t& covariance) const
  requires(N >= 2 && N <= s_nPars)
{
  auto firstIdx = toUnderlying(firstPar);
  assert(firstIdx + N <= s_nPars);
  UpdateStep retCode{UpdateStep::goodStep};
  // Current parameters mapped to an Eigen interface
  if constexpr (N == 3) {
    constexpr auto t0Idx = toUnderlying(FitParIndex::t0);
    assert(firstIdx + N - 1 == 2);
    std::swap(currentPars[2], currentPars[t0Idx]);
    std::swap(cache.gradient[2], cache.gradient[t0Idx]);
    cache.hessian.row(2).swap(cache.hessian.row(t0Idx));
    cache.hessian.col(2).swap(cache.hessian.col(t0Idx));
  }
  Eigen::Map<Vector<N>> miniPars{currentPars.data() + firstIdx};
  ACTS_VERBOSE(__func__ << "<" << N << ">() - " << __LINE__
                        << ": Current parameters " << toString(miniPars)
                        << " with chi2: " << cache.chi2 << ",  gradient: "
                        << toString(cache.gradient) << ", hessian: \n"
                        << toString(cache.hessian));

  // Take out the filled block from the gradient
  Eigen::Map<const Vector<N>> miniGradient{cache.gradient.data() + firstIdx};
  // The gradient is already small enough
  if (miniGradient.norm() < m_cfg.precCutOff) {
    ACTS_DEBUG(__func__ << "<" << N << ">() - " << __LINE__
                        << ": Gradient is small enough");
    retCode = UpdateStep::converged;
  }
  // Take out the filled block from the hessian
  Acts::SquareMatrix<N> miniHessian{
      cache.hessian.block<N, N>(firstIdx, firstIdx)};
  ACTS_VERBOSE(__func__ << "<" << N << ">() - " << __LINE__
                        << ": Projected parameters: " << toString(miniPars)
                        << " gradient: " << toString(miniGradient)
                        << ", hessian: \n"
                        << toString(miniHessian)
                        << ", determinant: " << miniHessian.determinant());

  auto inverseH = safeInverse(miniHessian);
  // The Hessian can safely be inverted
  if (inverseH) {
    const Vector<N> update{(*inverseH) * miniGradient};
    // We compute also the normalized update, defined as the parameter
    // update expressed in units of the parameter uncertainties. This quantifies
    // the significance of the update relative to the estimated errors.
    double normUpdate{0.};
    for (unsigned i = 0; i < N; ++i) {
      normUpdate += Acts::square(update[i]) / (*inverseH)(i, i);
    }

    if (std::sqrt(normUpdate) < m_cfg.normPrecCutOff ||
        update.norm() < m_cfg.precCutOff) {
      ACTS_DEBUG(__func__ << "<" << N << ">() - " << __LINE__
                          << ": Update/ Normalized update " << toString(update)
                          << "/ " << normUpdate << " are negligible small.");
      retCode = UpdateStep::converged;
    }

    // If we converged, set the covariance and return (only if Hessian is
    // invertible)
    if (retCode == UpdateStep::converged) {
      covariance.block<N, N>(firstIdx, firstIdx) = *inverseH;

      // swap back t0 component if needed
      if constexpr (N == 3) {
        constexpr auto t0Idx = toUnderlying(FitParIndex::t0);
        covariance.row(2).swap(covariance.row(t0Idx));
        covariance.col(2).swap(covariance.col(t0Idx));
        std::swap(currentPars[2], currentPars[t0Idx]);
      }
      ACTS_DEBUG(__func__ << "<" << N << ">() - " << __LINE__
                          << ": Evaluated covariance: \n"
                          << toString(covariance));
      return retCode;
    }

    ACTS_VERBOSE(__func__ << "<" << N << ">() - " << __LINE__
                          << ": Inverted Hessian \n"
                          << toString(*inverseH) << "\n-> Update parameters by "
                          << toString(update));
    miniPars -= update;

  } else if (retCode != UpdateStep::converged) {
    // Fall back to gradient decent with a fixed damping factor
    const Vector<N> update{std::min(m_cfg.gradientStep, miniGradient.norm()) *
                           miniGradient.normalized()};

    ACTS_VERBOSE(__func__ << "<" << N << ">() - " << __LINE__
                          << ": Update parameters by " << toString(update));
    miniPars -= update;
  }
  // swap back t0 component if needed
  if constexpr (N == 3) {
    std::swap(currentPars[2], currentPars[toUnderlying(FitParIndex::t0)]);
  }
  return retCode;
}
}  // namespace Acts::Experimental
