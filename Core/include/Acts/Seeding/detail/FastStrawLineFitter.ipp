// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/detail/FastStrawLineFitter.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Enumerate.hpp"

#include <format>
namespace Acts::Experimental::detail {

template <CompositeSpacePointContainer StrawCont_t>
std::optional<FastStrawLineFitter::FitResult> FastStrawLineFitter::fit(
    const StrawCont_t& measurements, const std::vector<int>& signs) const {
  if (measurements.size() != signs.size()) {
    ACTS_WARNING(
        __func__ << "() - " << __LINE__
                 << ": Not all measurements are associated with a drift sign");
    return std::nullopt;
  }

  auto result = fit(fillAuxiliaries(measurements, signs));
  if (!result) {
    return std::nullopt;
  }
  // Calculate the chi2
  calcPostFitChi2(measurements, *result);
  return result;
}

template <CompositeSpacePointContainer StrawCont_t>
void FastStrawLineFitter::calcPostFitChi2(const StrawCont_t& measurements,
                                          FitResult& result) const {
  const TrigonomHelper angles{result.theta};
  result.chi2 = 0.;
  for (const auto& strawMeas : measurements) {
    result.chi2 += chi2Term(angles, result.y0, *strawMeas);
  }
  ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Overall chi2: "
                      << result.chi2 << ", nDoF: " << result.nDoF
                      << ", redChi2: " << (result.chi2 / result.nDoF));
}

template <CompositeSpacePoint Point_t>
double FastStrawLineFitter::chi2Term(const TrigonomHelper& angle,
                                     const double y0, const Point_t& strawMeas,
                                     std::optional<double> r) const {
  if (!strawMeas.isStraw()) {
    return 0.;
  }
  const double cov = strawMeas.covariance()[s_covIdx];
  if (cov < std::numeric_limits<double>::epsilon()) {
    return 0.;
  }
  const Vector& pos = strawMeas.localPosition();
  const double y = pos.dot(strawMeas.toNextSensor());
  const double z = pos.dot(strawMeas.planeNormal());
  const double dist = Acts::abs((y - y0) * angle.cosTheta - z * angle.sinTheta);
  ACTS_VERBOSE(__func__ << "() - " << __LINE__ << ": Distance straw (" << y
                        << ", " << z
                        << "), r: " << r.value_or(strawMeas.driftRadius())
                        << " - track: " << dist);
  return Acts::pow(dist - r.value_or(strawMeas.driftRadius()), 2) / cov;
}

template <CompositeSpacePointContainer StripCont_t>
std::optional<FastStrawLineFitter::FitResult> FastStrawLineFitter::fit(
    const StripCont_t& measurements, const ResidualIdx projection) const {
  if (projection == ResidualIdx::time) {
    ACTS_WARNING(__func__ << "() - " << __LINE__
                          << ": Only spatial projections, "
                          << "i.e. nonBending / bending are sensible");
    return std::nullopt;
  }

  FitAuxiliaries auxVars{};
  Vector centerOfGravity{Vector::Zero()};
  // @brief Selector function to check that the current strip provides a constraint in the projetor direction
  auto select = [&projection](const auto& strip) -> bool {
    // Skip straw measurements that are not twins &
    // don't measure non-bending coordinate
    if (strip->isStraw()) {
      return strip->measuresLoc0() && projection == ResidualIdx::nonBending;
    }
    // Check that the strip is actually measuring the projection
    return (strip->measuresLoc0() && projection == ResidualIdx::nonBending) ||
           (strip->measuresLoc1() && projection == ResidualIdx::bending);
  };

  auxVars.invCovs.resize(measurements.size());
  for (const auto& [sIdx, strip] : enumerate(measurements)) {
    if (!select(strip)) {
      ACTS_VERBOSE(__func__ << "() - " << __LINE__
                            << ": Skip strip measurement " << toString(*strip));
      continue;
    }
    const auto& invCov =
        (auxVars.invCovs[sIdx] =
             1. / strip->covariance()[toUnderlying(projection)]);
    auxVars.covNorm += invCov;
    centerOfGravity += invCov * strip->localPosition();
    ++auxVars.nDoF;
  }
  // To little information provided
  if (auxVars.nDoF < 3) {
    return std::nullopt;
  }
  // Reduce the number of degrees of freedom by 2 to account
  // for the two free parameters to fit
  auxVars.nDoF -= 2u;
  auxVars.covNorm = 1. / auxVars.covNorm;
  centerOfGravity *= auxVars.covNorm;

  bool centerSet{false};
  for (const auto& [sIdx, strip] : enumerate(measurements)) {
    if (!select(strip)) {
      continue;
    }
    const Vector pos = strip->localPosition() - centerOfGravity;
    const Vector& measDir{
        (projection == ResidualIdx::nonBending && strip->measuresLoc1()) ||
                strip->isStraw()
            ? strip->sensorDirection()
            : strip->toNextSensor()};

    if (!centerSet) {
      auxVars.centerY = centerOfGravity.dot(measDir);
      auxVars.centerZ = centerOfGravity.dot(strip->planeNormal());
      centerSet = true;
    }

    const double y = pos.dot(measDir);
    const double z = pos.dot(strip->planeNormal());

    const auto& invCov = auxVars.invCovs[sIdx];
    auxVars.T_zzyy += invCov * (Acts::square(z) - Acts::square(y));
    auxVars.T_yz += invCov * z * y;
  }
  return fit(auxVars);
}

template <CompositeSpacePointContainer StrawCont_t,
          CompositeSpacePointFastCalibrator<
              Acts::RemovePointer_t<typename StrawCont_t::value_type>>
              Calibrator_t>
void FastStrawLineFitter::calcPostFitChi2(const Acts::CalibrationContext& ctx,
                                          const StrawCont_t& measurements,
                                          const Calibrator_t& calibrator,
                                          FitResultT0& result) const {
  const TrigonomHelper angles{result.theta};
  result.chi2 = 0.;
  for (const auto& strawMeas : measurements) {
    result.chi2 += chi2Term(angles, result.y0, *strawMeas,
                            calibrator.driftRadius(ctx, *strawMeas, result.t0));
  }
  ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": Overall chi2: "
                      << result.chi2 << ", nDoF: " << result.nDoF
                      << ", redChi2: " << (result.chi2 / result.nDoF));
}

template <CompositeSpacePointContainer StrawCont_t>
FastStrawLineFitter::FitAuxiliaries FastStrawLineFitter::fillAuxiliaries(
    const StrawCont_t& measurements, const std::vector<int>& signs) const {
  FitAuxiliaries auxVars{};
  auxVars.invCovs.resize(signs.size(), -1.);
  Vector centerOfGravity{Vector::Zero()};

  // Calculate first the center of gravity
  for (const auto& [sIdx, strawMeas] : enumerate(measurements)) {
    if (!strawMeas->isStraw()) {
      ACTS_DEBUG(__func__ << "() - " << __LINE__ << ": The measurement "
                          << toString(*strawMeas) << " is not a straw");
      continue;
    }
    const double cov = strawMeas->covariance()[s_covIdx];
    if (cov < std::numeric_limits<double>::epsilon()) {
      ACTS_WARNING(__func__ << "() - " << __LINE__ << ": The covariance ("
                            << cov << ") of the measurement "
                            << toString(*strawMeas) << " is invalid.");
      continue;
    }
    ACTS_VERBOSE(__func__ << "() - " << __LINE__ << ": Fill "
                          << toString(*strawMeas) << ".");

    auto& invCov = (auxVars.invCovs[sIdx] = 1. / cov);
    auxVars.covNorm += invCov;
    centerOfGravity += invCov * strawMeas->localPosition();
    ++auxVars.nDoF;
  }
  if (auxVars.nDoF < 3) {
    std::stringstream sstr{};
    for (const auto& [sIdx, strawMeas] : enumerate(measurements)) {
      sstr << " --- " << (sIdx + 1) << ") " << toString(*strawMeas)
           << ", weight: " << auxVars.invCovs[sIdx] << std::endl;
    }
    ACTS_WARNING(__func__ << "() - " << __LINE__
                          << ": At least 3 measurements are required to "
                             "perform the straw line fit\n"
                          << sstr.str());
    auxVars.nDoF = 0u;
    return auxVars;
  }
  // Reduce the number of degrees of freedom by 2 to account
  // for the two free parameters to fit
  auxVars.nDoF -= 2u;
  auxVars.covNorm = 1. / auxVars.covNorm;
  centerOfGravity *= auxVars.covNorm;

  // Now calculate the fit constants
  bool centerSet{false};
  for (const auto& [sIdx, strawMeas] : enumerate(measurements)) {
    const auto& invCov = auxVars.invCovs[sIdx];
    // Invalid measurements were marked
    if (invCov < 0.) {
      continue;
    }
    if (!centerSet) {
      auxVars.centerY = centerOfGravity.dot(strawMeas->toNextSensor());
      auxVars.centerZ = centerOfGravity.dot(strawMeas->planeNormal());
      centerSet = true;
    }
    const Vector pos = strawMeas->localPosition() - centerOfGravity;
    const double y = pos.dot(strawMeas->toNextSensor());
    const double z = pos.dot(strawMeas->planeNormal());
    const double r = strawMeas->driftRadius();

    auxVars.T_zzyy += invCov * (Acts::square(z) - Acts::square(y));
    auxVars.T_yz += invCov * z * y;
    const double sInvCov = -invCov * signs[sIdx];
    auxVars.T_rz += sInvCov * z * r;
    auxVars.T_ry += sInvCov * y * r;
    auxVars.fitY0 += sInvCov * r;
  }
  auxVars.fitY0 *= auxVars.covNorm;

  return auxVars;
}

template <CompositeSpacePointContainer StrawCont_t,
          CompositeSpacePointFastCalibrator<
              Acts::RemovePointer_t<typename StrawCont_t::value_type>>
              Calibrator_t>
std::optional<FastStrawLineFitter::FitResultT0> FastStrawLineFitter::fit(
    const Acts::CalibrationContext& ctx, const Calibrator_t& calibrator,
    const StrawCont_t& measurements, const std::vector<int>& signs,
    std::optional<double> startT0) const {
  using namespace Acts::UnitLiterals;
  if (measurements.size() != signs.size()) {
    ACTS_WARNING(
        __func__ << "() - " << __LINE__
                 << ": Not all measurements are associated with a drift sign");
    return std::nullopt;
  }

  FitResultT0 result{};
  result.t0 = startT0.value_or(0.);

  FitAuxiliariesWithT0 fitPars{
      fillAuxiliaries(ctx, calibrator, measurements, signs, result.t0)};
  result.theta = startTheta(fitPars);
  result.nDoF = fitPars.nDoF;
  ACTS_DEBUG(__func__ << "() - " << __LINE__
                      << ": Initial fit parameters: " << result);
  UpdateStatus iterStatus{UpdateStatus::GoodStep};

  while ((iterStatus = updateIteration(fitPars, result)) !=
         UpdateStatus::Exceeded) {
    if (iterStatus == UpdateStatus::Converged) {
      calcPostFitChi2(ctx, measurements, calibrator, result);
      return result;
    }
    fitPars = fillAuxiliaries(ctx, calibrator, measurements, signs, result.t0);
  }
  if (logger().doPrint(Logging::VERBOSE)) {
    ACTS_VERBOSE("Fit failed, printing all measurements:");
    for (const auto& meas : measurements) {
      ACTS_VERBOSE(toString(*meas)
                   << ", t0: " << result.t0 / 1._ns
                   << ", truthR, RecoR: " << meas->driftRadius() << ", "
                   << calibrator.driftRadius(ctx, *meas, result.t0)
                   << ", velocity: "
                   << calibrator.driftVelocity(ctx, *meas, result.t0) * 1._ns
                   << ", acceleration: "
                   << calibrator.driftAcceleration(ctx, *meas, result.t0) *
                          Acts::square(1._ns));
    }
    ACTS_VERBOSE("Result: " << result);
  }
  return std::nullopt;
}

template <CompositeSpacePointContainer StrawCont_t,
          CompositeSpacePointFastCalibrator<
              Acts::RemovePointer_t<typename StrawCont_t::value_type>>
              Calibrator_t>
FastStrawLineFitter::FitAuxiliariesWithT0 FastStrawLineFitter::fillAuxiliaries(
    const CalibrationContext& ctx, const Calibrator_t& calibrator,
    const StrawCont_t& measurements, const std::vector<int>& signs,
    const double t0) const {
  using namespace Acts::UnitLiterals;
  FitAuxiliariesWithT0 auxVars{fillAuxiliaries(measurements, signs)};
  if (auxVars.nDoF <= 1) {
    auxVars.nDoF = 0;
    return auxVars;
  }
  // Account for the time offset as extra degree of freedom
  --auxVars.nDoF;
  // Fill the new extra variables
  auxVars.T_rz = 0.;
  auxVars.T_ry = 0.;
  auxVars.fitY0 = 0.;
  for (const auto& [spIdx, strawMeas] : enumerate(measurements)) {
    const double& invCov = auxVars.invCovs[spIdx];
    // Invalid (non)-straw measurements
    if (invCov < 0.) {
      continue;
    }
    const double sInvCov = -invCov * signs[spIdx];
    const double r = calibrator.driftRadius(ctx, *strawMeas, t0);
    const double v = calibrator.driftVelocity(ctx, *strawMeas, t0);
    const double a = calibrator.driftAcceleration(ctx, *strawMeas, t0);
    const double y = strawMeas->localPosition().dot(strawMeas->toNextSensor()) -
                     auxVars.centerY;
    const double z = strawMeas->localPosition().dot(strawMeas->planeNormal()) -
                     auxVars.centerZ;

    ACTS_VERBOSE(__func__ << "() - " << __LINE__ << ": # " << (spIdx + 1)
                          << ") " << toString(*strawMeas) << ", t0: "
                          << t0 / 1._ns << " r: " << r << ", v: " << v * 1._ns
                          << ", a: " << a * Acts::square(1._ns));
    auxVars.fitY0 += sInvCov * r;
    auxVars.R_v += sInvCov * v;
    auxVars.R_a += sInvCov * a;

    auxVars.T_rz += sInvCov * z * r;
    auxVars.T_ry += sInvCov * y * r;

    auxVars.T_vy += sInvCov * v * y;
    auxVars.T_vz += sInvCov * v * z;

    auxVars.R_vr += invCov * r * v;
    auxVars.R_vv += invCov * v * v;

    auxVars.T_ay += sInvCov * a * y;
    auxVars.T_az += sInvCov * a * z;

    auxVars.R_ar += invCov * a * r;
  }
  auxVars.fitY0 *= auxVars.covNorm;
  ACTS_DEBUG(__func__ << "() - " << __LINE__ << " Fit constants calculated \n"
                      << auxVars);
  return auxVars;
}

}  // namespace Acts::Experimental::detail
