// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/EventData/MuonSpacePointCalibrator.hpp"

#include "Acts/Surfaces/detail/LineHelper.hpp"
#include "Acts/Utilities/detail/Polynomials.hpp"

namespace ActsExamples {

using TechField = MuonSpacePoint::MuonId::TechField;

MuonSpacePointCalibrator::MuonSpacePointCalibrator(
    const Config& cfg, std::unique_ptr<const Acts::Logger> logger)
    : m_logger{std::move(logger)}, m_cfg{cfg} {}
MuonSpacePointCalibrator::CalibSpCont_t MuonSpacePointCalibrator::calibrate(
    const Acts::CalibrationContext& ctx, const Acts::Vector3& trackPos,
    const Acts::Vector3& trackDir, const double trackT0,
    const UnCalibSpVec_t& uncalibCont) const {
  CalibSpCont_t outContainer{};
  outContainer.reserve(uncalibCont.size());
  for (const MuonSpacePoint* calibMe : uncalibCont) {
    calibrate(ctx, trackPos, trackDir, trackT0, *calibMe, outContainer);
  }
  return outContainer;
}

void MuonSpacePointCalibrator::calibrate(
    const Acts::CalibrationContext& /*ctx*/, const Acts::Vector3& trackPos,
    const Acts::Vector3& trackDir, const double trackT0,
    const MuonSpacePoint& spacePoint, CalibSpCont_t& outContainer) const {
  using namespace Acts::detail::LineHelper;
  auto calibSp = std::make_unique<MuonSpacePoint>(spacePoint);
  const Acts::Intersection3D closePoint =
      lineIntersect<3>(trackPos, trackDir, spacePoint.localPosition(),
                       spacePoint.sensorDirection());

  switch (spacePoint.id().technology()) {
    using enum TechField;
    case Mdt: {
      if (closePoint.isValid()) {
        calibSp->defineCoordinates(Acts::Vector3{closePoint.position()},
                                   Acts::Vector3{spacePoint.sensorDirection()},
                                   Acts::Vector3{spacePoint.toNextSensor()});

        double driftT = driftTime(spacePoint.driftRadius());
        if (m_cfg.includeTrackT0) {
          driftT -= trackT0;
        }
        calibSp->setRadius(driftRadius(driftT));
      }
      break;
    }
    /** Micromega & sTGC technology don't re-adjust the time */
    case Tgc:
    case Mm:
    case sTgc: {
      if (closePoint.isValid() &&
          (!spacePoint.id().measuresEta() || !spacePoint.id().measuresPhi())) {
        calibSp->defineCoordinates(Acts::Vector3{closePoint.position()},
                                   Acts::Vector3{spacePoint.sensorDirection()},
                                   Acts::Vector3{spacePoint.toNextSensor()});
      }
      break;
    }
    case Rpc: {
      /** TODO: Add time adjustment taking the signal propagation time into
       * account */
      if (closePoint.isValid() &&
          (!spacePoint.id().measuresEta() || !spacePoint.id().measuresPhi())) {
        calibSp->defineCoordinates(Acts::Vector3{closePoint.position()},
                                   Acts::Vector3{spacePoint.sensorDirection()},
                                   Acts::Vector3{spacePoint.toNextSensor()});
      }
      break;
    }
    default:
      break;
  }
  outContainer.push_back(std::move(calibSp));
}

std::optional<double> MuonSpacePointCalibrator::expandPolySeries(
    const double x, unsigned derivative, const CalibPolyType polyType,
    const double upperBound, const double lowerBound,
    const std::vector<double>& coeffs) const {
  if (x < lowerBound || x > upperBound) {
    ACTS_VERBOSE("Parsed x:" << x << " is outside of the configured range ["
                             << lowerBound << ", " << upperBound << "]");
    return std::nullopt;
  }
  const double intervalNorm = 2. / (upperBound - lowerBound);
  const double xNorm = intervalNorm * (x - 0.5 * (upperBound + lowerBound));
  const double chainFactor = Acts::pow(intervalNorm, derivative);
  double result{0.};
  const auto N = static_cast<std::uint32_t>(coeffs.size());
  for (std::uint32_t k = 0; k < N; ++k) {
    switch (polyType) {
      case CalibPolyType::Chebychev:
        result +=
            Acts::detail::chebychevPolyTn(xNorm, k, derivative) * coeffs[k];
        break;
      case CalibPolyType::Legendre:
        result += Acts::detail::legendrePoly(xNorm, k, derivative) * coeffs[k];
        break;
    }
  }
  return result * chainFactor;
}

double MuonSpacePointCalibrator::driftRadius(const double driftTime) const {
  return expandPolySeries(driftTime, 0u, m_cfg.rtPolyType, m_cfg.maxDriftT,
                          m_cfg.minDriftT, m_cfg.rtCoefficients)
      .value_or(0.);
}
double MuonSpacePointCalibrator::driftRadiusUncert(
    const double driftRadius) const {
  return expandPolySeries(driftTime(driftRadius), 0u, m_cfg.rtPolyType,
                          m_cfg.maxDriftT, m_cfg.minDriftT,
                          m_cfg.rtCoefficients)
      .value_or(0.);
}
double MuonSpacePointCalibrator::driftVelocity(const double driftTime) const {
  return expandPolySeries(driftTime, 1u, m_cfg.rtPolyType, m_cfg.maxDriftT,
                          m_cfg.minDriftT, m_cfg.rtCoefficients)
      .value_or(0.);
}
double MuonSpacePointCalibrator::driftAcceleration(
    const double driftTime) const {
  return expandPolySeries(driftTime, 2u, m_cfg.rtPolyType, m_cfg.maxDriftT,
                          m_cfg.minDriftT, m_cfg.rtCoefficients)
      .value_or(m_cfg.minTubeR);
}
double MuonSpacePointCalibrator::driftTime(const double driftRadius) const {
  return expandPolySeries(driftRadius, 0, m_cfg.trPolyType, m_cfg.maxTubeR,
                          m_cfg.minTubeR, m_cfg.trCoefficients)
      .value_or(m_cfg.minDriftT);
}

double MuonSpacePointCalibrator::driftVelocity(
    const Acts::CalibrationContext& /*ctx*/, const MuonSpacePoint& sp) const {
  return sp.id().technology() == MuonSpacePoint::MuonId::TechField::Mdt
             ? driftVelocity(driftTime(sp.driftRadius()))
             : 0.;
}

double MuonSpacePointCalibrator::driftAcceleration(
    const Acts::CalibrationContext& /*ctx*/, const MuonSpacePoint& sp) const {
  return sp.id().technology() == MuonSpacePoint::MuonId::TechField::Mdt
             ? driftAcceleration(driftTime(sp.driftRadius()))
             : 0.;
}

}  // namespace ActsExamples
