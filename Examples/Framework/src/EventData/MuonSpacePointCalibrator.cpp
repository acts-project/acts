// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/EventData/MuonSpacePointCalibrator.hpp"

#include "Acts/Surfaces/detail/LineHelper.hpp"

namespace {

constexpr double unitIntervalPrime(const double lowerEdge,
                                   const double upperEdge) {
  return 2. / (upperEdge - lowerEdge);
}

constexpr double unitInterval(const double x, const double lowerEdge,
                              const double upperEdge) {
  return unitIntervalPrime(lowerEdge, upperEdge) *
         (x - 0.5 * (upperEdge + lowerEdge));
}

}  // namespace

namespace ActsExamples {



using TechField = MuonSpacePoint::MuonId::TechField;

MuonSpacePointCalibrator::MuonSpacePointCalibrator(
    Config&& cfg, std::unique_ptr<const Acts::Logger> logger)
    : m_cfg{std::move(cfg)}, m_logger{std::move(logger)} {}

double MuonSpacePointCalibrator::reducedTime(const double t) const {
  return unitInterval(t, m_cfg.minDriftT, m_cfg.maxDriftT);
}
double MuonSpacePointCalibrator::reducedRadius(const double r) const {
  return unitInterval(r, m_cfg.minTubeR, m_cfg.maxTubeR);
}

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
                                   Acts::Vector3{spacePoint.sensorDirection()});
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
                                   Acts::Vector3{spacePoint.sensorDirection()});
      }
      break;
    }
    case Rpc: {
      /** TODO: Add time adjustment taking the signal propagation time into
       * account */
      if (closePoint.isValid() &&
          (!spacePoint.id().measuresEta() || !spacePoint.id().measuresPhi())) {
        calibSp->defineCoordinates(Acts::Vector3{closePoint.position()},
                                   Acts::Vector3{spacePoint.sensorDirection()});
      }
      break;
    }
    default:
      break;
  }
  outContainer.push_back(std::move(calibSp));
}



double MuonSpacePointCalibrator::driftRadius(const Acts::CalibrationContext& ctx,
                                             const double driftTime) const {
  if (driftTime < m_cfg.minDriftT || driftTime > m_cfg.maxDriftT) {
     ACTS_VERBOSE("Drift time " << driftTime
                   << " is outside of the configured range ["
                   << m_cfg.minDriftT << ", " << m_cfg.maxDriftT << "]");
    return 0.;
  }
  const double reducedT = reducedTime(driftTime);
  double result{0.};
  switch (m_cfg.rtPolyType) {
    using enum CalibPolyType;
    case Chebychev: {
      break;
    } case Legendre: {
      break;
  }}
    return result;
  }
double MuonSpacePointCalibrator::driftVelocity(
    const Acts::CalibrationContext& /*ctx*/,
    const MuonSpacePoint& /*sp*/) const {
  return 0.;
}

double MuonSpacePointCalibrator::driftAcceleration(
    const Acts::CalibrationContext& /*ctx*/,
    const MuonSpacePoint& /*sp*/) const {
  return 0.;
}

//
//    double getReducedTime(const double  t) const {
//              return mapToUnitInterval(t, tLower(), tUpper());
//           }
//           /* @brief Returns the derivative term of the reduced time w.r.t. the time*/
//           double dReducedTimeDt() const {
//               return unitIntervalPrime(tLower(), tUpper());
//           }
//
//         private:
//

}  // namespace ActsExamples
