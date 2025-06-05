// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/EventData/MuonSpacePointCalibrator.hpp"

#include "Acts/Surfaces/detail/LineHelper.hpp"

namespace ActsExamples {
using TechField = MuonSpacePoint::MuonId::TechField;

MuonSpacePointCalibrator::MuonSpacePointCalibrator(
    Config&& cfg, std::unique_ptr<const Acts::Logger> logger)
    : m_cfg{std::move(cfg)}, m_logger{std::move(logger)} {}
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
  auto calibSp = std::make_unique<MuonSpacePoint>(spacePoint);
  const Acts::Intersection3D closePoint = Acts::LineHelper::lineIntersect<3>(
      trackPos, trackDir, spacePoint.localPosition(),
      spacePoint.sensorDirection());

  switch (spacePoint.id().technology()) {
    case TechField::Mdt: {
      if (closePoint.isValid()) {
        calibSp->defineCoordinates(Acts::Vector3{closePoint.position()},
                                   Acts::Vector3{spacePoint.sensorDirection()});
      }
      break;
    }
    /** Micromega & sTGC technology don't readjust the time */
    case TechField::Mm:
    case TechField::Tgc:
    case TechField::sTgc:{
        if (closePoint.isValid() && (!spacePoint.id().measuresEta() || !spacePoint.id().measuresPhi())) {
            calibSp->defineCoordinates(Acts::Vector3{closePoint.position()},
                                   Acts::Vector3{spacePoint.sensorDirection()});
        }
        break;
    }
    case TechField::Rpc:{
        /** TODO: Add time adjustment taking the signal propagation time into account */
        if (closePoint.isValid() && (!spacePoint.id().measuresEta() || !spacePoint.id().measuresPhi())) {
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

}  // namespace ActsExamples