// This file is part of the ACTS project.
//
// Copyright (C) 2025 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/EventData/MuonSpacePointCalibrator.hpp"


namespace ActsExamples{
    using TechField = MuonSpacePoint::MuonId::TechField;


    MuonSpacePointCalibrator::MuonSpacePointCalibrator( Config&& cfg,
                                                std::unique_ptr<const Acts::Logger> logger):
        m_cfg{std::move(cfg)},
        m_logger{std::move(logger)} {}
    MuonSpacePointCalibrator::CalibSpCont_t 
        MuonSpacePointCalibrator::calibrate(const Acts::CalibrationContext& ctx,
                                            const Acts::Vector3& trackPos,
                                            const Acts::Vector3& trackDir,
                                            const double trackT0,
                                            const UnCalibSpVec_t& uncalibCont) const {
        CalibSpCont_t outContainer{};
        outContainer.reserve(uncalibCont.size());
        for (const MuonSpacePoint* calibMe : uncalibCont){
            calibrate(ctx, trackPos, trackDir, trackT0, *calibMe, outContainer);
        }
        return outContainer;
    }


    void MuonSpacePointCalibrator::calibrate(const Acts::CalibrationContext& /*ctx*/,
                                             const Acts::Vector3& trackPos,
                                             const Acts::Vector3& trackDir,
                                             const double trackT0,
                                             const MuonSpacePoint& spacePoint,
                                             CalibSpCont_t& outContainer) const {
        
        auto calibSp = std::make_unique<MuonSpacePoint>(spacePoint);
        switch (spacePoint.id().technology()) {
            case TechField::Mdt:{
                break;
            }

        }
    }

double MuonSpacePointCalibrator::driftVelocity(const Acts::CalibrationContext& /*ctx*/,
                                               const MuonSpacePoint& sp) const {
    return 0.;
}

double MuonSpacePointCalibrator::driftAcceleration(const Acts::CalibrationContext& /*ctx*/,
                                                    const MuonSpacePoint& sp) const {
    return 0.;
}
   
}