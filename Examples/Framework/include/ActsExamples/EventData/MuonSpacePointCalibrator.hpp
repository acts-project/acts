// This file is part of the ACTS project.
//
// Copyright (C) 2025 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#pragma once

#include "ActsExamples/EventData/MuonSpacePoint.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <memory>
namespace ActsExamples{

    class MuonSpacePointCalibrator {
        public:
            /** @brief Calibration configuration struct */
            struct Config{
                /** @brief Coefficients of the r-t chebychev polynomial */
                std::vector<double> rtChebyCoeff{};
                /** @brief Minimum drift time covered by the straw calibration */
                double minDriftT{0.};
                /** @brief Maximum drift time covered by the straw calibration */
                double maxDriftT{0.};
            };

            /** @brief default constructor */
            MuonSpacePointCalibrator( Config&& cfg,
                                      std::unique_ptr<const Acts::Logger> logger);
            /** @brief Abbrivation of the uncalibrated space point container type */
            using UnCalibSpVec_t = MuonSpacePointSorter::SpVec_t;
            /** @brief Abbrivation to access a single element in the uncalibrated spacepoint container */
            using UnCalibSp_t = UnCalibSpVec_t::value_type;
            using CalibSp_t = std::unique_ptr<MuonSpacePoint>;
            using CalibSpCont_t = std::vector<CalibSp_t>;
            
            
            CalibSpCont_t calibrate(const Acts::CalibrationContext& ctx,
                                    const Acts::Vector3& trackPos,
                                    const Acts::Vector3& trackDir,
                                    const double trackT0,
                                    const UnCalibSpVec_t& uncalibCont) const;
        

            void calibrate(const Acts::CalibrationContext& ctx,
                                    const Acts::Vector3& trackPos,
                                    const Acts::Vector3& trackDir,
                                    const double trackT0,
                                    const MuonSpacePoint& spacePoint,
                                    CalibSpCont_t& outContainer) const;
            
            double driftVelocity(const Acts::CalibrationContext& ctx,
                                 const MuonSpacePoint& sp) const;

            double driftAcceleration(const Acts::CalibrationContext& ctx,
                                     const MuonSpacePoint& sp) const;
        private:
            const Config m_cfg{};
            std::unique_ptr<const Acts::Logger> m_logger{};


    };
}