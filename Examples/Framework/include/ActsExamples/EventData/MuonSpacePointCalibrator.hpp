// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/MuonSPLayerSplitter.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"

#include <memory>
namespace ActsExamples {

class MuonSpacePointCalibrator {
 public:
  enum class CalibPolyType { Chebychev, Legendre };
  /// @brief Calibration configuration struct */
  struct Config {
    /// @brief Coefficients of the r-t relation parametrised as othogonal polynomial
    std::vector<double> rtCoefficients{};
    ///  @brief Coefficients of the t-r relation as othogonal polynomial
    std::vector<double> trCoefficients{};
    /// @brief Coefficients of the drift radius uncertainty parameterized as an othogonal polynomial
    ///        as a function of the drift time
    std::vector<double> rtUncertCoefficients{};
    
    /// @brief Type of the calibration polynomial
    CalibPolyType rtPolyType{CalibPolyType::Chebychev};
    /// @brief Type of the calibration polynomial
    CalibPolyType trPolyType{CalibPolyType::Chebychev};
    /// @brief Type of the calibration polynomial
    CalibPolyType rtUncertPolyType{CalibPolyType::Chebychev};
    


    /// @brief Minimum drift time covered by the r-t relation
    double minDriftT{0.};
    /// @brief Maximum drift time covered by the r-t relation
    double maxDriftT{0.};
    /// @brief Minimum radius mapped by the t-r relation
    double minTubeR{0.};
    /// @brief Maximum radius mapped by the t-r relation
    double maxTubeR{0.};
    /// @brief Propagation velocity of the electronic signal along the wire
    double propagationV{0.};
  };

  /// @brief default constructor */
  MuonSpacePointCalibrator(Config&& cfg,
                           std::unique_ptr<const Acts::Logger> logger);
  /// @brief Abbrivation of the uncalibrated space point container type
  using UnCalibSpVec_t = MuonSPLayerSplitter::SpVec_t;
  /// @brief Abbrivation to access a single element in the uncalibrated spacepoint container
  using UnCalibSp_t = UnCalibSpVec_t::value_type;
  /// @brief Use unique pointer for the calibrated measurements
  using CalibSp_t = std::unique_ptr<MuonSpacePoint>;
  /// @brief Use a simple vector of unique_ptrs as container type
  using CalibSpCont_t = std::vector<CalibSp_t>;

  CalibSpCont_t calibrate(const Acts::CalibrationContext& ctx,
                          const Acts::Vector3& trackPos,
                          const Acts::Vector3& trackDir, const double trackT0,
                          const UnCalibSpVec_t& uncalibCont) const;

  CalibSpCont_t calibrate(const Acts::CalibrationContext& ctx,
                          const Acts::Vector3& trackPos,
                          const Acts::Vector3& trackDir, const double trackT0,
                          const CalibSpCont_t& uncalibCont) const;

  void calibrate(const Acts::CalibrationContext& ctx,
                 const Acts::Vector3& trackPos, const Acts::Vector3& trackDir,
                 const double trackT0, const MuonSpacePoint& spacePoint,
                 CalibSpCont_t& outContainer) const;

  double driftVelocity(const Acts::CalibrationContext& ctx,
                       const MuonSpacePoint& sp) const;

  double driftAcceleration(const Acts::CalibrationContext& ctx,
                           const MuonSpacePoint& sp) const;

  /// @brief Calculates the drift radius from the drift time
  /// @param ctx: Calibration context
  /// @param driftTime: Drift time in nanoseconds
  double driftRadius(const Acts::CalibrationContext& ctx,
                     const double driftTime) const;
 private:
  /// @brief Maps the drift time to the interval [-1; 1.]
  ///        w.r.t configured min & max drift time interval
  /// @param t: Drift time in nano seconds
  double reducedTime(const double t) const;
  /// @brief Maps the drift radius to the interval [-1;1.]
  ///        w.r.t configured min & max drift radius interval
  /// @param r: Drift radius of the measurement
  double reducedRadius(const double r) const;
  const Config m_cfg{};
  std::unique_ptr<const Acts::Logger> m_logger{};

  const Acts::Logger& logger() const { return *m_logger; }
};
}  // namespace ActsExamples
