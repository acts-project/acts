// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#pragma once

#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/EventData/MuonSpacePoint.hpp"

#include <memory>
namespace ActsExamples {

using namespace Acts::UnitLiterals;

/// @brief  Calibration class of M
class MuonSpacePointCalibrator {
 public:
  enum class CalibPolyType { Chebychev, Legendre };
  /// @brief Calibration configuration struct */
  struct Config {
    /// @brief Coefficients of the r-t relation parametrised as othogonal polynomial
    std::vector<double> rtCoefficients{9.0077,  6.7416,    -1.5805,
                                       0.56482, -0.1911,   0.018885,
                                       0.03027, -0.034393, 0.023613};
    ///  @brief Coefficients of the t-r relation as othogonal polynomial
    std::vector<double> trCoefficients{256.328, 348.773, 118.157, 18.667,
                                       -6.4899};
    /// @brief Coefficients of the drift radius uncertainty parameterized as an othogonal polynomial
    ///        as a function of the drift time
    std::vector<double> rtUncertCoefficients{0.10833, -0.071956, 0.037702,
                                             -0.011763};
    /// @brief Type of the calibration polynomial
    CalibPolyType rtPolyType{CalibPolyType::Chebychev};
    /// @brief Type of the calibration polynomial
    CalibPolyType trPolyType{CalibPolyType::Legendre};
    /// @brief Type of the calibration polynomial
    CalibPolyType rtUncertPolyType{CalibPolyType::Chebychev};
    /// @brief Minimum drift time covered by the r-t relation
    double minDriftT{0.};
    /// @brief Maximum drift time covered by the r-t relation
    double maxDriftT{740.203_ns};
    /// @brief Minimum radius mapped by the t-r relation
    double minTubeR{0.065987_mm};
    /// @brief Maximum radius mapped by the t-r relation
    double maxTubeR{14.567_mm};
    /// @brief Propagation velocity of the electronic signal along the wire
    double propagationV{0.};
    /// @brief Include trackT0 & time offset in the calibration
    bool includeTrackT0{false};
    /// @brief
    double rpcEtaStripPitch{2.6_cm};

    double rpcPhiStripPitch{3.5_cm};
  };

  /// @brief default constructor */
  MuonSpacePointCalibrator(const Config& cfg,
                           std::unique_ptr<const Acts::Logger> logger);
  /// @brief Abbrivation of the uncalibrated space point container type
  using UnCalibSpVec_t = std::vector<const MuonSpacePoint*>;
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

  /// @brief Calculates the drift radius from the drift time of a Mdt
  ///        measurement
  /// @param driftTime: Drift time in nanoseconds
  double driftRadius(const double driftTime) const;
  /// @brief Calculates the drift radius from the drift time
  /// @param driftTime: Drift time in nanoseconds
  double driftVelocity(const double driftTime) const;
  /// @brief Calculates the drift radius from the drift time
  /// @param driftTime: Drift time in nanoseconds
  double driftAcceleration(const double driftTime) const;
  /// @brief Calculate the drift time from the drift radius
  /// @param driftRadius: Drift radius in mm
  double driftTime(const double driftRadius) const;
  /// @brief Calculates the uncertainty on the drift radius measurement
  ///        as a function of the drift radius
  double driftRadiusUncert(const double driftRadius) const;
  /// @brief  Const access to the config
  const Config& config() const { return m_cfg; }

 private:
  /// @brief Evaluate the polynomial or its derivative which is either passed as a
  ///        series of Legendre or Cheby polynomials
  /// @param xValue: Unbounded value where the polynomial shall be evaluated
  /// @param derivative: Order of the derivative to evaluate. 0 is the bare
  ///                    polynomial itself
  /// @param polyType: Flag toglging whether the polynomial is Legendre or Cheby
  /// @param upperBound: Upper boundary of the passable xValues. If exceeded a nullopt
  ///                    is returned
  /// @param lowerBound: Lower boundary of the possible xValues. If exceeded a nullopt is returned.
  ///                    lowerBound & upperBound are used to map the xValue into
  ///                    the interval [-1;1]
  /// @param coeffs: List of the coefficient expansion
  std::optional<double> expandPolySeries(
      const double xValue, unsigned derivative, const CalibPolyType polyType,
      const double upperBound, const double lowerBound,
      const std::vector<double>& coeffs) const;

  const Acts::Logger& logger() const { return *m_logger; }

  std::unique_ptr<const Acts::Logger> m_logger{};
  Config m_cfg{};
  /// @brief Boundary values on the drift radius uncertainty
  double m_minDriftRUncert{driftRadius(m_cfg.minDriftT)};
  double m_maxDriftRUncert{driftRadius(m_cfg.maxDriftT)};
};
}  // namespace ActsExamples
