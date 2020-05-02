// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <optional>
#include <utility>

#include "ACTFW/EventData/SimSourceLink.hpp"
#include "ACTFW/Validation/ProtoTrackClassification.hpp"
#include "Acts/EventData/MultiTrajectory.hpp"
#include "Acts/EventData/TrackParameters.hpp"

namespace FW {

/// @brief Struct for truth fitting/finding result
///
/// It contains a MultiTrajectory with a list of entry indices for different
/// tracks and fitted parameters for those tracks
struct SimMultiTrajectory {
 public:
  // Default constructor
  SimMultiTrajectory() = default;

  /// Constructor from fitted trajectory and fitted track parameter
  ///
  /// @param parameters The fitted track parameter
  /// @param trajectory The fitted multiTrajectory
  /// @param tTip The fitted multiTrajectory entry point
  SimMultiTrajectory(
      std::optional<Acts::BoundParameters> parameters,
      std::optional<Acts::MultiTrajectory<SimSourceLink>> trajectory,
      size_t tTip = SIZE_MAX)
      : m_trackTip(tTip) {
    if (parameters) {
      m_trackParameters = std::move(*parameters);
    }
    if (trajectory) {
      m_trajectory = std::move(*trajectory);
    }
  }

  /// @brief Copy constructor
  ///
  /// @param rhs is the source SimMultiTrajectory
  SimMultiTrajectory(const SimMultiTrajectory& rhs)
      : m_trackParameters(rhs.m_trackParameters),
        m_trajectory(rhs.m_trajectory),
        m_trackTip(rhs.m_trackTip) {}

  /// Copy move constructor
  ///
  /// @param rhs is the source SimMultiTrajectory
  SimMultiTrajectory(SimMultiTrajectory&& rhs)
      : m_trackParameters(std::move(rhs.m_trackParameters)),
        m_trajectory(std::move(rhs.m_trajectory)),
        m_trackTip(std::move(rhs.m_trackTip)) {}

  /// @brief Default destructor
  ///
  ~SimMultiTrajectory() = default;

  /// @brief assignment operator
  ///
  /// @param rhs is the source SimMultiTrajectory
  SimMultiTrajectory& operator=(const SimMultiTrajectory& rhs) {
    m_trackParameters = rhs.m_trackParameters;
    m_trajectory = rhs.m_trajectory;
    m_trackTip = rhs.m_trackTip;
    return *this;
  }

  /// @brief assignment move operator
  ///
  /// @param rhs is the source SimMultiTrajectory
  SimMultiTrajectory& operator=(SimMultiTrajectory&& rhs) {
    m_trackParameters = std::move(rhs.m_trackParameters);
    m_trajectory = std::move(rhs.m_trajectory);
    m_trackTip = std::move(rhs.m_trackTip);
    return *this;
  }

  /// Get fitted track parameter
  const Acts::BoundParameters& trackParameters() const {
    if (m_trackParameters) {
      return *m_trackParameters;
    } else {
      throw std::runtime_error(
          "No fitted track parameter for this trajectory!");
    }
  }

  /// Get trajectory along with the entry point
  const std::pair<size_t, Acts::MultiTrajectory<SimSourceLink>> trajectory()
      const {
    if (m_trajectory) {
      return std::make_pair(m_trackTip, *m_trajectory);
    } else {
      throw std::runtime_error("No fitted states on this trajectory!");
    };
  }

  /// Get number of track states
  size_t numStates() const {
    size_t nStates = 0;
    if (m_trajectory) {
      (*m_trajectory).visitBackwards(m_trackTip, [&](const auto&) {
        nStates++;
      });
    }
    return nStates;
  }

  /// Get number of track states that have measurements
  size_t numMeasurements() const {
    size_t nMeasurements = 0;
    if (m_trajectory) {
      (*m_trajectory).visitBackwards(m_trackTip, [&](const auto& state) {
        if (state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
          nMeasurements++;
        }
      });
    }
    return nMeasurements;
  }

  /// Indicator for having fitted track parameter or not
  bool hasTrackParameters() const { return m_trackParameters != std::nullopt; }

  /// Indicator for having fitted trajectory or not
  bool hasTrajectory() const { return m_trajectory != std::nullopt; }

  /// Get the truth particle counts to help identify majority particle
  std::vector<ParticleHitCount> identifyMajorityParticle() const {
    std::vector<ParticleHitCount> particleHitCount;

    if (m_trajectory) {
      (*m_trajectory).visitBackwards(m_trackTip, [&](const auto& state) {
        // No truth info with non-measurement state
        if (not state.typeFlags().test(Acts::TrackStateFlag::MeasurementFlag)) {
          return true;
        }
        // Find the truth particle associated with this state
        const auto particleId = state.uncalibrated().truthHit().particleId();
        // Find if the particle already exists
        auto it = std::find_if(particleHitCount.begin(), particleHitCount.end(),
                               [=](const ParticleHitCount& phc) {
                                 return phc.particleId == particleId;
                               });

        // Either increase count if we saw the particle before or add it
        if (it != particleHitCount.end()) {
          it->hitCount += 1;
        } else {
          particleHitCount.push_back({particleId, 1u});
        }
        return true;
      });
    }

    if (not particleHitCount.empty()) {
      // sort by hit count, i.e. majority particle first
      std::sort(particleHitCount.begin(), particleHitCount.end(),
                [](const ParticleHitCount& lhs, const ParticleHitCount& rhs) {
                  return lhs.hitCount > rhs.hitCount;
                });
    }

    return particleHitCount;
  }

 private:
  // The optional Parameters at the provided surface
  std::optional<Acts::BoundParameters> m_trackParameters{std::nullopt};

  // The optional fitted multitrajectory
  std::optional<Acts::MultiTrajectory<SimSourceLink>> m_trajectory{
      std::nullopt};

  // This is the index of the 'tip' of the track stored in multitrajectory.
  size_t m_trackTip = SIZE_MAX;
};

}  // namespace FW
