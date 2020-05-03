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
using IndexedParams = std::unordered_map<size_t, Acts::BoundParameters>;

/// @brief Struct for truth track fitting/finding result with
/// Acts::KalmanFitter/Acts::CombinatorialKalmanFilter
///
/// It contains a MultiTrajectory with a vector of entry indices for individual
/// trajectories, and a map of fitted parameters indexed by the entry index.
/// In case of track fitting, there is at most one trajectory in the
/// MultiTrajectory; In case of track finding, there could be multiple
/// trajectories in the MultiTrajectory.
struct SimMultiTrajectory {
 public:
  // Default constructor
  SimMultiTrajectory() = default;

  /// Constructor from multiTrajectory and fitted track parameters
  ///
  /// @param multiTraj The multiTrajectory
  /// @param tTips The entry indices for trajectories in multiTrajectory
  /// @param parameters The fitted track parameters indexed by trajectory entry
  /// index
  SimMultiTrajectory(
      std::optional<Acts::MultiTrajectory<SimSourceLink>> multiTraj,
      const std::vector<size_t>& tTips, const IndexedParams& parameters)
      : m_trackTips(tTips), m_trackParameters(parameters) {
    if (multiTraj) {
      m_multiTrajectory = std::move(*multiTraj);
    }
  }

  /// @brief Copy constructor
  ///
  /// @param rhs The source SimMultiTrajectory
  SimMultiTrajectory(const SimMultiTrajectory& rhs)
      : m_multiTrajectory(rhs.m_multiTrajectory),
        m_trackTips(rhs.m_trackTips),
        m_trackParameters(rhs.m_trackParameters) {}

  /// Copy move constructor
  ///
  /// @param rhs The source SimMultiTrajectory
  SimMultiTrajectory(SimMultiTrajectory&& rhs)
      : m_multiTrajectory(std::move(rhs.m_multiTrajectory)),
        m_trackTips(std::move(rhs.m_trackTips)),
        m_trackParameters(std::move(rhs.m_trackParameters)) {}

  /// @brief Default destructor
  ///
  ~SimMultiTrajectory() = default;

  /// @brief assignment operator
  ///
  /// @param rhs The source SimMultiTrajectory
  SimMultiTrajectory& operator=(const SimMultiTrajectory& rhs) {
    m_multiTrajectory = rhs.m_multiTrajectory;
    m_trackTips = rhs.m_trackTips;
    m_trackParameters = rhs.m_trackParameters;
    return *this;
  }

  /// @brief assignment move operator
  ///
  /// @param rhs The source SimMultiTrajectory
  SimMultiTrajectory& operator=(SimMultiTrajectory&& rhs) {
    m_multiTrajectory = std::move(rhs.m_multiTrajectory);
    m_trackTips = std::move(rhs.m_trackTips);
    m_trackParameters = std::move(rhs.m_trackParameters);
    return *this;
  }

  /// @brief Indicator of multiTrajectory
  ///
  /// @return Whether having multiTrajectory or not
  bool hasTrajectory() const { return m_multiTrajectory != std::nullopt; }

  /// @brief Indicator of fitted track parameters for one trajectory
  ///
  /// @param entryIndex The trajectory entry index
  ///
  /// @return Whether having fitted track parameters or not
  bool hasTrackParameters(const size_t& entryIndex) const {
    return m_trackParameters.count(entryIndex) > 0;
  }

  /// @brief Getter for multiTrajectory
  ///
  /// @return The multiTrajectory with trajectory entry indices
  std::pair<std::vector<size_t>, Acts::MultiTrajectory<SimSourceLink>>
  trajectory() const {
    if (m_multiTrajectory and not m_trackTips.empty()) {
      return std::make_pair(m_trackTips, *m_multiTrajectory);
    } else {
      throw std::runtime_error("No multiTrajectory available!");
    };
  }

  /// @brief Getter of fitted track parameters for one trajectory
  ///
  /// @param entryIndex The trajectory entry index
  ///
  /// @return The fitted track parameters of the trajectory
  const Acts::BoundParameters& trackParameters(const size_t& entryIndex) const {
    auto it = m_trackParameters.find(entryIndex);
    if (it != m_trackParameters.end()) {
      return it->second;
    } else {
      throw std::runtime_error(
          "No fitted track parameters for trajectory with entry index = " +
          std::to_string(entryIndex));
    }
  }

  /// @brief Counter of associated truth particles for one trajectory
  ///
  /// @param entryIndex The trajectory entry index
  ///
  /// @return The truth particle counts in ascending order
  std::vector<ParticleHitCount> identifyMajorityParticle(
      const size_t& entryIndex) const {
    std::vector<ParticleHitCount> particleHitCount;
    particleHitCount.reserve(10);
    if (m_multiTrajectory) {
      (*m_multiTrajectory).visitBackwards(entryIndex, [&](const auto& state) {
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
  // The optional fitted multiTrajectory
  std::optional<Acts::MultiTrajectory<SimSourceLink>> m_multiTrajectory{
      std::nullopt};

  // The entry indices of trajectories stored in multiTrajectory
  std::vector<size_t> m_trackTips = {};

  // The optional Parameters at the provided surface
  IndexedParams m_trackParameters = {};
};

}  // namespace FW
