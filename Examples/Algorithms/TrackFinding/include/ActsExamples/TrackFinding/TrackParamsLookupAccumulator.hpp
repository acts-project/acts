// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/TrackFinding/TrackParamsLookupTable.hpp"

namespace ActsExamples {

/// @brief Class to accumulate and average track lookup tables
///
/// This class is used to accumulate track parameters in
/// reference layer grids and average them to create a lookup
/// table for track parameter estimation in seeding
class TrackParamsLookupAccumulator {
 public:
  /// @brief Nested configuration struct
  struct Config {
    /// Axis generator
    LookupAxisGen axisGen;
  };

  /// @brief Constructor
  TrackParamsLookupAccumulator(const Config& config)
      : m_cfg(std::move(config)),
        m_ipGrid(m_cfg.axisGen()),
        m_refGrid(m_cfg.axisGen()) {}

  /// @brief Add track parameters to the accumulator
  ///
  /// @param ipTrackParameters the track parameters at the IP
  /// @param refTrackParameters the track parameters at the reference layer
  /// @param position local position of the track hit on the reference layer
  void addTrack(const Acts::CurvilinearTrackParameters& ipTrackParameters,
                const Acts::CurvilinearTrackParameters& refTrackParameters,
                const Acts::Vector2& position);

  /// @brief Finalize the lookup table
  ///
  /// Return the grid with the bin track parameters averaged
  LookupGrid finalizeLookup();

  /// Get readonly access to the config parameters
  const Config& config() const { return m_cfg; }

 private:
  /// Configuration
  Config m_cfg;

  /// Mutex for modifying the grid
  std::mutex m_writeMutex;

  /// Grids to accumulate IP and reference
  /// layer track parameters
  LookupAccumGrid m_ipGrid;
  LookupAccumGrid m_refGrid;
};

}  // namespace ActsExamples
