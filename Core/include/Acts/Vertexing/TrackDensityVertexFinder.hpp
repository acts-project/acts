// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexFinderOptions.hpp"

namespace Acts {

/// @class TrackDensityVertexFinder
///
/// TODO: docs
template <typename vfitter_t>
class TrackDensityVertexFinder {


  /// @struct Config Configuration struct
  struct Config {
    
  };


  /// @brief TODO
  Result<std::vector<Vertex<InputTrack_t>>> find(
      const std::vector<InputTrack_t>& trackVector,
      const VertexFinderOptions<InputTrack_t>& vFinderOptions) const;

 private:
  Config m_cfg;

  /// @brief Function to extract track parameters,
  /// InputTrack_t objects are BoundParameters by default, function to be
  /// overwritten to return BoundParameters for other InputTrack_t objects.
  ///
  /// @param InputTrack_t object to extract track parameters from
  std::function<BoundParameters(InputTrack_t)> m_extractParameters;

  /// Logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private access to logging instance
  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts

#include "TrackDensityVertexFinder.ipp"
