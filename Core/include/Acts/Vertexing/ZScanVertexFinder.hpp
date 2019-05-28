// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <unordered_map>

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"

#include "Acts/Vertexing/FsmwMode1dFinder.hpp"
#include "Acts/Vertexing/TrackToVertexIPEstimator.hpp"
#include "Acts/Vertexing/Vertex.hpp"

#include "Acts/Vertexing/VertexFinderOptions.hpp"

namespace Acts {

/// @class ZScanVertexFinder
///
/// @brief Implements a vertex finder based on the mode of z0 values:
/// 1. Determines the mode value of all input track z0 values
/// 2. If no contraint is given, returns (0,0, z0_mode) as vertex position
/// 3. If vertex contraint is given with x=x_constr and y=y_constr,
///    the returned vertex position will be (x_constr, y_constr, z0_mode).
template <typename bfield_t, typename input_track_t, typename propagator_t>
class ZScanVertexFinder {
 private:
  // functor to compare two unordered_map Key values for equality
  struct pred_perigee {
    bool operator()(const BoundParameters& left,
                    const BoundParameters& right) const {
      return (left.parameters()[ParID_t::eLOC_D0] ==
              right.parameters()[ParID_t::eLOC_D0]) &&
             (left.parameters()[ParID_t::eLOC_Z0] ==
              right.parameters()[ParID_t::eLOC_Z0]) &&
             (left.parameters()[ParID_t::ePHI] ==
              right.parameters()[ParID_t::ePHI]) &&
             (left.parameters()[ParID_t::eTHETA] ==
              right.parameters()[ParID_t::eTHETA]) &&
             (left.parameters()[ParID_t::eQOP] ==
              right.parameters()[ParID_t::eQOP]);
    }
  };

  // functor to hash key for unordered_map
  struct hash_perigee {
    size_t operator()(const BoundParameters& perigee) const {
      return std::hash<double>()(perigee.parameters()[ParID_t::eLOC_D0]) ^
             std::hash<double>()(perigee.parameters()[ParID_t::eLOC_Z0]) ^
             std::hash<double>()(perigee.parameters()[ParID_t::ePHI]) ^
             std::hash<double>()(perigee.parameters()[ParID_t::eTHETA]) ^
             std::hash<double>()(perigee.parameters()[ParID_t::eQOP]);
    }
  };

 public:
  /// @struct Config Configuration struct
  struct Config {
    // TrackToVertexIPEstimator
    TrackToVertexIPEstimator<input_track_t, propagator_t> ipEstimator;

    // FsmwMode1dFinder
    FsmwMode1dFinder mode1dFinder;

    // disables all weights, set all weights to 1.
    bool disableAllWeights = false;
    // constraint parameters
    float constraintcutoff = 9.;
    float constrainttemp = 1.;
    // use LogPt for weighting
    bool useLogPt = true;
    // use pt for weighting
    bool usePt = false;
    // minimum pt
    double minPt = 0.4 * units::_GeV;
    // exponent used for weighting if usePt
    double expPt = 1.;
    // cache weights
    bool cacheWeights = true;
    // minimum required weight
    double minWeight = 0.01;
  };

  struct State {
    // Empty vertex collection, to be filled by finder
    // ZScanVertexFinder always returns only one single seed,
    // hence this collection will always only be filled with
    // a single vertex candidate
    std::vector<Vertex<input_track_t>> vertexCollection;

    // hashtable to avoid computing perigee more than once per track
    std::unordered_map<BoundParameters, std::pair<double, double>, hash_perigee,
                       pred_perigee>
        weightMap;
  };

  /// @brief Constructor used if input_track_t type == BoundParameters
  ///
  /// @param cfg Configuration object
  /// @param logger The logging instance
  template <typename T = input_track_t,
            std::enable_if_t<std::is_same<T, BoundParameters>::value, int> = 0>
  ZScanVertexFinder(Config cfg,
                    std::unique_ptr<const Logger> logger =
                        getDefaultLogger("ZScanVertexFinder", Logging::INFO))
      : m_cfg(std::move(cfg)),
        m_extractParameters([](T params) { return params; }),
        m_logger(std::move(logger)) {}

  /// @brief Constructor for user-defined input_track_t type =! BoundParameters
  ///
  /// @param cfg Configuration object
  /// @param func Function extracting BoundParameters from input_track_t object
  /// @param logger Logging instance
  ZScanVertexFinder(Config cfg,
                    std::function<BoundParameters(input_track_t)> func,
                    std::unique_ptr<const Logger> logger =
                        getDefaultLogger("ZScanVertexFinder", Logging::INFO))
      : m_cfg(std::move(cfg)),
        m_extractParameters(func),
        m_logger(std::move(logger)) {}

  /// @brief Function that determines single vertex,
  /// based on z0 values of input tracks,
  /// using a Half Sample Mode algorithm
  ///
  /// @param trackVector Input track collection
  /// @param state The state object, to be filled with the determined vertex in
  /// vertexCollection
  /// @param propagator Propagator
  /// @param vFinderOptions Vertex finder options
  ///
  /// @return Result void object, empty in case of success,
  /// otherwise holding error
  Result<void> find(
      const std::vector<input_track_t>& trackVector, State& state,
      const propagator_t& propagator,
      const VertexFinderOptions<input_track_t>& vFinderOptions) const;

 private:
  Config m_cfg;

  /// @brief Function to extract track parameters,
  /// input_track_t objects are BoundParameters by default, function to be
  /// overwritten to return BoundParameters for other input_track_t objects.
  ///
  /// @param input_track_t object to extract track parameters from
  std::function<BoundParameters(input_track_t)> m_extractParameters;

  /// Logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private access to logging instance
  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts

#include "ZScanVertexFinder.ipp"