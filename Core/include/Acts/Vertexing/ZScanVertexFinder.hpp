// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/FsmwMode1dFinder.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexFitterConcept.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"

#include <unordered_map>

namespace Acts {

/// @class ZScanVertexFinder
///
/// @brief Implements a vertex finder based on the mode of z0 values:
/// 1. Determines the mode value of all input track z0 values
/// 2. If no constraint is given, returns (0,0, z0_mode) as vertex position
/// 3. If vertex constraint is given with x=x_constr and y=y_constr,
///    the returned vertex position will be (x_constr, y_constr, z0_mode).
template <typename vfitter_t>
class ZScanVertexFinder {
  static_assert(VertexFitterConcept<vfitter_t>,
                "Vertex fitter does not fulfill vertex fitter concept.");
  using Propagator_t = typename vfitter_t::Propagator_t;

 public:
  using InputTrack_t = typename vfitter_t::InputTrack_t;

  /// Configuration struct
  struct Config {
    /// @brief Finder configuration
    ///
    /// @param ipEst ImpactPointEstimator
    Config(const ImpactPointEstimator<InputTrack_t, Propagator_t>& ipEst)
        : ipEstimator(ipEst) {}

    // ImpactPointEstimator
    ImpactPointEstimator<InputTrack_t, Propagator_t> ipEstimator;

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
    double minPt = 0.4 * UnitConstants::GeV;
    // exponent used for weighting if usePt
    double expPt = 1.;
    // minimum required weight
    double minWeight = 0.01;
  };

  /// State struct for fulfilling interface
  struct State {};

  /// @brief Constructor used if InputTrack_t type == BoundTrackParameters
  ///
  /// @param cfg Configuration object
  /// @param logger The logging instance
  template <
      typename T = InputTrack_t,
      std::enable_if_t<std::is_same<T, BoundTrackParameters>::value, int> = 0>

  ZScanVertexFinder(const Config& cfg,
                    std::unique_ptr<const Logger> logger =
                        getDefaultLogger("ZScanVertexFinder", Logging::INFO))
      : m_cfg(std::move(cfg)),
        m_extractParameters([](T params) { return params; }),
        m_logger(std::move(logger)) {}

  /// @brief Constructor for user-defined InputTrack_t type =!
  /// BoundTrackParameters
  ///
  /// @param cfg Configuration object
  /// @param func Function extracting BoundTrackParameters from InputTrack_t
  /// object
  /// @param logger Logging instance
  ZScanVertexFinder(const Config& cfg,
                    std::function<BoundTrackParameters(InputTrack_t)> func,
                    std::unique_ptr<const Logger> logger =
                        getDefaultLogger("ZScanVertexFinder", Logging::INFO))
      : m_cfg(cfg), m_extractParameters(func), m_logger(std::move(logger)) {}

  /// @brief Function that determines single vertex,
  /// based on z0 values of input tracks,
  /// using a Half Sample Mode algorithm
  ///
  /// @param trackVector Input track collection
  /// @param vertexingOptions Vertexing options
  /// @param state State for fulfilling correct interface
  ///
  /// @return Vector of vertices, filled with a single
  ///         vertex (for consistent interfaces)
  Result<std::vector<Vertex<InputTrack_t>>> find(
      const std::vector<const InputTrack_t*>& trackVector,
      const VertexingOptions<InputTrack_t>& vertexingOptions,
      State& state) const;

 private:
  Config m_cfg;

  /// @brief Function to extract track parameters,
  /// InputTrack_t objects are BoundTrackParameters by default, function to be
  /// overwritten to return BoundTrackParameters for other InputTrack_t objects.
  std::function<BoundTrackParameters(InputTrack_t)> m_extractParameters;

  /// Logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private access to logging instance
  const Logger& logger() const { return *m_logger; }
};

}  // namespace Acts

#include "ZScanVertexFinder.ipp"
