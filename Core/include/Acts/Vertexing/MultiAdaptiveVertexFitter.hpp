// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Vertexing/Chi2TrackCompatibilityEstimator.hpp"
#include "Acts/Vertexing/ImpactPoint3dEstimator.hpp"
#include "Acts/Vertexing/KalmanVertexUpdator.hpp"
#include "Acts/Vertexing/LinearizedTrackFactory.hpp"
#include "Acts/Vertexing/MAVFInfo.hpp"
#include "Acts/Vertexing/SequentialVertexSmoother.hpp"
#include "Acts/Vertexing/VertexAnnealingTool.hpp"
#include "Acts/Vertexing/VertexFitterOptions.hpp"

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Result.hpp"

#include <functional>

namespace Acts {

// TODO: add docs
template <typename bfield_t, typename input_track_t, typename propagator_t>
class MultiAdaptiveVertexFitter {
 public:
  struct State {
    // Vertex collection to be fitted
    std::vector<Vertex<input_track_t>> vertexCollection;

    // Annealing state
    VertexAnnealingTool::State annealingState;
  };

  struct Config {
    Config(const bfield_t& bIn)
        : linFactory(
              typename LinearizedTrackFactory<bfield_t, propagator_t>::Config(
                  bIn)) {}

    /// Linearized track factory
    LinearizedTrackFactory<bfield_t, propagator_t> linFactory;

    /// Track compatibility estimator
    Chi2TrackCompatibilityEstimator<input_track_t> trackCompEst;

    /// ImpactPoint3dEstimator
    ImpactPoint3dEstimator<input_track_t> ipEst;

    /// Vertex updator
    KalmanVertexUpdator<input_track_t> vertexUpdator;

    /// SequentialVertexSmoother
    SequentialVertexSmoother<input_track_t> vertexSmoother;
    /// Annealing tool
    VertexAnnealingTool annealingTool;

    /// Number of max iterations
    int maxIterations{30};

    /// Max distance to linearization point allowed
    /// without relinearization
    double maxDistToLinPoint{0.5};
  };

  /// @brief Constructor used if input_track_t type == BoundParameters
  ///
  /// @param cfg Configuration object
  template <typename T = input_track_t,
            std::enable_if_t<std::is_same<T, BoundParameters>::value, int> = 0>
  MultiAdaptiveVertexFitter(const Config& config)
      : m_cfg(config), m_extractParameters([](T params) { return params; }) {}

  /// @brief Constructor for user-defined input_track_t type =! BoundParameters
  ///
  /// @param cfg Configuration object
  /// @param func Function extracting BoundParameters from input_track_t object
  MultiAdaptiveVertexFitter(const Config& config,
                            std::function<BoundParameters(input_track_t)> func)
      : m_cfg(config), m_extractParameters(func) {}

  /// @brief The actual fit function
  ///
  /// @param state The state object
  /// @param vFitterOptions Vertex fitter options
  ///
  /// @return Result<void> object
  Result<void> fit(
      State& state,
      const VertexFitterOptions<input_track_t>& vFitterOptions) const;

  /// @brief Adds new vertex to a previous multi-vertex fit
  /// and fits everything together:
  /// 1. The new vertex is added to the fit (all the tracks get initialized,
  /// so that the plane through their IP point and the seed vertex
  /// (IP3dAtAPlane) is created, to be later able to estimate in a fast way the
  /// compatibility of the tracks to their respective vertices.
  /// 2. All tracks belonging to the new vertex are scanned and all the vertices
  ///  which shares tracks with the new vertex to be fit are also added to the
  ///  fit.
  /// 3. The multivertex fit is performed with all involved vertices.
  ///
  /// This has the advantage that only vertices that are affected by adding the
  /// new vertex get refitted.
  ///
  /// Note: newVertex has to be properly initialized (seed vertex,
  /// constraint vertex, list of MAV)
  ///
  /// @param state The state object
  /// @param newVertex New vertex to be added to fit
  /// @param MAVFTrackAtVtxInfo TrackAtVertex info object
  /// @param vFitterOptions Vertex fitter options
  ///
  /// @return Result<void> object
  Result<void> addVertexToFit(
      State& state, Vertex<input_track_t> newVertex,
      MAVFTrackAtVtxInfo<input_track_t> trackAtVtxInfo,
      const VertexFitterOptions<input_track_t>& vFitterOptions) const;

 private:
  /// Configuration object
  const Config m_cfg;

  /// @brief Function to extract track parameters,
  /// input_track_t objects are BoundParameters by default, function to be
  /// overwritten to return BoundParameters for other input_track_t objects.
  ///
  /// @param input_track_t object to extract track parameters from
  const std::function<BoundParameters(input_track_t)> m_extractParameters;

  /// @brief Tests if vertex is already in list of vertices or not
  ///
  /// @param vtx Vertex to test
  /// @param verticesVec Vector of vertices to search
  ///
  /// @return True if vtx is already in verticesVec
  bool isAlreadyInList(
      Vertex<input_track_t>* vtx,
      const std::vector<Vertex<input_track_t>*>& verticesVec) const;
};

}  // namespace Acts

#include "Acts/Vertexing/MultiAdaptiveVertexFitter.ipp"
