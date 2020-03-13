// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/AnnealingUtility.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/AMVFInfo.hpp"
#include "Acts/Vertexing/ImpactPoint3dEstimator.hpp"
#include "Acts/Vertexing/LinearizerConcept.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexFitterOptions.hpp"

#include <functional>

namespace Acts {

/// @class AdaptiveMultiVertexFitter
/// @brief Implements an adaptive multi-vertex fitter as described
///   in detail in Section 5.3.5 in:
///   Ref. (1): CERN-THESIS-2010-027, Author: Piacquadio, Giacinto:
///   `Identification of b-jets and investigation of the discovery potential
///   of a Higgs boson in the WH−−>lvbb¯ channel with the ATLAS experiment`
///
///////////////////////////////////////////////////////////////////////////
///
/// @tparam input_track_t Track object type
/// @tparam linearizer_t Track linearizer type
template <typename input_track_t, typename linearizer_t>
class AdaptiveMultiVertexFitter {
  static_assert(LinearizerConcept<linearizer_t>,
                "Linearizer does not fulfill linearizer concept.");

 public:
  using InputTrack_t = input_track_t;
  using Propagator_t = typename linearizer_t::Propagator_t;
  using Linearizer_t = linearizer_t;

 private:
  using ImpactPointEstimator =
      ImpactPoint3dEstimator<InputTrack_t, Propagator_t>;

 public:
  /// @brief The fitter state
  struct State {
    // Vertex collection to be fitted
    std::vector<Vertex<InputTrack_t>*> vertexCollection;

    // Annealing state
    AnnealingUtility::State annealingState;

    // Map to store vertices information
    std::map<Vertex<InputTrack_t>*, VertexInfo<InputTrack_t>> vtxInfoMap;

    std::multimap<const InputTrack_t*, Vertex<InputTrack_t>*>
        trackToVerticesMultiMap;

    std::map<std::pair<const InputTrack_t*, Vertex<InputTrack_t>*>,
             TrackAtVertex<InputTrack_t>>
        tracksAtVerticesMap;

    /// @brief Default State constructor
    State() = default;

    // Adds a vertex to trackToVerticesMultiMap
    void addVertexToMultiMap(Vertex<InputTrack_t>& vtx) {
      for (auto trk : vtxInfoMap[&vtx].trackLinks) {
        trackToVerticesMultiMap.emplace(trk, &vtx);
      }
    }

    // Removes a vertex from trackToVerticesMultiMap
    void removeVertexFromMultiMap(Vertex<InputTrack_t>& vtx) {
      for (auto iter = trackToVerticesMultiMap.begin();
           iter != trackToVerticesMultiMap.end();) {
        if (iter->second == &vtx) {
          iter = trackToVerticesMultiMap.erase(iter);
        } else {
          ++iter;
        }
      }
    }
  };

  struct Config {
    /// @brief Config constructor
    ///
    /// @param est ImpactPoint3dEstimator
    Config(ImpactPointEstimator est) : ipEst(std::move(est)) {}

    // ImpactPoint3dEstimator
    ImpactPointEstimator ipEst;

    /// Annealing tool used for a thermodynamic annealing scheme for the
    /// track weight factors in such a way that with high temperature values
    /// (at the beginning) only a slight preference is given to tracks
    /// compatible with the estimated vertex position. With lower temperatures
    /// the weighting get stricter such that all incompatible tracks will be
    /// dropped at the end while keeping all compatible tracks with a weight=1.
    ///   Ref. (1): CERN-THESIS-2010-027, Author: Piacquadio, Giacinto:
    ///   `Identification of b-jets and investigation of the discovery potential
    ///   of a Higgs boson in the WH−−>lvbb¯ channel with the ATLAS experiment`
    AnnealingUtility annealingTool;

    // Number of max iterations
    unsigned int maxIterations{30};

    // Max distance to linearization point allowed
    // without relinearization
    double maxDistToLinPoint{0.5};

    // Minimum track weight needed for track to be considered
    double minWeight{0.0001};

    // Max relative shift of vertex during one iteration
    double maxRelativeShift{0.01};

    // Do smoothing after multivertex fit
    bool doSmoothing{false};
  };

  /// @brief Constructor used if InputTrack_t type == BoundParameters
  ///
  /// @param cfg Configuration object
  /// @param logger The logging instance
  template <typename T = InputTrack_t,
            std::enable_if_t<std::is_same<T, BoundParameters>::value, int> = 0>
  AdaptiveMultiVertexFitter(Config& cfg,
                            std::unique_ptr<const Logger> logger =
                                getDefaultLogger("AdaptiveMultiVertexFitter",
                                                 Logging::INFO))
      : m_cfg(std::move(cfg)),
        m_extractParameters([](T params) { return params; }),
        m_logger(std::move(logger)) {}

  /// @brief Constructor for user-defined InputTrack_t type != BoundParameters
  ///
  /// @param cfg Configuration object
  /// @param func Function extracting BoundParameters from InputTrack_t object
  /// @param logger The logging instance
  AdaptiveMultiVertexFitter(Config& cfg,
                            std::function<BoundParameters(InputTrack_t)> func,
                            std::unique_ptr<const Logger> logger =
                                getDefaultLogger("AdaptiveMultiVertexFitter",
                                                 Logging::INFO))
      : m_cfg(std::move(cfg)),
        m_extractParameters(func),
        m_logger(std::move(logger)) {}

  /// @brief The actual fit function, performs a simultaneous
  /// fit of all vertices in `verticesToFit` by invoking `fitImpl`
  ///
  /// @param state The state object
  /// @param verticesToFit Vector containing all vertices to be fitted
  /// @param linearizer The track linearizer
  /// @param vFitterOptions Vertex fitter options
  ///
  /// @return Result<void> object
  Result<void> fit(
      State& state, const std::vector<Vertex<InputTrack_t>*>& verticesToFit,
      const Linearizer_t& linearizer,
      const VertexFitterOptions<InputTrack_t>& vFitterOptions) const;

  /// @brief Adds new vertex to an existing multi-vertex fit
  /// and fits everything together (by invoking the fit_impl method):
  /// 1. The new vertex is added to the fit: all associated tracks get
  /// initialized, i.e. ParamsAtIP3d are created (from ImpactPoint3dEstimator)
  /// to be later able to estimate in a fast way the compatibility of the tracks
  /// to their respective vertices.
  /// 2. All tracks belonging to the new vertex are scanned and all the vertices
  ///  which share tracks with the new vertex to be fit are also added to the
  ///  fit.
  /// 3. The multivertex fit is performed with all involved vertices.
  ///
  /// This has the advantage that only vertices that are affected by adding the
  /// new vertex are refitted.
  ///
  /// Note: newVertex has to be properly initialized (seed vertex,
  /// constraint vertex, list of MAV)
  ///
  /// @param state The state object
  /// @param newVertex New vertex to be added to fit
  /// @param linearizer The track linearizer
  /// @param vFitterOptions Vertex fitter options
  ///
  /// @return Result<void> object
  Result<void> addVtxToFit(
      State& state, Vertex<InputTrack_t>& newVertex,
      const Linearizer_t& linearizer,
      const VertexFitterOptions<InputTrack_t>& vFitterOptions) const;

 private:
  /// Configuration object
  const Config m_cfg;

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

  /// @brief The actual fit function, performs a simultaneous
  ///   fit of all vertices in state.vertexCollection
  ///
  /// @param state The state object
  /// @param linearizer The track linearizer
  /// @param vFitterOptions Vertex fitter options
  ///
  /// @return Result<void> object
  Result<void> fitImpl(
      State& state, const Linearizer_t& linearizer,
      const VertexFitterOptions<InputTrack_t>& vFitterOptions) const;

  /// @brief Tests if vertex is already in list of vertices or not
  ///
  /// @param vtx Vertex to test
  /// @param verticesVec Vector of vertices to search
  ///
  /// @return True if vtx is already in verticesVec
  bool isAlreadyInList(
      Vertex<InputTrack_t>* vtx,
      const std::vector<Vertex<InputTrack_t>*>& verticesVec) const;

  /// @brief Prepares vertex object for the actual fit, i.e.
  /// all TrackAtVertex objects at current vertex will obtain
  /// `ip3dParams` from ImpactPoint3dEstimator::getParamsAtClosestApproach
  /// in order to later faster estimate compatibilities of track
  /// with different vertices
  ///
  /// @param vtx The vertex object
  /// @param vFitterOptions Vertex fitter options
  Result<void> prepareVertexForFit(
      State& state, Vertex<InputTrack_t>* vtx,
      const VertexFitterOptions<InputTrack_t>& vFitterOptions) const;

  /// @brief Sets vertexCompatibility for all TrackAtVertex objects
  /// at current vertex
  ///
  /// @param state The state object
  /// @param geoContext The geometry context
  /// @param currentVtx Current vertex
  Result<void> setAllVertexCompatibilities(
      State& state, const GeometryContext& geoContext,
      Vertex<InputTrack_t>* currentVtx) const;

  /// @brief Sets weights to the track according to Eq.(5.46) in Ref.(1)
  ///  and updates the vertices by calling the VertexUpdater
  ///
  /// @param state The state object
  /// @param linearizer The track linearizer
  Result<void> setWeightsAndUpdate(State& state,
                                   const Linearizer_t& linearizer) const;

  /// @brief Collects all compatibility values of the track `trk`
  /// at all vertices it is currently attached to and outputs
  /// these values in a vector
  ///
  /// @param state The state object
  /// @param trk The track
  ///
  /// @return Vector of compatibility values
  std::vector<double> collectTrackToVertexCompatibilities(
      State& state, const InputTrack_t* trk) const;

  /// @brief Determines if vertex position has shifted more than
  /// m_cfg.maxRelativeShift in last iteration
  ///
  /// @param state The state object
  ///
  /// @return False if shift was larger than maxRelativeShift
  bool checkSmallShift(State& state) const;

  /// @brief Updates tracks for current vertex with knowledge
  /// of current vertex position
  ///
  /// @param state The state object
  /// @param geoContext The geometry context
  void doVertexSmoothing(State& state, const GeometryContext& geoContext) const;
};

}  // namespace Acts

#include "Acts/Vertexing/AdaptiveMultiVertexFitter.ipp"
