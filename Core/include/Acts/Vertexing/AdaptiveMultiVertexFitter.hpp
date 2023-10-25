// This file is part of the Acts project.
//
// Copyright (C) 2019-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/AnnealingUtility.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/AMVFInfo.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/LinearizerConcept.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"

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
  using IPEstimator = ImpactPointEstimator<InputTrack_t, Propagator_t>;

 public:
  /// @brief The fitter state
  struct State {
    State(const MagneticFieldProvider& field,
          const Acts::MagneticFieldContext& magContext)
        : ipState(field.makeCache(magContext)),
          linearizerState(field.makeCache(magContext)) {}
    // Vertex collection to be fitted
    std::vector<Vertex<InputTrack_t>*> vertexCollection;

    // Annealing state
    AnnealingUtility::State annealingState;

    // IPEstimator state
    typename IPEstimator::State ipState;

    // Linearizer state
    typename Linearizer_t::State linearizerState;

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
    /// @param est ImpactPointEstimator
    Config(const IPEstimator& est) : ipEst(est) {}

    // ImpactPointEstimator
    IPEstimator ipEst;

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

    // Use time information when calculating the vertex compatibility
    bool useTime{false};
  };

  /// @brief Constructor used if InputTrack_t type == BoundTrackParameters
  ///
  /// @param cfg Configuration object
  /// @param logger The logging instance
  template <
      typename T = InputTrack_t,
      std::enable_if_t<std::is_same<T, BoundTrackParameters>::value, int> = 0>
  AdaptiveMultiVertexFitter(Config& cfg,
                            std::unique_ptr<const Logger> logger =
                                getDefaultLogger("AdaptiveMultiVertexFitter",
                                                 Logging::INFO))
      : m_cfg(std::move(cfg)),
        m_extractParameters([](T params) { return params; }),
        m_logger(std::move(logger)) {}

  /// @brief Constructor for user-defined InputTrack_t type !=
  /// BoundTrackParameters
  ///
  /// @param cfg Configuration object
  /// @param func Function extracting BoundTrackParameters from InputTrack_t
  /// object
  /// @param logger The logging instance
  AdaptiveMultiVertexFitter(
      Config& cfg, std::function<BoundTrackParameters(InputTrack_t)> func,
      std::unique_ptr<const Logger> logger =
          getDefaultLogger("AdaptiveMultiVertexFitter", Logging::INFO))
      : m_cfg(std::move(cfg)),
        m_extractParameters(func),
        m_logger(std::move(logger)) {}

  /// @brief Performs a simultaneous fit of all vertices in `verticesToFit`
  /// by invoking `fitImpl`.
  ///
  /// @param state Fitter state
  /// @param verticesToFit Vector containing all vertices to be fitted
  /// @param linearizer Track linearizer
  /// @param vertexingOptions Vertexing options
  ///
  /// @return Result<void> object
  Result<void> fit(
      State& state, const std::vector<Vertex<InputTrack_t>*>& verticesToFit,
      const Linearizer_t& linearizer,
      const VertexingOptions<InputTrack_t>& vertexingOptions) const;

  /// @brief Adds a new vertex to an existing multi-vertex fit.
  /// 1. The 3D impact parameters are calculated for all tracks associated
  /// with newVertex.
  /// 2. A list of all vertices that are connected with newVertex via shared
  /// tracks is created. It also considers indirect connections (e.g., vertex A
  /// which shares a track with vertex B which, in turn, shares a track with
  /// newVertex).
  /// 3. The multivertex fit is performed for all vertices on said list.
  ///
  /// @param state Fitter state
  /// @param newVertex Vertex to be added to fit
  /// @param linearizer Track linearizer
  /// @param vertexingOptions Vertexing options
  ///
  /// @return Result<void> object
  Result<void> addVtxToFit(
      State& state, Vertex<InputTrack_t>& newVertex,
      const Linearizer_t& linearizer,
      const VertexingOptions<InputTrack_t>& vertexingOptions) const;

 private:
  /// Configuration object
  const Config m_cfg;

  /// @brief Function to extract track parameters,
  /// InputTrack_t objects are BoundTrackParameters by default, function to be
  /// overwritten to return BoundTrackParameters for other InputTrack_t objects.
  std::function<BoundTrackParameters(InputTrack_t)> m_extractParameters;

  /// Logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private access to logging instance
  const Logger& logger() const { return *m_logger; }

  /// @brief Performs a simultaneous fit of all vertices in
  /// state.vertexCollection
  ///
  /// @param state Fitter state
  /// @param linearizer Track linearizer
  /// @param vertexingOptions Vertexing options
  ///
  /// @return Result<void> object
  Result<void> fitImpl(
      State& state, const Linearizer_t& linearizer,
      const VertexingOptions<InputTrack_t>& vertexingOptions) const;

  /// @brief Tests if vertex is already in list of vertices or not
  ///
  /// @param vtx Vertex to test
  /// @param verticesVec Vector of vertices to search
  ///
  /// @return True if vtx is already in verticesVec
  bool isAlreadyInList(
      Vertex<InputTrack_t>* vtx,
      const std::vector<Vertex<InputTrack_t>*>& verticesVec) const;

  /// @brief Checks whether the impact parameters of the associated
  /// tracks were calculated wrt the vertex position. Updates them
  /// if needed.
  ///
  /// @param state Vertex fitter state
  /// @param vtx Vertex object
  /// @param vertexingOptions Vertexing options
  Result<void> updateImpactParams3D(
      State& state, Vertex<InputTrack_t>* vtx,
      const VertexingOptions<InputTrack_t>& vertexingOptions) const;

  /// @brief Sets the vertexCompatibility for all TrackAtVertex objects
  /// at the current vertex
  ///
  /// @param state Fitter state
  /// @param currentVtx Current vertex
  /// @param vertexingOptions Vertexing options
  Result<void> setAllVertexCompatibilities(
      State& state, Vertex<InputTrack_t>* currentVtx,
      const VertexingOptions<input_track_t>& vertexingOptions) const;

  /// @brief Sets weights to the track according to Eq.(5.46) in Ref.(1)
  ///  and updates the vertices by calling the VertexUpdater
  ///
  /// @param state Fitter state
  /// @param linearizer The track linearizer
  /// @param vertexingOptions Vertexing options
  Result<void> setWeightsAndUpdate(
      State& state, const Linearizer_t& linearizer,
      const VertexingOptions<input_track_t>& vertexingOptions) const;

  /// @brief Collects the compatibility values of the track `trk`
  /// wrt to all of its associated vertices
  ///
  /// @param state Fitter state
  /// @param trk Track
  ///
  /// @return Vector of compatibility values
  std::vector<double> collectTrackToVertexCompatibilities(
      State& state, const InputTrack_t* trk) const;

  /// @brief Determines if any vertex position has shifted more than
  /// m_cfg.maxRelativeShift in the last iteration
  ///
  /// @param state Fitter state
  ///
  /// @return False if at least one shift was larger than maxRelativeShift
  bool checkSmallShift(State& state) const;

  /// @brief Updates tracks for current vertex with knowledge
  /// of current vertex position
  ///
  /// @param state Fitter state
  void doVertexSmoothing(State& state) const;
};

}  // namespace Acts

#include "Acts/Vertexing/AdaptiveMultiVertexFitter.ipp"
