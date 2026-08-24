// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Utilities/AnnealingUtility.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/AMVFInfo.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/TrackAtVertex.hpp"
#include "Acts/Vertexing/TrackLinearizer.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexFitProblem.hpp"
#include "Acts/Vertexing/VertexingError.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"

#include <algorithm>

namespace Acts {

/// @class AdaptiveMultiVertexFitter
/// @brief Implements an adaptive multi-vertex fitter as described
///   in detail in Section 5.3.5 in:
///   Ref. (1): CERN-THESIS-2010-027, Author: Piacquadio, Giacinto:
///   `Identification of b-jets and investigation of the discovery potential
///   of a Higgs boson in the WH−−>lvbb¯ channel with the ATLAS experiment`
///
class AdaptiveMultiVertexFitter {
 public:
  /// @brief The fitter-private cache
  ///
  /// Holds everything that is scratch data of a running fit: the annealing
  /// state, the impact point estimator state, the magnetic field cache and the
  /// per-vertex linearization bookkeeping. The data describing *what* is being
  /// fitted lives in @c VertexFitProblem instead.
  struct Cache {
    /// Constructor for the multi-vertex fitter cache
    /// @param field Magnetic field provider for track extrapolation
    /// @param magContext Magnetic field context for field evaluations
    Cache(const MagneticFieldProvider& field,
          const Acts::MagneticFieldContext& magContext)
        : ipState{field.makeCache(magContext)},
          fieldCache(field.makeCache(magContext)) {}

    /// Annealing state for thermodynamic track weighting
    AnnealingUtility::State annealingState;

    /// Impact point estimator state for track parameter calculations
    ImpactPointEstimator::State ipState;

    /// Magnetic field cache for field evaluations during fitting
    MagneticFieldProvider::Cache fieldCache;

    /// Per-vertex scratch data, lazily created by @c scratchFor
    std::map<const Vertex*, VertexScratch> vertexScratch;

    /// Drop the scratch data associated with @p vtx
    ///
    /// Must be called whenever a vertex's candidate description is re-seeded
    /// after it has already been fitted, since the cached linearization point
    /// and impact parameters then refer to a stale seed position.
    /// @param vtx Vertex whose scratch data is to be discarded
    void invalidate(const Vertex& vtx) { vertexScratch.erase(&vtx); }

   private:
    /// Construct a cache adopting an already existing impact point estimator
    /// state and magnetic field cache. Used by the deprecated @c State based
    /// entry points, which own those two caches themselves.
    /// @param ipStateIn Impact point estimator state to adopt
    /// @param fieldCacheIn Magnetic field cache to adopt
    Cache(ImpactPointEstimator::State ipStateIn,
          MagneticFieldProvider::Cache fieldCacheIn)
        : ipState(std::move(ipStateIn)), fieldCache(std::move(fieldCacheIn)) {}

    friend class AdaptiveMultiVertexFitter;
  };

  /// @brief Configuration options for the adaptive multi-vertex fit.
  struct Config {
    /// @brief Config constructor
    ///
    /// @param est ImpactPointEstimator
    explicit Config(ImpactPointEstimator est) : ipEst(std::move(est)) {}

    /// ImpactPointEstimator
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

    /// Number of max iterations
    unsigned int maxIterations{30};

    /// Max distance to linearization point allowed
    /// without relinearization
    double maxDistToLinPoint{0.5};

    /// Minimum track weight needed for track to be considered
    double minWeight{0.0001};

    /// Max relative shift of vertex during one iteration
    double maxRelativeShift{0.01};

    /// Do smoothing after multivertex fit
    bool doSmoothing{false};

    /// Use time information when calculating the vertex compatibility
    bool useTime{false};

    /// Function to extract parameters from InputTrack
    InputTrack::Extractor extractParameters;

    /// TrackLinearizer
    TrackLinearizer trackLinearizer;
  };

  /// @brief Constructor for user-defined InputTrack_t type !=
  /// BoundTrackParameters
  ///
  /// @param cfg Configuration object
  /// object
  /// @param logger The logging instance
  explicit AdaptiveMultiVertexFitter(
      Config cfg, std::unique_ptr<const Logger> logger = getDefaultLogger(
                      "AdaptiveMultiVertexFitter", Logging::INFO));

  /// @brief Adds a new vertex to an existing multi-vertex fit.
  /// 1. The 3D impact parameters are calculated for all tracks associated
  /// with newVertex.
  /// 2. A list of all vertices that are connected with newVertex via shared
  /// tracks is created. It also considers indirect connections (e.g., vertex A
  /// which shares a track with vertex B which, in turn, shares a track with
  /// newVertex).
  /// 3. The multivertex fit is performed for all vertices on said list.
  ///
  /// @param problem The multi-vertex fit problem
  /// @param newVertices Vertex to be added to fit
  /// @param vertexingOptions Vertexing options
  /// @param cache Fitter cache
  ///
  /// @return Result<void> object
  Result<void> addVtxToFit(VertexFitProblem& problem,
                           const std::vector<Vertex*>& newVertices,
                           const VertexingOptions& vertexingOptions,
                           Cache& cache) const;

  /// @brief Performs a simultaneous fit of all vertices in
  /// problem.vertices
  ///
  /// @param problem The multi-vertex fit problem
  /// @param vertexingOptions Vertexing options
  /// @param cache Fitter cache
  ///
  /// @return Result<void> object
  Result<void> fit(VertexFitProblem& problem,
                   const VertexingOptions& vertexingOptions,
                   Cache& cache) const;

  // Backwards compatibility: the pre-split fitter state and the entry points
  // taking it. Both are adapters on top of the interface above, which they
  // reproduce exactly: the state is decomposed into a problem and a cache on
  // the way in and recomposed from them on the way out.
  //
  // Only the entry points carry the deprecation attribute. The types are
  // documented as deprecated but not attributed, so that naming them in the
  // adapter's own declarations does not warn -- a caller reaching the legacy
  // interface has to go through fit() or addVtxToFit() anyway.

  /// @brief The fitter state
  ///
  /// @deprecated Split into @c VertexFitProblem, which describes the vertices
  /// to be fitted and is owned by the caller, and @c Cache, which holds the
  /// fitter's scratch data. Construct those two separately and use the
  /// corresponding @c fit and @c addVtxToFit overloads instead.
  struct State {
    /// Constructor for multi-vertex fitter state
    /// @param field Magnetic field provider for track extrapolation
    /// @param magContext Magnetic field context for field evaluations
    State(const MagneticFieldProvider& field,
          const Acts::MagneticFieldContext& magContext)
        : ipState{field.makeCache(magContext)},
          fieldCache(field.makeCache(magContext)) {}

    /// Vertex collection to be fitted
    std::vector<Vertex*> vertexCollection;

    /// Annealing state for thermodynamic track weighting
    AnnealingUtility::State annealingState;

    /// Impact point estimator state for track parameter calculations
    ImpactPointEstimator::State ipState;

    /// Magnetic field cache for field evaluations during fitting
    MagneticFieldProvider::Cache fieldCache;

    /// Map storing vertex information for each vertex in the fit
    std::map<Vertex*, VertexInfo> vtxInfoMap;

    /// Multimap connecting tracks to their associated vertices
    std::multimap<InputTrack, Vertex*> trackToVerticesMultiMap;

    /// Map storing track-at-vertex information for each track-vertex pair
    std::map<std::pair<InputTrack, Vertex*>, TrackAtVertex> tracksAtVerticesMap;

    /// Adds a vertex to trackToVerticesMultiMap
    /// @param vtx Vertex to add to the multimap with its associated tracks
    void addVertexToMultiMap(Vertex& vtx) {
      for (auto trk : vtxInfoMap[&vtx].trackLinks) {
        trackToVerticesMultiMap.emplace(trk, &vtx);
      }
    }

    /// Removes a vertex from trackToVerticesMultiMap
    /// @param vtx Vertex to remove from the multimap along with its track associations
    void removeVertexFromMultiMap(const Vertex& vtx) {
      for (auto iter = trackToVerticesMultiMap.begin();
           iter != trackToVerticesMultiMap.end();) {
        if (iter->second == &vtx) {
          iter = trackToVerticesMultiMap.erase(iter);
        } else {
          ++iter;
        }
      }
    }

    /// Remove a vertex from the vertex collection
    /// @param vtxToRemove Vertex to remove from the collection
    /// @param logger Logger for diagnostic messages
    /// @return Result indicating success or failure of the removal operation
    Result<void> removeVertexFromCollection(Vertex& vtxToRemove,
                                            const Logger& logger) {
      auto it = std::ranges::find(vertexCollection, &vtxToRemove);
      // Check if the value was found before erasing
      if (it == vertexCollection.end()) {
        ACTS_ERROR("vtxToRemove is not part of vertexCollection.");
        return VertexingError::ElementNotFound;
      }
      // Erase the element if found
      vertexCollection.erase(it);
      return {};
    }
  };

  /// @brief Adds a new vertex to an existing multi-vertex fit
  ///
  /// @deprecated Use the overload taking a @c VertexFitProblem and a @c Cache
  ///
  /// @param state Fitter state
  /// @param newVertices Vertex to be added to fit
  /// @param vertexingOptions Vertexing options
  ///
  /// @return Result<void> object
  [[deprecated(
      "Use addVtxToFit(VertexFitProblem&, const std::vector<Vertex*>&, const "
      "VertexingOptions&, Cache&)")]]
  Result<void> addVtxToFit(State& state,
                           const std::vector<Vertex*>& newVertices,
                           const VertexingOptions& vertexingOptions) const;

  /// @brief Performs a simultaneous fit of all vertices in
  /// state.vertexCollection
  ///
  /// @deprecated Use the overload taking a @c VertexFitProblem and a @c Cache
  ///
  /// @param state Fitter state
  /// @param vertexingOptions Vertexing options
  ///
  /// @return Result<void> object
  [[deprecated("Use fit(VertexFitProblem&, const VertexingOptions&, Cache&)")]]
  Result<void> fit(State& state,
                   const VertexingOptions& vertexingOptions) const;

 private:
  /// Configuration object
  const Config m_cfg;

  /// Logging instance
  std::unique_ptr<const Logger> m_logger;

  /// Private access to logging instance
  const Logger& logger() const { return *m_logger; }

  /// @brief Tests if vertex is already in list of vertices or not
  ///
  /// @param vtx Vertex to test
  /// @param verticesVec Vector of vertices to search
  ///
  /// @return True if vtx is already in verticesVec
  bool isAlreadyInList(Vertex* vtx,
                       const std::vector<Vertex*>& verticesVec) const;

  /// @brief Returns the scratch data for @p vtx, creating it on first access
  ///
  /// A newly created entry is seeded with the vertex's seed position, so that
  /// the first iteration measures the distance to the seed rather than to the
  /// origin.
  ///
  /// @param problem The multi-vertex fit problem
  /// @param vtx Vertex object
  /// @param cache Fitter cache
  ///
  /// @return Reference to the scratch data of @p vtx
  static VertexScratch& scratchFor(const VertexFitProblem& problem, Vertex* vtx,
                                   Cache& cache);

  /// @brief 1) Calls ImpactPointEstimator::estimate3DImpactParameters
  /// for all tracks that are associated with vtx (i.e., all elements
  /// of the trackLinks vector in the candidate of vtx).
  /// 2) Saves the 3D impact parameters in the scratch data of vtx.
  ///
  /// @param problem The multi-vertex fit problem
  /// @param vtx Vertex object
  /// @param vertexingOptions Vertexing options
  /// @param cache Fitter cache
  Result<void> prepareVertexForFit(const VertexFitProblem& problem, Vertex* vtx,
                                   const VertexingOptions& vertexingOptions,
                                   Cache& cache) const;

  /// @brief Sets the vertexCompatibility for all TrackAtVertex objects
  /// at the current vertex
  ///
  /// @param problem The multi-vertex fit problem
  /// @param currentVtx Current vertex
  /// @param vertexingOptions Vertexing options
  /// @param cache Fitter cache
  Result<void> setAllVertexCompatibilities(
      VertexFitProblem& problem, Vertex* currentVtx,
      const VertexingOptions& vertexingOptions, Cache& cache) const;

  /// @brief Sets weights to the track according to Eq.(5.46) in Ref.(1)
  ///  and updates the vertices by calling the VertexUpdater
  ///
  /// @param problem The multi-vertex fit problem
  /// @param vertexingOptions Vertexing options
  /// @param cache Fitter cache
  Result<void> setWeightsAndUpdate(VertexFitProblem& problem,
                                   const VertexingOptions& vertexingOptions,
                                   Cache& cache) const;

  /// @brief Collects the compatibility values of the track `trk`
  /// wrt to all of its associated vertices
  ///
  /// @param problem The multi-vertex fit problem
  /// @param trk Track
  ///
  /// @return Vector of compatibility values
  std::vector<double> collectTrackToVertexCompatibilities(
      const VertexFitProblem& problem, const InputTrack& trk) const;

  /// @brief Determines if any vertex position has shifted more than
  /// m_cfg.maxRelativeShift in the last iteration
  ///
  /// @param problem The multi-vertex fit problem
  /// @param cache Fitter cache
  ///
  /// @return False if at least one shift was larger than maxRelativeShift
  bool checkSmallShift(const VertexFitProblem& problem, Cache& cache) const;

  /// @brief Updates tracks for current vertex with knowledge
  /// of current vertex position
  ///
  /// @param problem The multi-vertex fit problem
  void doVertexSmoothing(VertexFitProblem& problem) const;

  /// @brief Logs vertices in problem.vertices and associated tracks
  ///
  /// @param problem The multi-vertex fit problem
  /// @param geoContext Geometry context
  /// @param cache Fitter cache
  void logDebugData(const VertexFitProblem& problem,
                    const GeometryContext& geoContext, Cache& cache) const;

  /// @brief Moves the contents of a deprecated state into a problem and a cache
  ///
  /// @param state Fitter state to take the contents from
  /// @param[out] problem The multi-vertex fit problem
  /// @param[out] cache Fitter cache
  static void decomposeState(State& state, VertexFitProblem& problem,
                             Cache& cache);

  /// @brief Moves a problem and a cache back into a deprecated state
  ///
  /// @param[out] state Fitter state to write the contents to
  /// @param problem The multi-vertex fit problem
  /// @param cache Fitter cache
  static void recomposeState(State& state, VertexFitProblem& problem,
                             Cache& cache);
};

}  // namespace Acts
