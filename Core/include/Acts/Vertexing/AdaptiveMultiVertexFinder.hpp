// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Vertexing/AMVFInfo.hpp"
#include "Acts/Vertexing/TrackToVertexIPEstimator.hpp"
#include "Acts/Vertexing/VertexFinderOptions.hpp"
#include "Acts/Vertexing/VertexFitterOptions.hpp"

namespace Acts {
/// @class AdaptiveMultiVertexFinder
///
/// @brief Implements an iterative vertex finder
///
////////////////////////////////////////////////////////////
///
/// Brief description of the algorithm implemented:
/// TODO
///
////////////////////////////////////////////////////////////
///
/// @tparam vfitter_t Vertex fitter type
/// @tparam sfinder_t Seed finder type
template <typename vfitter_t, typename sfinder_t>
class AdaptiveMultiVertexFinder {
  using Propagator_t = typename vfitter_t::Propagator_t;
  using InputTrack_t = typename vfitter_t::InputTrack_t;
  using Linearizer_t = typename vfitter_t::Linearizer_t;
  using FitterState_t = typename vfitter_t::State;

 public:
  /// @struct Config Configuration struct
  struct Config {
    /// @brief Config constructor
    ///
    /// @param fitter The vertex fitter
    /// @param sfinder The seed finder
    /// @param ipEst TrackToVertexIPEstimator
    /// @param lin Track linearizer
    Config(vfitter_t fitter, sfinder_t sfinder,
           TrackToVertexIPEstimator<InputTrack_t, Propagator_t> ipEst,
           Linearizer_t lin)
        : vertexFitter(std::move(fitter)),
          seedFinder(std::move(sfinder)),
          ipEstimator(std::move(ipEst)),
          linearizer(std::move(lin)) {}

    // Vertex fitter
    vfitter_t vertexFitter;

    // Vertex seed finder
    sfinder_t seedFinder;

    // TrackToVertexIPEstimator
    TrackToVertexIPEstimator<InputTrack_t, Propagator_t> ipEstimator;

    // Track linearizer
    Linearizer_t linearizer;

    // Use a beam spot constraint, vertexConstraint in VertexFinderOptions
    // has to be set in this case
    bool useBeamSpotConstraint = true;

    // Max z interval used for adding tracks to fit:
    // When adding a new vertex to the multi vertex fit,
    // only the tracks whose z at PCA is closer
    // to the seeded vertex than tracksMaxZinterval
    // are added to this new vertex.
    //
    // Note: If you cut too hard, you cut out
    // the good cases where the seed finder is not
    // reliable, but the fit would be still able to converge
    // towards the right vertex. If you cut too soft, you
    // consider a lot of tracks which just slow down the fit.
    double tracksMaxZinterval = 3. * Acts::UnitConstants::mm;

    // Maximum allowed significance of track position to vertex seed
    // to consider track as compatible track for vertex fit
    double tracksMaxSignificance = 5.;

    // Max chi2 value for which tracks are considered compatible with
    // the fitted vertex. These tracks are removed from the seedTracks
    // after the fit has been performed.
    double maxVertexChi2 = 18.42;

    // Perform a 'real' multi-vertex fit as intended by the algorithm.
    // If switched to true, always all (!) tracks are considered to be
    // added to the new vertex candidate after seeding. If switched to
    // false, only the seedTracks, i.e. all tracks that are considered
    // as outliers of previously fitted vertices, are used.
    bool doRealMultiVertex = true;

    // Decides if you want to use the ```vertexCompatibility``` of the
    //  track (set to true) or the ```chi2Track``` (set to false) as an
    // estimate for a track being an outlier or not.
    // In case the track refitting is switched on in the AMVFitter, you
    // may want to use the refitted ```chi2Track```.
    bool useFastCompatibility = true;

    // Maximum significance on the distance between two vertices
    // to allow merging of two vertices.
    double maxMergeVertexSignificance = 3.;

    // Minimum weight a track has to have to be considered a compatible
    // track with a vertex candidate.
    //
    // Note: This value has to be the same as the one in the AMVFitter.
    double minWeight = 0.0001;

    // Maximal number of iterations in the finding procedure
    int maxIterations = 100;

    // Include also single track vertices
    bool addSingleTrackVertices = false;

    // Use 3d information fo evaluating the vertex distance significance
    // for vertex merging/splitting
    bool do3dSplitting = false;

    // Maximum vertex contamination value
    double maximumVertexContamination = 0.5;

    // Use seed vertex as a constraint for the fit
    bool useSeedConstraint = true;

    // Diagonal constraint covariance entries in case
    // no beamspot constraint is provided
    double looseConstrValue = 1e+8;

    // Default fitQuality for constraint vertex in case no beamspot
    // constraint is provided
    std::pair<double, double> defaultConstrFitQuality{0., -3.};

    // Do an adaptive multi vertex fit after
    // a bad vertex was removed.
    // If false, the old fitter state is just copied,
    // this should give the same results with better
    // performance. To be further investigated.
    bool refitAfterBadVertex = true;

    // Use the full available vertex covariance information after
    // seeding for the IP estimation. In original implementation
    // this is not (!) done, however, this is probably not correct.
    // So definitely consider setting this to true.
    bool useVertexCovForIPEstimation = false;

  };  // Config struct

  /// @brief Constructor used if InputTrack_t type == BoundParameters
  ///
  /// @param cfg Configuration object
  /// @param logger The logging instance
  template <typename T = InputTrack_t,
            std::enable_if_t<std::is_same<T, BoundParameters>::value, int> = 0>
  AdaptiveMultiVertexFinder(Config& cfg,
                            std::unique_ptr<const Logger> logger =
                                getDefaultLogger("AdaptiveMultiVertexFinder",
                                                 Logging::INFO))
      : m_cfg(std::move(cfg)),
        m_extractParameters([](T params) { return params; }),
        m_logger(std::move(logger)) {}

  /// @brief Constructor for user-defined InputTrack_t type != BoundParameters
  ///
  /// @param cfg Configuration object
  /// @param func Function extracting BoundParameters from InputTrack_t object
  /// @param logger The logging instance
  AdaptiveMultiVertexFinder(Config& cfg,
                            std::function<BoundParameters(InputTrack_t)> func,
                            std::unique_ptr<const Logger> logger =
                                getDefaultLogger("AdaptiveMultiVertexFinder",
                                                 Logging::INFO))
      : m_cfg(std::move(cfg)),
        m_extractParameters(func),
        m_logger(std::move(logger)) {}

  /// @brief Function that performs the adaptive
  /// multi-vertex finding
  ///
  /// @param allTracks Input track collection
  /// @param vFinderOptions Vertex finder options
  ///
  /// @return Vector of all reconstructed vertices
  Result<std::vector<Vertex<InputTrack_t>>> find(
      const std::vector<const InputTrack_t*>& allTracks,
      const VertexFinderOptions<InputTrack_t>& vFinderOptions) const;

 private:
  /// Configuration object
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

  /// @brief Calls the seed finder and sets constraints on the found seed
  /// vertex if desired
  ///
  /// @param trackVector All tracks to be used for seeding
  /// @param currentConstraint Vertex constraint
  /// @param vFinderOptions Vertex finder options
  ///
  /// @return The seed vertex
  Result<Vertex<InputTrack_t>> doSeeding(
      const std::vector<const InputTrack_t*>& trackVector,
      Vertex<InputTrack_t>& currentConstraint,
      const VertexFinderOptions<InputTrack_t>& vFinderOptions) const;

  /// @brief Estimates delta Z between a track and a vertex position
  ///
  /// @param track The track
  /// @param vtxPos The vertex position
  ///
  /// @return The delta Z estimate
  double estimateDeltaZ(const BoundParameters& track,
                        const Vector3D& vtxPos) const;

  /// @brief Calculates the IP significance of a track to a given vertex
  ///
  /// @param track The track
  /// @param vtx The vertex
  ///
  /// @return The IP significance
  Result<double> getIPSignificance(const InputTrack_t* track,
                                   const Vertex<InputTrack_t>& vtx) const;

  /// @brief Adds compatible track to vertex candidate
  ///
  /// @param tracks The tracks
  /// @param vtx The vertex candidate
  /// @param[out] fitterState The vertex fitter state
  Result<void> addCompatibleTracksToVertex(
      const std::vector<const InputTrack_t*>& tracks, Vertex<InputTrack_t>& vtx,
      FitterState_t& fitterState) const;

  /// @brief Method that tries to recover from cases where no tracks
  /// were added to the vertex candidate after seeding
  ///
  /// @param allTracks The tracks to be considered (either origTrack or
  /// seedTracks)
  /// @param seedTracks The seed tracks
  /// @param[out] vtx The vertex candidate
  /// @param currentConstraint Vertex constraint
  /// @param[out] fitterState The vertex fitter state
  ///
  /// return True if recovery was successful, false otherwise
  Result<bool> canRecoverFromNoCompatibleTracks(
      const std::vector<const InputTrack_t*>& allTracks,
      const std::vector<const InputTrack_t*>& seedTracks,
      Vertex<InputTrack_t>& vtx, const Vertex<InputTrack_t>& currentConstraint,
      FitterState_t& fitterState) const;

  /// @brief Method that tries to prepare the vertex for the fit
  ///
  /// @param allTracks The tracks to be considered (either origTrack or
  /// seedTracks)
  /// @param seedTracks The seed tracks
  /// @param[out] vtx The vertex candidate
  /// @param currentConstraint Vertex constraint
  /// @param[out] fitterState The vertex fitter state
  ///
  /// @return True if preparation was successful, false otherwise
  Result<bool> canPrepareVertexForFit(
      const std::vector<const InputTrack_t*>& allTracks,
      const std::vector<const InputTrack_t*>& seedTracks,
      Vertex<InputTrack_t>& vtx, const Vertex<InputTrack_t>& currentConstraint,
      FitterState_t& fitterState) const;

  /// @brief Method that checks if vertex is a good vertex and if
  /// compatible tracks are available
  ///
  /// @param vtx The vertex candidate
  /// @param seedTracks The seed tracks
  /// @param fitterState The vertex fitter state
  ///
  /// @return pair(nCompatibleTracks, isGoodVertex)
  std::pair<int, bool> checkVertexAndCompatibleTracks(
      Vertex<InputTrack_t>& vtx,
      const std::vector<const InputTrack_t*>& seedTracks,
      FitterState_t& fitterState) const;

  /// @brief Method that removes all tracks that are compatible with
  /// current vertex from seedTracks
  ///
  /// @param vtx The vertex candidate
  /// @param[out] seedTracks The seed tracks
  /// @param fitterState The vertex fitter state
  void removeCompatibleTracksFromSeedTracks(
      Vertex<InputTrack_t>& vtx, std::vector<const InputTrack_t*>& seedTracks,
      FitterState_t& fitterState) const;

  /// @brief Method that tries to remove a non-compatible track
  /// from seed tracks after removing a compatible track failed.
  ///
  /// @param vtx The vertex candidate
  /// @param[out] seedTracks The seed tracks
  /// @param fitterState The vertex fitter state
  ///
  /// @return Non-compatible track was removed
  bool canRemoveNonCompatibleTrackFromSeedTracks(
      Vertex<InputTrack_t>& vtx, std::vector<const InputTrack_t*>& seedTracks,
      FitterState_t& fitterState) const;

  /// @brief Method that evaluates if the new vertex candidate should
  /// be kept, i.e. saved, or not
  ///
  /// @param vtx The vertex candidate
  /// @param allVertices All so far found vertices
  /// @param fitterState The vertex fitter state
  ///
  /// @return Keep new vertex
  bool keepNewVertex(Vertex<InputTrack_t>& vtx,
                     const std::vector<Vertex<InputTrack_t>*>& allVertices,
                     FitterState_t& fitterState) const;

  /// @brief Method that evaluates if the new vertex candidate is
  /// merged with one of the previously found vertices
  ///
  /// @param vtx The vertex candidate
  /// @param allVertices All so far found vertices
  ///
  /// @return Vertex is merged
  bool isMergedVertex(
      const Vertex<InputTrack_t>& vtx,
      const std::vector<Vertex<InputTrack_t>*>& allVertices) const;

  /// @brief Method that deletes last vertex from list of all vertices
  /// and either refits all vertices afterwards (if refitAfterBadVertex)
  /// of reverts to the old state of the vertex fitter before the bad
  /// vertex was added to the fit (if not refitAfterBadVertex).
  ///
  /// @param vtx The last added vertex which will be removed
  /// @param allVertices Vector containing the unique_ptr to vertices
  /// @param allVerticesPtr Vector containing the actual addresses
  /// @param fitterState The current vertex fitter state
  /// @param oldFitterState The old vertex fitter state
  /// @param vFitterOptions The vertex fitter options
  Result<void> deleteLastVertex(
      Vertex<InputTrack_t>& vtx,
      std::vector<std::unique_ptr<Vertex<InputTrack_t>>>& allVertices,
      std::vector<Vertex<InputTrack_t>*>& allVerticesPtr,
      FitterState_t& fitterState, FitterState_t& oldFitterState,
      const VertexFitterOptions<InputTrack_t>& vFitterOptions) const;

  /// @brief Prepares the output vector of vertices
  ///
  /// @param allVerticesPtr Vector of pointers to vertices
  /// @param fitterState The vertex fitter state
  ///
  /// @return The output vertex collection
  std::vector<Vertex<InputTrack_t>> getVertexOutputList(
      const std::vector<Vertex<InputTrack_t>*>& allVerticesPtr,
      FitterState_t& fitterState) const;
};

}  // namespace Acts

#include "Acts/Vertexing/AdaptiveMultiVertexFinder.ipp"
