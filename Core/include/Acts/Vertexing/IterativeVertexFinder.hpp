// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/FsmwMode1dFinder.hpp"
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexFitterConcept.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"
#include "Acts/Vertexing/ZScanVertexFinder.hpp"

namespace Acts {

/// @class IterativeVertexFinder
///
/// @brief Implements an iterative vertex finder
///
////////////////////////////////////////////////////////////
///
/// Brief description of the algorithm implemented:
/// Iterative vertex finder which iteratively finds and fits vertices:
/// 1. A list of seed tracks (`seedTracks`, which is the same as the
///   input track list to the finder at the very first iteration) is used
///   to retrieve a single vertex seed using the ZScanVertexFinder.
/// 2. All tracks compatible with the current vertex seed are kept and used
///   for fitting the single vertex.
/// 3.1 If the vertex is a 'good' vertex (i.e. meets requirements) and no
///   track reassignment after first fit is required, go to step 4. If vertex
///   is not a good vertex, remove all tracks in tracksToFit from seedTracks.
/// 3.2 If vertex meets requirements and track reassignment after first fit
///   is required, iterate over all previously found vertices ("old vertex")
///   and over all their tracksAtVertex. Compare compatibility of each track
///   with old vertex and current vertex. If track is more compatible with
///   current vertex, remove track from old vertex, put track back to
///   tracksToFit and refit current vertex with additional track.
/// 4. If good vertex, `removeUsedCompatibleTracks` method is called, which
///   removes all used tracks that are compatible with the fitted vertex
///   from `tracksToFit` and `seedTracks`. It also removes outliers tracks
///   from tracksAtVertex if not compatible.
/// 5. Add vertex to vertexCollection
/// 6. Repeat until no seedTracks are left or max. number of vertices found
///
////////////////////////////////////////////////////////////
///
/// @tparam vfitter_t Vertex fitter type
/// @tparam sfinder_t Seed finder type
template <typename vfitter_t, typename sfinder_t>
class IterativeVertexFinder {
  static_assert(VertexFitterConcept<vfitter_t>,
                "Vertex fitter does not fulfill vertex fitter concept.");
  using Propagator_t = typename vfitter_t::Propagator_t;
  using Linearizer_t = typename vfitter_t::Linearizer_t;

 public:
  using InputTrack_t = typename vfitter_t::InputTrack_t;
  using IPEstimator = ImpactPointEstimator<InputTrack_t, Propagator_t>;

  /// Configuration struct
  struct Config {
    /// @brief Config constructor
    ///
    /// @param fitter Vertex fitter
    /// @param lin Track linearizer
    /// @param sfinder The seed finder
    /// @param est ImpactPointEstimator
    Config(const vfitter_t& fitter, Linearizer_t lin, sfinder_t sfinder,
           const IPEstimator& est)
        : vertexFitter(fitter),
          linearizer(std::move(lin)),
          seedFinder(std::move(sfinder)),
          ipEst(est) {}

    /// Vertex fitter
    vfitter_t vertexFitter;

    /// Linearized track factory
    Linearizer_t linearizer;

    /// Vertex seed finder
    sfinder_t seedFinder;

    /// ImpactPointEstimator
    IPEstimator ipEst;

    /// Vertex finder configuration variables.
    /// Tracks that are within a distance of
    ///
    /// significanceCutSeeding * sqrt(sigma(d0)^2+sigma(z0)^2)
    ///
    /// are considered compatible with the vertex.
    double significanceCutSeeding = 10;
    double maximumChi2cutForSeeding = 36.;
    int maxVertices = 50;

    /// Assign a certain fraction of compatible tracks to a different (so-called
    /// split) vertex if boolean is set to true.
    bool createSplitVertices = false;
    /// Inverse of the fraction of tracks that will be assigned to the split
    /// vertex. E.g., if splitVerticesTrkInvFraction = 2, about 50% of
    /// compatible tracks will be assigned to the split vertex.
    int splitVerticesTrkInvFraction = 2;
    bool reassignTracksAfterFirstFit = false;
    bool doMaxTracksCut = false;
    int maxTracks = 5000;
    double cutOffTrackWeight = 0.01;
    /// If `reassignTracksAfterFirstFit` is set this threshold will be used to
    /// decide if a track should be checked for reassignment to other vertices
    double cutOffTrackWeightReassign = 1;
  };

  /// State struct
  struct State {
    State(const MagneticFieldProvider& field,
          const Acts::MagneticFieldContext& magContext)
        : ipState(field.makeCache(magContext)),
          linearizerState(field.makeCache(magContext)),
          fitterState(field.makeCache(magContext)) {}
    /// The IP estimator state
    typename IPEstimator::State ipState;
    /// The inearizer state
    typename Linearizer_t::State linearizerState;
    /// The fitter state
    typename vfitter_t::State fitterState;
  };

  /// @brief Constructor used if InputTrack_t type == BoundTrackParameters
  ///
  /// @param cfg Configuration object
  /// @param logger The logging instance
  template <
      typename T = InputTrack_t,
      std::enable_if_t<std::is_same<T, BoundTrackParameters>::value, int> = 0>
  IterativeVertexFinder(Config& cfg,
                        std::unique_ptr<const Logger> logger = getDefaultLogger(
                            "IterativeVertexFinder", Logging::INFO))
      : m_cfg(std::move(cfg)),
        m_extractParameters([](T params) { return params; }),
        m_logger(std::move(logger)) {}

  /// @brief Constructor for user-defined InputTrack_t type =!
  /// BoundTrackParameters
  ///
  /// @param cfg Configuration object
  /// @param func Function extracting BoundTrackParameters from InputTrack_t
  /// object
  /// @param logger The logging instance
  IterativeVertexFinder(Config& cfg,
                        std::function<BoundTrackParameters(InputTrack_t)> func,
                        std::unique_ptr<const Logger> logger = getDefaultLogger(
                            "IterativeVertexFinder", Logging::INFO))
      : m_cfg(std::move(cfg)),
        m_extractParameters(func),
        m_logger(std::move(logger)) {}

  /// @brief Finds vertices corresponding to input trackVector
  ///
  /// @param trackVector Input tracks
  /// @param vertexingOptions Vertexing options
  /// @param state State for fulfilling interfaces
  ///
  /// @return Collection of vertices found by finder
  Result<std::vector<Vertex<InputTrack_t>>> find(
      const std::vector<const InputTrack_t*>& trackVector,
      const VertexingOptions<InputTrack_t>& vertexingOptions,
      State& state) const;

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

  /// @brief Method that calls seed finder to retrieve a vertex seed
  ///
  /// @param seedTracks Seeding tracks
  /// @param vertexingOptions Vertexing options
  Result<Vertex<InputTrack_t>> getVertexSeed(
      const std::vector<const InputTrack_t*>& seedTracks,
      const VertexingOptions<InputTrack_t>& vertexingOptions) const;

  /// @brief Removes all tracks in tracksToRemove from seedTracks
  ///
  /// @param tracksToRemove Tracks to be removed from seedTracks
  /// @param seedTracks List to remove tracks from
  void removeTracks(const std::vector<const InputTrack_t*>& tracksToRemove,
                    std::vector<const InputTrack_t*>& seedTracks) const;

  /// @brief Function for calculating how compatible
  /// a given track is to a given vertex
  ///
  /// @param params Track parameters
  /// @param vertex The vertex
  /// @param perigeeSurface The perigee surface at vertex position
  /// @param vertexingOptions Vertexing options
  /// @param state The state object
  Result<double> getCompatibility(
      const BoundTrackParameters& params, const Vertex<InputTrack_t>& vertex,
      const Surface& perigeeSurface,
      const VertexingOptions<InputTrack_t>& vertexingOptions,
      State& state) const;

  /// @brief Function that removes used tracks compatible with
  /// current vertex (`vertex`) from `tracksToFit` and `seedTracks`
  /// as well as outliers from vertex.tracksAtVertex
  ///
  /// @param vertex Current vertex
  /// @param tracksToFit Tracks used to fit `vertex`
  /// @param seedTracks Tracks used for vertex seeding
  /// @param vertexingOptions Vertexing options
  /// @param state The state object
  Result<void> removeUsedCompatibleTracks(
      Vertex<InputTrack_t>& vertex,
      std::vector<const InputTrack_t*>& tracksToFit,
      std::vector<const InputTrack_t*>& seedTracks,
      const VertexingOptions<InputTrack_t>& vertexingOptions,
      State& state) const;

  /// @brief Function that fills vector with tracks compatible with seed vertex
  ///
  /// @param seedTracks List of all available tracks used for seeding
  /// @param seedVertex Seed vertex
  /// @param tracksToFitOut Tracks to fit
  /// @param tracksToFitSplitVertexOut Tracks to fit to split vertex
  /// @param vertexingOptions Vertexing options
  /// @param state The state object
  Result<void> fillTracksToFit(
      const std::vector<const InputTrack_t*>& seedTracks,
      const Vertex<InputTrack_t>& seedVertex,
      std::vector<const InputTrack_t*>& tracksToFitOut,
      std::vector<const InputTrack_t*>& tracksToFitSplitVertexOut,
      const VertexingOptions<InputTrack_t>& vertexingOptions,
      State& state) const;

  /// @brief Function that reassigns tracks from other vertices
  ///        to the current vertex if they are more compatible
  ///
  /// @param vertexCollection Collection of vertices
  /// @param currentVertex Current vertex to assign tracks to
  /// @param tracksToFit Tracks to fit vector
  /// @param seedTracks Seed tracks vector
  /// @param origTracks Vector of original track objects
  /// @param vertexingOptions Vertexing options
  /// @param state The state object
  ///
  /// @return Bool if currentVertex is still a good vertex
  Result<bool> reassignTracksToNewVertex(
      std::vector<Vertex<InputTrack_t>>& vertexCollection,
      Vertex<InputTrack_t>& currentVertex,
      std::vector<const InputTrack_t*>& tracksToFit,
      std::vector<const InputTrack_t*>& seedTracks,
      const std::vector<const InputTrack_t*>& origTracks,
      const VertexingOptions<InputTrack_t>& vertexingOptions,
      State& state) const;

  /// @brief Counts all tracks that are significant for a vertex
  ///
  /// @param vtx The vertex
  ///
  /// @return Number of significant tracks
  int countSignificantTracks(const Vertex<InputTrack_t>& vtx) const;
};

}  // namespace Acts

#include "IterativeVertexFinder.ipp"
