// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/Result.hpp"
#include "Acts/Vertexing/ImpactPoint3dEstimator.hpp"
#include "Acts/Vertexing/KalmanVertexUpdator.hpp"
#include "Acts/Vertexing/LinearizedTrackFactory.hpp"
#include "Acts/Vertexing/MAVFInfo.hpp"
#include "Acts/Vertexing/SequentialVertexSmoother.hpp"
#include "Acts/Vertexing/VertexAnnealingTool.hpp"
#include "Acts/Vertexing/VertexFitterOptions.hpp"

#include <functional>

namespace Acts {

// TODO: add docs
///   Ref. (1): CERN-THESIS-2010-027, Author: Piacquadio, Giacinto:
///   `Identification of b-jets and investigation of the discovery potential
///   of a Higgs boson in the WH−−>lvbb¯ channel with the ATLAS experiment`
template <typename bfield_t, typename input_track_t, typename propagator_t>
class MultiAdaptiveVertexFitter {
 public:
  struct State {
    // Vertex collection to be fitted
    std::vector<Vertex<input_track_t>*> vertexCollection;

    // Annealing state
    VertexAnnealingTool::State annealingState;

    // map to store vertices information
    std::map<Vertex<input_track_t>*, MAVFVertexInfo<input_track_t>> vtxInfoMap;

    // map to store tracks information
    std::map<unsigned long, MAVFTrackAtVtxInfo<input_track_t>> trkInfoMap;
  };

  struct Config {

    /// @brief Config constructor
    ///
    /// @param bIn The magnetic field
    /// @param propagatorIn The propagator
    Config(const bfield_t& bIn, const propagator_t& propagatorIn)
        : propagator(propagatorIn),

          linFactory(
              typename LinearizedTrackFactory<bfield_t, propagator_t>::Config(
                  bIn)),
          ipEst(typename ImpactPoint3dEstimator<
                bfield_t, input_track_t, propagator_t>::Config(bIn,
                                                               propagatorIn)) {}
    /// Propagator
    propagator_t propagator;

    /// Linearized track factory
    LinearizedTrackFactory<bfield_t, propagator_t> linFactory;

    /// ImpactPoint3dEstimator
    ImpactPoint3dEstimator<input_track_t> ipEst;

    /// Vertex updator
    KalmanVertexUpdator<input_track_t> vertexUpdator;

    /// SequentialVertexSmoother
    SequentialVertexSmoother<input_track_t> vertexSmoother;

    /// Annealing tool
    VertexAnnealingTool annealingTool;

    /// Number of max iterations
    int maxIterations{50};

    /// Max distance to linearization point allowed
    /// without relinearization
    double maxDistToLinPoint{0.5};

    /// Minimum track weight needed for track to be considered
    double minWeight{0.001};

    /// Max relative shift of vertex during one iteration
    double maxRelativeShift{0.01};

    /// Do smoothing after multivertex fit
    bool doSmoothing{false};
  };

  /// @brief Constructor used if input_track_t type == BoundParameters
  ///
  /// @param cfg Configuration object
  /// @param logger The logging instance
  template <typename T = input_track_t,
            std::enable_if_t<std::is_same<T, BoundParameters>::value, int> = 0>
  MultiAdaptiveVertexFitter(Config& cfg,
                            std::unique_ptr<const Logger> logger =
                                getDefaultLogger("MultiAdaptiveVertexFitter",
                                                 Logging::INFO))
      : m_cfg(std::move(cfg)),
        m_extractParameters([](T params) { return params; }),
        m_logger(std::move(logger)) {}

  /// @brief Constructor for user-defined input_track_t type =! BoundParameters
  ///
  /// @param cfg Configuration object
  /// @param func Function extracting BoundParameters from input_track_t object
  MultiAdaptiveVertexFitter(Config& cfg,
                            std::function<BoundParameters(input_track_t)> func,
                            std::unique_ptr<const Logger> logger =
                                getDefaultLogger("MultiAdaptiveVertexFitter",
                                                 Logging::INFO))
      : m_cfg(std::move(cfg)),
        m_extractParameters(func),
        m_logger(std::move(logger)) {}

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
      State& state, Vertex<input_track_t>& newVertex,
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
  bool isAlreadyInList(
      Vertex<input_track_t>* vtx,
      const std::vector<Vertex<input_track_t>*>& verticesVec) const;

  /// @brief Prepares vertex object for the actual fit, i.e.
  /// all TrackAtVertex objects at current vertex will obtain
  /// `ip3dParams` from ImpactPoint3dEstimator::getParamsAtIP3d
  /// in order to later faster estimate compatibilities of track
  /// with different vertices
  ///
  /// @param vtx The vertex object
  /// @param vFitterOptions Vertex fitter options
  Result<void> prepareVtxForFit(
      State& state, Vertex<input_track_t>* vtx,
      const VertexFitterOptions<input_track_t>& vFitterOptions) const;

  /// @brief Sets vertexCompatibility for all TrackAtVertex objects
  /// at current vertex
  ///
  /// @param state The state object
  /// @param geoContext The geometry context
  /// @param mfContext The magnetic field context
  /// @param currentVtx Current vertex
  Result<void> setAllVtxCompatibilities(
      State& state, const GeometryContext& geoContext,
      const MagneticFieldContext& mfContext,
      Vertex<input_track_t>* currentVtx) const;

  /// @brief Sets weights to the track according to Eq.(5.46) in Ref.(1)
  ///  and updates the vertices by calling the VertexUpdator
  ///
  /// @param state The state object
  /// @param vFitterOptions Vertex fitter options
  Result<void> setWeightsAndUpdate(
      State& state,
      const VertexFitterOptions<input_track_t>& vFitterOptions) const;

  /// @brief Collects all compatibility values of the track `trk`
  /// at all vertices it is currently attached to and outputs
  /// these values in a vector
  ///
  /// @param state The state object
  /// @param trk The track
  ///
  /// @return Vector of compatibility values
  Result<std::vector<double>> collectTrkToVtxCompatibilities(
      State& state, const TrackAtVertex<input_track_t>& trk) const;

  /// @brief Determines if vertex position has shifted more than
  /// m_cfg.maxRelativeShift in last iteration
  ///
  /// @param state The state object
  ///
  /// @return False if shift was larger than maxRelativeShift
  bool checkSmallShift(State& state) const;
};

}  // namespace Acts

#include "Acts/Vertexing/MultiAdaptiveVertexFitter.ipp"
