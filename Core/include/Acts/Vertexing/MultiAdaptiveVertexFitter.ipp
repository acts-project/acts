// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/VertexingError.hpp"

template <typename bfield_t, typename input_track_t, typename propagator_t>
Acts::Result<void>
Acts::MultiAdaptiveVertexFitter<bfield_t, input_track_t, propagator_t>::fit(
    State& state,
    const VertexFitterOptions<input_track_t>& vFitterOptions) const {
  auto& geoContext = vFitterOptions.geoContext;
  auto& mfContext = vFitterOptions.magFieldContext;

  // Reset annealing tool
  m_cfg.annealingTool.reset(state.annealingState);

  // Indicates how much the vertex positions have shifted
  // in last fit iteration. Will be false if vertex position
  // shift was too big. Needed if equilibrium is reached in
  // annealing procedure but fitter has not fully converged
  // yet and needs some more iterations until vertex position
  // shifts between iterations are small (converged).
  bool isSmallShift = true;

  // Number of iterations counter
  unsigned int nIter = 0;

  // Start iterating
  while (nIter < m_cfg.maxIterations &&
         (!state.annealingState.equilibriumReached || !isSmallShift)) {
    // Initial loop over all vertices in state.vertexCollection
    for (auto currentVtx : state.vertexCollection) {
      MAVFVertexInfo<input_track_t>& currentVtxInfo =
          state.vtxInfoMap[currentVtx];
      currentVtxInfo.relinearize = false;

      // Store old position of vertex, i.e. seed position
      // in case of first iteration or position determined
      // in previous iteration afterwards
      currentVtxInfo.oldPosition = currentVtx->fullPosition();

      // Determine if relinearization is needed
      if ((currentVtxInfo.oldPosition - currentVtxInfo.linPoint).norm() >
          m_cfg.maxDistToLinPoint) {
        // Relinearization needed, distance too big
        currentVtxInfo.relinearize = true;
        // Prepare for fit with new vertex position
        prepareVtxForFit(state, currentVtx, vFitterOptions);
      }

      // Determine if constraint vertex exist
      if (state.vtxInfoMap[currentVtx]
              .constraintVertex.fullCovariance()
              .determinant() != 0) {
        currentVtx->setFullPosition(
            state.vtxInfoMap[currentVtx].constraintVertex.fullPosition());
        currentVtx->setFitQuality(
            state.vtxInfoMap[currentVtx].constraintVertex.fitQuality());
        currentVtx->setFullCovariance(
            state.vtxInfoMap[currentVtx].constraintVertex.fullCovariance());
      }

      else if (currentVtx->fullCovariance().determinant() == 0.) {
        return VertexingError::NoCovariance;
      }

      currentVtx->setFullCovariance(
          currentVtx->fullCovariance() * 1. /
          m_cfg.annealingTool.getWeight(state.annealingState, 1.));

      // Set vertexCompatibility for all TrackAtVertex objects
      // at current vertex
      setAllVtxCompatibilities(state, geoContext, mfContext, currentVtx);
    }  // End loop over vertex collection

    // Now after having estimated all compatibilities of all tracks at
    // all vertices, run again over all vertices to set track weights
    // and update the vertex
    setWeightsAndUpdate(state, vFitterOptions);

    if (!state.annealingState.equilibriumReached) {
      m_cfg.annealingTool.anneal(state.annealingState);
    }

    isSmallShift = checkSmallShift(state);

    ++nIter;
  }
  // Multivertex fit is finished

  // Check if smoothing is required
  if (m_cfg.doSmoothing) {
    for (auto vtx : state.vertexCollection) {
      // Smooth all tracks at vertex `vtx`
      m_cfg.vertexSmoother.smooth(geoContext, vtx);
    }
  }

  return {};
}

template <typename bfield_t, typename input_track_t, typename propagator_t>
Acts::Result<void>
Acts::MultiAdaptiveVertexFitter<bfield_t, input_track_t, propagator_t>::
    addVertexToFit(
        State& state, Vertex<input_track_t>& newVertex,
        const VertexFitterOptions<input_track_t>& vFitterOptions) const {
  if (newVertex.tracks().empty()) {
    return VertexingError::EmptyInput;
  }

  std::vector<Vertex<input_track_t>*> verticesToFit;

  // Prepares vtx and tracks for fast estimation method of their
  // compatibility with vertex
  auto res = prepareVtxForFit(state, &newVertex, vFitterOptions);
  if (!res.ok()) {
    return res.error();
  }

  // List of vertices added in last iteration
  std::vector<Vertex<input_track_t>*> lastIterAddedVertices = {&newVertex};
  // List of vertices added in current iteration
  std::vector<Vertex<input_track_t>*> currentIterAddedVertices;

  // Loop as long as new vertices are found that share tracks with
  // previously added vertices
  while (!lastIterAddedVertices.empty()) {
    for (auto& lastVtxIter : lastIterAddedVertices) {
      // Loop over all track at current lastVtxIter
      for (const TrackAtVertex<input_track_t>& trackIter :
           lastVtxIter->tracks()) {
        // Retrieve list of links to all vertices that currently use the current
        // track
        std::vector<Vertex<input_track_t>*>& linksToVertices =
            state.trkInfoMap[trackIter.id].linksToVertices;

        // Loop over all attached vertices and add those to vertex fit
        // which are not already in `verticesToFit`
        for (auto newVtxIter : linksToVertices) {
          if (!isAlreadyInList(newVtxIter, verticesToFit)) {
            // Add newVtxIter to verticesToFit
            verticesToFit.push_back(newVtxIter);

            // Add newVtxIter vertex to currentIterAddedVertices
            // if vertex != lastVtxIter
            if (newVtxIter != lastVtxIter) {
              currentIterAddedVertices.push_back(newVtxIter);
            }
          }
        }  // End for loop over linksToVertices
      }
    }  // End loop over lastIterAddedVertices

    lastIterAddedVertices = currentIterAddedVertices;
    currentIterAddedVertices.clear();
  }  // End while loop

  state.vertexCollection = verticesToFit;

  // Perform fit on all added vertices
  auto fitRes = fit(state, vFitterOptions);

  if (!fitRes.ok()) {
    return fitRes.error();
  }

  return {};
}

template <typename bfield_t, typename input_track_t, typename propagator_t>
bool Acts::MultiAdaptiveVertexFitter<bfield_t, input_track_t, propagator_t>::
    isAlreadyInList(
        Vertex<input_track_t>* vtx,
        const std::vector<Vertex<input_track_t>*>& verticesVec) const {
  return std::find_if(verticesVec.begin(), verticesVec.end(),
                      [vtx](Vertex<input_track_t>* v) { return v == vtx; }) !=
         verticesVec.end();
}

template <typename bfield_t, typename input_track_t, typename propagator_t>
Acts::Result<void> Acts::MultiAdaptiveVertexFitter<
    bfield_t, input_track_t,
    propagator_t>::prepareVtxForFit(State& state, Vertex<input_track_t>* vtx,
                                    const VertexFitterOptions<input_track_t>&
                                        vFitterOptions) const {
  const Vector3D& refPos = vtx->position();
  auto& geoContext = vFitterOptions.geoContext;
  auto& mfContext = vFitterOptions.magFieldContext;

  // Loop over all tracks at current vertex
  for (const auto& trkAtVtx : vtx->tracks()) {
    auto res = m_cfg.ipEst.getParamsAtIP3d(
        geoContext, mfContext, m_extractParameters(trkAtVtx.originalTrack),
        refPos);
    if (!res.ok()) {
      return res.error();
    }
    // Set ip3dParams for current trackAtVertex
    state.trkInfoMap[trkAtVtx.id].ip3dParams = std::move(res.value());
  }
  return {};
}

template <typename bfield_t, typename input_track_t, typename propagator_t>
Acts::Result<void>
Acts::MultiAdaptiveVertexFitter<bfield_t, input_track_t, propagator_t>::
    setAllVtxCompatibilities(State& state, const GeometryContext& geoContext,
                             const MagneticFieldContext& mfContext,
                             Vertex<input_track_t>* currentVtx) const {
  MAVFVertexInfo<input_track_t>& currentVtxInfo = state.vtxInfoMap[currentVtx];
  // Create empty list of new TrackAtVertex objects
  // to be filled below. Needed due to constness of
  // tracksAtVertex list at vertex
  std::vector<TrackAtVertex<input_track_t>> newTracks;
  newTracks.reserve(currentVtx->tracks().size());

  // Loop over tracks at current vertex and
  // estimate compatibility with vertex
  for (auto& trkAtVtx : currentVtx->tracks()) {
    // Recover from cases where linearization point != 0 but
    // more tracks were added later on
    if (!state.trkInfoMap[trkAtVtx.id].ip3dParams) {
      auto res = m_cfg.ipEst.getParamsAtIP3d(
          geoContext, mfContext, m_extractParameters(trkAtVtx.originalTrack),
          VectorHelpers::position(currentVtxInfo.linPoint));
      if (!res.ok()) {
        return res.error();
      }
      // Set ip3dParams for current trackAtVertex
      state.trkInfoMap[trkAtVtx.id].ip3dParams = std::move(res.value());
    }

    // Create copy of current trackAtVertex in order
    // to modify it below
    newTracks.push_back(trkAtVtx);
    TrackAtVertex<input_track_t>* newTrkPtr = &(newTracks.back());

    // Set compatibility with current vertex
    newTrkPtr->vertexCompatibility =

        m_cfg.ipEst.getVtxCompatibility(
            geoContext, state.trkInfoMap[trkAtVtx.id].ip3dParams.get(),
            VectorHelpers::position(currentVtxInfo.oldPosition));
  }
  // Set list of updated tracks to current vertex
  currentVtx->setTracksAtVertex(newTracks);
  return {};
}

template <typename bfield_t, typename input_track_t, typename propagator_t>
Acts::Result<void> Acts::MultiAdaptiveVertexFitter<
    bfield_t, input_track_t,
    propagator_t>::setWeightsAndUpdate(State& state,
                                       const VertexFitterOptions<input_track_t>&
                                           vFitterOptions) const {
  for (auto vtx : state.vertexCollection) {
    // Create empty list of new TrackAtVertex objects
    // to be filled below. Needed due to constness of
    // tracksAtVertex list at vertex
    std::vector<TrackAtVertex<input_track_t>> newTracks;
    newTracks.reserve(vtx->tracks().size());

    for (const auto& trkAtVtx : vtx->tracks()) {
      // Create copy of current trackAtVertex in order
      // to modify it below
      newTracks.push_back(trkAtVtx);
      TrackAtVertex<input_track_t>* newTrkPtr = &(newTracks.back());

      // Get all compatibilities of track to all vertices it is attached to
      auto collectRes = collectTrkToVtxCompatibilities(state, trkAtVtx);
      if (!collectRes.ok()) {
        return collectRes.error();
      }

      // Set trackWeight for current track
      newTrkPtr->trackWeight = m_cfg.annealingTool.getWeight(
          state.annealingState, trkAtVtx.vertexCompatibility, *collectRes);

      if (newTrkPtr->trackWeight > m_cfg.minWeight) {
        // Check if linearization state exists or need to be relinearized
        if (newTrkPtr->linearizedState.covarianceAtPCA ==
                BoundSymMatrix::Zero() ||
            state.vtxInfoMap[vtx].relinearize) {
          const auto& origParams =
              m_extractParameters(newTrkPtr->originalTrack);
          auto result = m_cfg.linFactory.linearizeTrack(
              vFitterOptions.geoContext, vFitterOptions.magFieldContext,
              &origParams, state.vtxInfoMap[vtx].oldPosition, m_cfg.propagator);
          if (!result.ok()) {
            return result.error();
          }
          newTrkPtr->linearizedState = *result;
          state.vtxInfoMap[vtx].linPoint = state.vtxInfoMap[vtx].oldPosition;
        }
        // Update the vertex with the new track
        m_cfg.vertexUpdator.addAndUpdate(vtx, (*newTrkPtr));
      } else {
        ACTS_VERBOSE("Track weight too low. Skip track.");
      }

    }  // End loop over tracks at vertex

    // Update tracks at current vertex
    vtx->setTracksAtVertex(newTracks);

    ACTS_VERBOSE("New vertex position: " << vtx->fullPosition());
  }  // End loop over vertex collection

  return {};
}

template <typename bfield_t, typename input_track_t, typename propagator_t>
Acts::Result<std::vector<double>>
Acts::MultiAdaptiveVertexFitter<bfield_t, input_track_t, propagator_t>::
    collectTrkToVtxCompatibilities(
        State& state, const TrackAtVertex<input_track_t>& trk) const {
  // All vertices that currently hold the track `trk`
  std::vector<Vertex<input_track_t>*> vertices =
      state.trkInfoMap[trk.id].linksToVertices;

  // Vector to store all compatibility values, it will have
  // exactly the size of `vertices`(one value for each vertex
  // the track is attached to)
  std::vector<double> trkToVtxCompatibilities;
  trkToVtxCompatibilities.reserve(vertices.size());

  for (Vertex<input_track_t>* vtxPtr : vertices) {
    // find current track in list of tracks at vertex
    const auto& trkIter = std::find_if(
        vtxPtr->tracks().begin(), vtxPtr->tracks().end(),
        [&trk, this](auto& trkAtVtx) {
          return this->m_extractParameters(trkAtVtx.originalTrack) ==
                 this->m_extractParameters(trk.originalTrack);
        });
    if (trkIter == vtxPtr->tracks().end()) {
      return VertexingError::ElementNotFound;
    }
    // store vertexCompatibility of track to current vertex
    trkToVtxCompatibilities.push_back(trkIter->vertexCompatibility);
  }
  return trkToVtxCompatibilities;
}

template <typename bfield_t, typename input_track_t, typename propagator_t>
bool Acts::MultiAdaptiveVertexFitter<bfield_t, input_track_t, propagator_t>::
    checkSmallShift(State& state) const {
  for (auto vtx : state.vertexCollection) {
    SpacePointVector diff =
        state.vtxInfoMap[vtx].oldPosition - vtx->fullPosition();
    const SpacePointSymMatrix& vtxWgt = vtx->fullCovariance().inverse();
    double relativeShift = diff.dot(vtxWgt * diff);
    if (relativeShift > m_cfg.maxRelativeShift) {
      return false;
    }
  }
  return true;
}
