// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/KalmanVertexTrackUpdater.hpp"
#include "Acts/Vertexing/KalmanVertexUpdater.hpp"
#include "Acts/Vertexing/VertexingError.hpp"

template <typename input_track_t, typename linearizer_t>
Acts::Result<void>
Acts::AdaptiveMultiVertexFitter<input_track_t, linearizer_t>::fit(
    State& state, const std::vector<Vertex<input_track_t>*>& verticesToFit,
    const linearizer_t& linearizer,
    const VertexFitterOptions<input_track_t>& vFitterOptions) const {
  // Set all vertices to fit in the current state
  state.vertexCollection = verticesToFit;

  // Perform fit on all vertices simultaneously
  auto fitRes = fitImpl(state, linearizer, vFitterOptions);

  if (!fitRes.ok()) {
    return fitRes.error();
  }

  return {};
}

template <typename input_track_t, typename linearizer_t>
Acts::Result<void>
Acts::AdaptiveMultiVertexFitter<input_track_t, linearizer_t>::fitImpl(
    State& state, const linearizer_t& linearizer,
    const VertexFitterOptions<input_track_t>& vFitterOptions) const {
  auto& geoContext = vFitterOptions.geoContext;

  // Reset annealing tool
  state.annealingState = AnnealingUtility::State();

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
      VertexInfo<input_track_t>& currentVtxInfo = state.vtxInfoMap[currentVtx];
      currentVtxInfo.relinearize = false;
      // Store old position of vertex, i.e. seed position
      // in case of first iteration or position determined
      // in previous iteration afterwards
      currentVtxInfo.oldPosition = currentVtx->fullPosition();

      SpacePointVector dist =
          currentVtxInfo.oldPosition - currentVtxInfo.linPoint;
      double perpDist = std::sqrt(dist[0] * dist[0] + dist[1] * dist[1]);
      // Determine if relinearization is needed
      if (perpDist > m_cfg.maxDistToLinPoint) {
        // Relinearization needed, distance too big
        currentVtxInfo.relinearize = true;
        // Prepare for fit with new vertex position
        prepareVertexForFit(state, currentVtx, vFitterOptions);
      }
      // Determine if constraint vertex exist
      if (state.vtxInfoMap[currentVtx].constraintVertex.fullCovariance() !=
          SpacePointSymMatrix::Zero()) {
        currentVtx->setFullPosition(
            state.vtxInfoMap[currentVtx].constraintVertex.fullPosition());
        currentVtx->setFitQuality(
            state.vtxInfoMap[currentVtx].constraintVertex.fitQuality());
        currentVtx->setFullCovariance(
            state.vtxInfoMap[currentVtx].constraintVertex.fullCovariance());
      }

      else if (currentVtx->fullCovariance() == SpacePointSymMatrix::Zero()) {
        return VertexingError::NoCovariance;
      }
      double weight =
          1. / m_cfg.annealingTool.getWeight(state.annealingState, 1.);

      currentVtx->setFullCovariance(currentVtx->fullCovariance() * weight);

      // Set vertexCompatibility for all TrackAtVertex objects
      // at current vertex
      setAllVertexCompatibilities(state, geoContext, currentVtx);
    }  // End loop over vertex collection

    // Now after having estimated all compatibilities of all tracks at
    // all vertices, run again over all vertices to set track weights
    // and update the vertex
    setWeightsAndUpdate(state, linearizer);
    if (!state.annealingState.equilibriumReached) {
      m_cfg.annealingTool.anneal(state.annealingState);
    }
    isSmallShift = checkSmallShift(state);

    ++nIter;
  }
  // Multivertex fit is finished

  // Check if smoothing is required
  if (m_cfg.doSmoothing) {
    doVertexSmoothing(state, geoContext);
  }

  return {};
}

template <typename input_track_t, typename linearizer_t>
Acts::Result<void>
Acts::AdaptiveMultiVertexFitter<input_track_t, linearizer_t>::addVtxToFit(
    State& state, Vertex<input_track_t>& newVertex,
    const linearizer_t& linearizer,
    const VertexFitterOptions<input_track_t>& vFitterOptions) const {
  if (state.vtxInfoMap[&newVertex].trackLinks.empty()) {
    return VertexingError::EmptyInput;
  }

  std::vector<Vertex<input_track_t>*> verticesToFit;

  // Prepares vtx and tracks for fast estimation method of their
  // compatibility with vertex
  auto res = prepareVertexForFit(state, &newVertex, vFitterOptions);
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
      const std::vector<const input_track_t*>& trks =
          state.vtxInfoMap[lastVtxIter].trackLinks;
      for (const auto& trk : trks) {
        // Retrieve list of links to all vertices that currently use the current
        // track
        auto range = state.trackToVerticesMultiMap.equal_range(trk);

        // Loop over all attached vertices and add those to vertex fit
        // which are not already in `verticesToFit`
        for (auto vtxIter = range.first; vtxIter != range.second; ++vtxIter) {
          auto newVtxIter = vtxIter->second;
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
  auto fitRes = fitImpl(state, linearizer, vFitterOptions);
  if (!fitRes.ok()) {
    return fitRes.error();
  }

  return {};
}

template <typename input_track_t, typename linearizer_t>
bool Acts::AdaptiveMultiVertexFitter<input_track_t, linearizer_t>::
    isAlreadyInList(
        Vertex<input_track_t>* vtx,
        const std::vector<Vertex<input_track_t>*>& verticesVec) const {
  return std::find(verticesVec.begin(), verticesVec.end(), vtx) !=
         verticesVec.end();
}

template <typename input_track_t, typename linearizer_t>
Acts::Result<void> Acts::
    AdaptiveMultiVertexFitter<input_track_t, linearizer_t>::prepareVertexForFit(
        State& state, Vertex<input_track_t>* vtx,
        const VertexFitterOptions<input_track_t>& vFitterOptions) const {
  // The current vertex info object
  auto& currentVtxInfo = state.vtxInfoMap[vtx];
  // The seed position
  const Vector3D& seedPos = currentVtxInfo.seedPosition.template head<3>();
  auto& geoContext = vFitterOptions.geoContext;

  // Loop over all tracks at current vertex
  for (const auto& trk : currentVtxInfo.trackLinks) {
    auto res = m_cfg.ipEst.getParamsAtClosestApproach(
        geoContext, m_extractParameters(*trk), seedPos);
    if (!res.ok()) {
      return res.error();
    }
    // Set ip3dParams for current trackAtVertex
    currentVtxInfo.ip3dParams.emplace(trk, *(res.value()));
  }
  return {};
}

template <typename input_track_t, typename linearizer_t>
Acts::Result<void>
Acts::AdaptiveMultiVertexFitter<input_track_t, linearizer_t>::
    setAllVertexCompatibilities(State& state, const GeometryContext& geoContext,
                                Vertex<input_track_t>* currentVtx) const {
  VertexInfo<input_track_t>& currentVtxInfo = state.vtxInfoMap[currentVtx];

  // Loop over tracks at current vertex and
  // estimate compatibility with vertex
  for (const auto& trk : currentVtxInfo.trackLinks) {
    auto& trkAtVtx =
        state.tracksAtVerticesMap.at(std::make_pair(trk, currentVtx));
    // Recover from cases where linearization point != 0 but
    // more tracks were added later on
    if (currentVtxInfo.ip3dParams.find(trk) ==
        currentVtxInfo.ip3dParams.end()) {
      auto res = m_cfg.ipEst.getParamsAtClosestApproach(
          geoContext, m_extractParameters(*trk),
          VectorHelpers::position(currentVtxInfo.linPoint));
      if (!res.ok()) {
        return res.error();
      }
      // Set ip3dParams for current trackAtVertex
      currentVtxInfo.ip3dParams.emplace(trk, *(res.value()));
    }
    // Set compatibility with current vertex
    auto compRes = m_cfg.ipEst.getVertexCompatibility(
        geoContext, &(currentVtxInfo.ip3dParams.at(trk)),
        VectorHelpers::position(currentVtxInfo.oldPosition));
    if (!compRes.ok()) {
      return compRes.error();
    }
    trkAtVtx.vertexCompatibility = *compRes;
  }
  return {};
}

template <typename input_track_t, typename linearizer_t>
Acts::Result<void> Acts::AdaptiveMultiVertexFitter<
    input_track_t, linearizer_t>::setWeightsAndUpdate(State& state,
                                                      const linearizer_t&
                                                          linearizer) const {
  for (auto vtx : state.vertexCollection) {
    VertexInfo<input_track_t>& currentVtxInfo = state.vtxInfoMap[vtx];

    for (const auto& trk : currentVtxInfo.trackLinks) {
      auto& trkAtVtx = state.tracksAtVerticesMap.at(std::make_pair(trk, vtx));

      // Set trackWeight for current track
      double currentTrkWeight = m_cfg.annealingTool.getWeight(
          state.annealingState, trkAtVtx.vertexCompatibility,
          collectTrackToVertexCompatibilities(state, trk));
      trkAtVtx.trackWeight = currentTrkWeight;

      if (trkAtVtx.trackWeight > m_cfg.minWeight) {
        // Check if linearization state exists or need to be relinearized
        if (trkAtVtx.linearizedState.covarianceAtPCA ==
                BoundSymMatrix::Zero() ||
            state.vtxInfoMap[vtx].relinearize) {
          auto result = linearizer.linearizeTrack(
              m_extractParameters(*trk), state.vtxInfoMap[vtx].oldPosition);
          if (!result.ok()) {
            return result.error();
          }
          trkAtVtx.linearizedState = *result;
          state.vtxInfoMap[vtx].linPoint = state.vtxInfoMap[vtx].oldPosition;
        }
        // Update the vertex with the new track
        KalmanVertexUpdater::updateVertexWithTrack<input_track_t>(*vtx,
                                                                  trkAtVtx);
      } else {
        ACTS_VERBOSE("Track weight too low. Skip track.");
      }
    }  // End loop over tracks at vertex

    ACTS_VERBOSE("New vertex position: " << vtx->fullPosition());
  }  // End loop over vertex collection

  return {};
}

template <typename input_track_t, typename linearizer_t>
std::vector<double>
Acts::AdaptiveMultiVertexFitter<input_track_t, linearizer_t>::
    collectTrackToVertexCompatibilities(State& state,
                                        const input_track_t* trk) const {
  std::vector<double> trkToVtxCompatibilities;
  trkToVtxCompatibilities.reserve(state.vertexCollection.size());
  auto range = state.trackToVerticesMultiMap.equal_range(trk);

  for (auto vtxIter = range.first; vtxIter != range.second; ++vtxIter) {
    trkToVtxCompatibilities.push_back(
        state.tracksAtVerticesMap.at(std::make_pair(trk, vtxIter->second))
            .vertexCompatibility);
  }

  return trkToVtxCompatibilities;
}

template <typename input_track_t, typename linearizer_t>
bool Acts::AdaptiveMultiVertexFitter<
    input_track_t, linearizer_t>::checkSmallShift(State& state) const {
  for (auto vtx : state.vertexCollection) {
    auto diff = state.vtxInfoMap[vtx].oldPosition.template head<3>() -
                vtx->fullPosition().template head<3>();
    const auto& vtxWgt =
        (vtx->fullCovariance().template block<3, 3>(0, 0)).inverse();
    double relativeShift = diff.dot(vtxWgt * diff);
    if (relativeShift > m_cfg.maxRelativeShift) {
      return false;
    }
  }
  return true;
}

template <typename input_track_t, typename linearizer_t>
void Acts::AdaptiveMultiVertexFitter<input_track_t, linearizer_t>::
    doVertexSmoothing(State& state, const GeometryContext& geoContext) const {
  for (const auto vtx : state.vertexCollection) {
    for (const auto trk : state.vtxInfoMap[vtx].trackLinks) {
      KalmanVertexTrackUpdater::update<input_track_t>(
          geoContext, state.tracksAtVerticesMap.at(std::make_pair(trk, vtx)),
          *vtx);
    }
  }
}
