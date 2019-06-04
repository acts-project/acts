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

  // reset annealing tool
  m_cfg.annealingTool.reset(state.annealingState);

  bool isSmallShift = true;

  // number of iterations counter
  int nIter = 0;

  // start iterating
  while (nIter < m_cfg.maxIterations &&
         (!state.annealingState.equilibriumReached || !isSmallShift)) {
    // initial loop over all vertices in state.vertexCollection
    for (auto& currentVtx : state.vertexCollection) {
      MAVFVertexInfo<input_track_t>& currentVtxInfo =
          state.vtxInfoMap[&currentVtx];
      currentVtxInfo.relinearize = false;

      // TODO: where is this init set? in finder or in else directly below? (->
      // defaulted in Info)
      if (currentVtxInfo.isInitialized) {
        std::cout
            << "Candidate already has a position: storing it in oldpositions"
            << std::endl;
        currentVtxInfo.oldPosition = currentVtx.fullPosition();
      } else {
        // vertex will now be initialized
        currentVtxInfo.isInitialized = true;
        std::cout << "Candidate has no position so far: using as old position "
                     "the seedVertex"
                  << std::endl;
        if (currentVtxInfo.seedPos == SpacePointVector()) {
          std::cout << "seed not set.. something went wrong I guess. Why not "
                       "directly use the seed position from seeds in "
                       "vertexCollection?"
                    << std::endl;
        }
        currentVtxInfo.oldPosition = currentVtxInfo.seedPos;
      }  // end not isInitialized

      if (currentVtxInfo.linPoint == SpacePointVector()) {
        std::cout << "Candidate has no linearization point. where should this "
                     "be set?!"
                  << std::endl;
      }
      if ((currentVtxInfo.oldPosition - currentVtxInfo.linPoint).norm() >
          m_cfg.maxDistToLinPoint) {
        // relinearization needed, distance too big
        currentVtxInfo.relinearize = true;
        // prepare for fit with new vertex position
        prepareVtxForFit(state, currentVtx, vFitterOptions);
      }

      currentVtx.setFullPosition(
          currentVtxInfo.constraintVertex.fullPosition());
      currentVtx.setFullCovariance(
          currentVtxInfo.constraintVertex.fullCovariance() * 1. /
          m_cfg.annealingTool.getWeight(state.annealingState, 1.));
      currentVtx.setFitQuality(currentVtxInfo.constraintVertex.fitQuality());

      // set vertexCompatibility for all TrackAtVertex objects
      // at current vertex
      setAllVtxCompatibilities(state, geoContext, mfContext, currentVtx);
    }  // end loop over vertex collection

    // now after having estimated all compatibilities of all tracks at
    // all vertices, run again over all vertices to compute weights
    setAllTrackWeights(state);

    // TODO: dostuff
    ++nIter;
  }

  return {};
}

template <typename bfield_t, typename input_track_t, typename propagator_t>
Acts::Result<void> Acts::MultiAdaptiveVertexFitter<
    bfield_t, input_track_t,
    propagator_t>::addVertexToFit(State& state, Vertex<input_track_t> newVertex,
                                  const VertexFitterOptions<input_track_t>&
                                      vFitterOptions) const {
  if (newVertex.tracks().empty()) {
    return VertexingError::EmptyInput;
  }

  std::vector<Vertex<input_track_t>*> verticesToFit;

  // prepares vtx and tracks for fast estimation method of their
  // compatibility with vertex
  auto res = prepareVtxForFit(state, newVertex, vFitterOptions);
  if (!res.ok()) {
    return res.error();
  }

  // list of vertices added in last iteration
  std::vector<Vertex<input_track_t>*> lastIterAddedVertices = {&newVertex};
  // list of vertices added in current iteration
  std::vector<Vertex<input_track_t>*> currentIterAddedVertices;

  // loop as long as new vertices are found that share tracks with
  // previously added vertices
  while (!lastIterAddedVertices.empty()) {
    for (auto& lastVtxIter : lastIterAddedVertices) {
      // loop over all track at current lastVtxIter
      for (const TrackAtVertex<input_track_t>& trackIter :
           lastVtxIter->tracks()) {
        // retrieve list of links to all vertices that currently use the current
        // track
        std::vector<Vertex<input_track_t>*>& linksToVertices =
            state.trkInfoMap[&trackIter].linksToVertices;

        // loop over all attached vertices and add those to vertex fit
        // which are not already in `verticesToFit`
        for (auto newVtxIter : linksToVertices) {
          if (!isAlreadyInList(newVtxIter, verticesToFit)) {
            // add newVtxIter to verticesToFit
            verticesToFit.push_back(newVtxIter);

            // add newVtxIter vertex to currentIterAddedVertices
            // if vertex != lastVtxIter
            if (newVtxIter != lastVtxIter) {
              currentIterAddedVertices.push_back(newVtxIter);
            }
          }
        }  // end for loop over linksToVertices
      }
    }  // end loop over lastIterAddedVertices

    lastIterAddedVertices = currentIterAddedVertices;
    currentIterAddedVertices.clear();
  }  // end while loop

  // add vertices to fit to vertexCollection in state
  std::transform(verticesToFit.begin(), verticesToFit.end(),
                 state.vertexCollection.begin(),
                 [](Vertex<input_track_t>* vtx) { return *vtx; });

  // perform fit on all added vertices
  fit(state, vFitterOptions);

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
    propagator_t>::prepareVtxForFit(State& state, Vertex<input_track_t>& vtx,
                                    const VertexFitterOptions<input_track_t>&
                                        vFitterOptions) const {
  const Vector3D& refPos = vtx.position();
  auto& geoContext = vFitterOptions.geoContext;
  auto& mfContext = vFitterOptions.magFieldContext;

  // loop over all tracks at current vertex
  for (const auto& trkAtVtx : vtx.tracks()) {
    auto res = m_cfg.ipEst.getParamsAtIP3d(
        geoContext, mfContext, m_extractParameters(trkAtVtx.originalTrack),
        refPos);
    if (!res.ok()) {
      return res.error();
    }
    // set ip3dParams for current trackAtVertex
    state.trkInfoMap[&trkAtVtx].ip3dParams = std::move(res.value());
  }
  return {};
}

template <typename bfield_t, typename input_track_t, typename propagator_t>
Acts::Result<void>
Acts::MultiAdaptiveVertexFitter<bfield_t, input_track_t, propagator_t>::
    setAllVtxCompatibilities(State& state, const GeometryContext& geoContext,
                             const MagneticFieldContext& mfContext,
                             Vertex<input_track_t>& currentVtx) const {
  MAVFVertexInfo<input_track_t>& currentVtxInfo = state.vtxInfoMap[&currentVtx];
  // create empty list of new TrackAtVertex objects
  // to be filled below. Needed due to constness of
  // tracksAtVertex list at vertex
  std::vector<TrackAtVertex<input_track_t>> newTracks;
  newTracks.reserve(currentVtx.tracks().size());

  // loop over tracks at current vertex and
  // estimate compatibility with vertex
  for (auto& trkAtVtx : currentVtx.tracks()) {
    // TODO: DO I NEED THE ipEst.getParams stuff below?
    // recover from cases where linearization point != 0 but
    // more tracks were added later on
    if (!state.trkInfoMap[&trkAtVtx].ip3dParams) {
      auto res = m_cfg.ipEst.getParamsAtIP3d(
          geoContext, mfContext, m_extractParameters(trkAtVtx.originalTrack),
          VectorHelpers::position(currentVtxInfo.linPoint));
      if (!res.ok()) {
        return res.error();
      }
      // set ip3dParams for current trackAtVertex
      state.trkInfoMap[&trkAtVtx].ip3dParams = std::move(res.value());
    }

    // create copy of current trackAtVertex in order
    // to modify it below
    TrackAtVertex<input_track_t> newTrkAtVtx = trkAtVtx;

    // update trkInfoMap accordingly
    state.trkInfoMap[&newTrkAtVtx] = std::move(state.trkInfoMap[&trkAtVtx]);
    state.trkInfoMap.erase(&trkAtVtx);

    // set compatibility with current vertex
    // TODO: is that where the ip3d step beforehand is needed for?
    // not using any of these here
    m_cfg.trackCompEst.setTrackCompatibility(
        geoContext, newTrkAtVtx,
        VectorHelpers::position(currentVtxInfo.oldPosition),
        m_extractParameters);

    newTracks.push_back(newTrkAtVtx);
  }
  // set list of updated tracks to current vertex
  currentVtx.setTracksAtVertex(newTracks);
  return {};
}

template <typename bfield_t, typename input_track_t, typename propagator_t>
Acts::Result<void> Acts::MultiAdaptiveVertexFitter<
    bfield_t, input_track_t, propagator_t>::setAllTrackWeights(State& state)
    const {
  for (auto& vtx : state.vertexCollection) {
    // create empty list of new TrackAtVertex objects
    // to be filled below. Needed due to constness of
    // tracksAtVertex list at vertex
    std::vector<TrackAtVertex<input_track_t>> newTracks;
    newTracks.reserve(vtx.tracks().size());

    for (const auto& trkAtVtx : vtx.tracks()) {
      // create copy of current trackAtVertex in order
      // to modify it below
      TrackAtVertex<input_track_t> newTrkAtVtx = trkAtVtx;

      auto collectRes = collectTrkToVtxCompatibilities(state, trkAtVtx);
      if (!collectRes.ok()) {
        return collectRes.error();
      }

      // set vertexCompatibility for current track
      newTrkAtVtx.vertexCompatibility = m_cfg.annealingTool.getWeight(
          state.annealingState, trkAtVtx.vertexCompatibility, *collectRes);

      // update trkInfoMap accordingly
      state.trkInfoMap[&newTrkAtVtx] = std::move(state.trkInfoMap[&trkAtVtx]);
      state.trkInfoMap.erase(&trkAtVtx);

      newTracks.push_back(newTrkAtVtx);
    }  // end loop over tracks at vertex

    // update tracks at current vertex
    vtx.setTracksAtVertex(newTracks);
  }  // end loop over vertex collection

  return {};
}

template <typename bfield_t, typename input_track_t, typename propagator_t>
Acts::Result<std::vector<double>>
Acts::MultiAdaptiveVertexFitter<bfield_t, input_track_t, propagator_t>::
    collectTrkToVtxCompatibilities(
        State& state, const TrackAtVertex<input_track_t>& trk) const {
  // All vertices that currently hold the track `trk`
  std::vector<Vertex<input_track_t>*> vertices =
      state.trkInfoMap[&trk].linksToVertices;

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
