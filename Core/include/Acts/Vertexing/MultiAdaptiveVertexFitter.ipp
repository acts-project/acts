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
  // reset annealing tool
  m_cfg.annealingTool.reset(state.annealingState);

  // create some maps to store important information
  std::map<Vertex<input_track_t>*, MAVFVertexInfo<input_track_t>> infoMap;

  bool isSmallShift = true;

  // number of iterations counter
  int nIter = 0;

  // start iterating
  while (nIter < m_cfg.maxIterations &&
         (!state.annealingState.equilibriumReached || !isSmallShift)) {
    // initial loop over all vertices in state.vertexCollection
    for (auto& currentVtx : state.vertexCollection) {
      MAVFVertexInfo<input_track_t>& currentFitInfo = infoMap[*currentVtx];
      currentFitInfo.relinearize = false;

      // TODO: where is this init set? in finder or in else directly below? (->
      // defaulted in Info)
      if (currentFitInfo.isInitialized) {
        std::cout
            << "Candidate already has a position: storing it in oldpositions"
            << std::endl;
        currentFitInfo.oldPosition = currentVtx.position();
      } else {
        // vertex will now be initialized
        currentFitInfo.isInitialized = true;
        std::cout << "Candidate has no position so far: using as old position "
                     "the seedVertex"
                  << std::endl;
        if (currentFitInfo.seedPos == Vector3D()) {
          std::cout << "seed not set.. something went wrong I guess. Why not "
                       "directly use the seed position from seeds in "
                       "vertexCollection?"
                    << std::endl;
        }
        currentFitInfo.oldPosition = currentFitInfo.seedPos;
      }  // end not isInitialized

      if (currentFitInfo.linPoint == Vector3D()) {
        std::cout << "Candidate has no linearization point. where should this "
                     "be set?!"
                  << std::endl;
      }
      if ((currentFitInfo.oldPosition - currentFitInfo.linPoint).perp() >
          m_cfg.maxDistToLinPoint) {
        // relinearization needed, distance too big
        currentFitInfo.relinearize = true;
      }
    }

    // TODO: dostuff
    ++nIter;
  }

  return {};
}

template <typename bfield_t, typename input_track_t, typename propagator_t>
Acts::Result<void>
Acts::MultiAdaptiveVertexFitter<bfield_t, input_track_t, propagator_t>::
    addVertexToFit(
        State& state, Vertex<input_track_t> newVertex,
        // TODO: pass by const ref?
        MAVFTrackAtVtxInfo<input_track_t> trackAtVtxInfo,
        const VertexFitterOptions<input_track_t>& vFitterOptions) const {
  if (newVertex.tracks().empty()) {
    return VertexingError::EmptyInput;
  }

  std::vector<Vertex<input_track_t>*> verticesToFit;

  // prepares tracks for fast estimation method of their
  // compatibility with vertex
  // TODO: uncomment
  // prepareCompatibility(newVertex);

  // list of vertices added in last iteration
  std::vector<Vertex<input_track_t>*> lastIterAddedVertices = {&newVertex};
  // list of vertices added in current iteration
  std::vector<Vertex<input_track_t>*> currentIterAddedVertices;

  // loop as long as new vertices are found that share tracks with
  // previously added vertices
  while (!lastIterAddedVertices.empty()) {
    for (auto lastVtxIter : lastIterAddedVertices) {
      // loop over all track at current lastVtxIter
      for (auto& trackIter : lastVtxIter->tracks()) {
        // retrieve list of links to all vertices that currently use the current
        // track
        std::vector<Vertex<input_track_t>*>& linksToVertices =
            trackAtVtxInfo.linkMapToVertices[&trackIter];

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
                      [vtx](Vertex<input_track_t>* v) { return v == vtx; });
}
