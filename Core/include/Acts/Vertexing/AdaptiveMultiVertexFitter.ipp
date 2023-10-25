// This file is part of the Acts project.
//
// Copyright (C) 2019-2023 CERN for the benefit of the Acts project
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
    const VertexingOptions<input_track_t>& vertexingOptions) const {
  // Set all vertices to fit in the current state
  state.vertexCollection = verticesToFit;

  // Perform fit on all vertices simultaneously
  auto fitRes = fitImpl(state, linearizer, vertexingOptions);

  if (!fitRes.ok()) {
    return fitRes.error();
  }

  return {};
}

template <typename input_track_t, typename linearizer_t>
Acts::Result<void>
Acts::AdaptiveMultiVertexFitter<input_track_t, linearizer_t>::fitImpl(
    State& state, const linearizer_t& linearizer,
    const VertexingOptions<input_track_t>& vertexingOptions) const {
  // Reset annealing tool
  state.annealingState = AnnealingUtility::State();

  // Boolean indicating whether any of the vertices has moved more than
  // m_cfg.maxRelativeShift during the last iteration. We will keep iterating
  // until the equilibrium (i.e., the lowest temperature) is reached in
  // the annealing procedure and isSmallShift is true (or until the maximum
  // number of iterations is exceeded).
  bool isSmallShift = true;

  // Number of iterations counter
  unsigned int nIter = 0;

  // Start iterating
  while (nIter < m_cfg.maxIterations &&
         (!state.annealingState.equilibriumReached || !isSmallShift)) {
    // Initial loop over all vertices in state.vertexCollection
    for (auto vtx : state.vertexCollection) {
      VertexInfo<input_track_t>& vtxInfo = state.vtxInfoMap[vtx];
      vtxInfo.relinearize = false;
      // Store old position of vertex, i.e. seed position
      // in case of first iteration or position determined
      // in previous iteration afterwards
      vtxInfo.oldPosition = vtx->fullPosition();

      // Calculate the x-y-distance between the current vertex position
      // and the linearization point of the tracks. If it is too large,
      // we relinearize the tracks and recalculate their 3D impact
      // parameters.
      ActsVector<2> xyDiff = vtxInfo.oldPosition.template head<2>() -
                             vtxInfo.linPoint.template head<2>();
      if (xyDiff.norm() > m_cfg.maxDistToLinPoint) {
        // Set flag for relinearization
        vtxInfo.relinearize = true;
      }

      // Check if we use the constraint during the vertex fit
      if (state.vtxInfoMap[vtx].constraint.fullCovariance() !=
          SquareMatrix4::Zero()) {
        const Acts::Vertex<input_track_t>& constraint =
            state.vtxInfoMap[vtx].constraint;
        vtx->setFullPosition(constraint.fullPosition());
        vtx->setFitQuality(constraint.fitQuality());
        vtx->setFullCovariance(constraint.fullCovariance());
      } else if (vtx->fullCovariance() == SquareMatrix4::Zero()) {
        return VertexingError::NoCovariance;
      }

      // Set vertexCompatibility for all TrackAtVertex objects
      // at the current vertex
      setAllVertexCompatibilities(state, vtx, vertexingOptions);
    }  // End loop over vertex collection

    // Recalculate all track weights and update vertices
    setWeightsAndUpdate(state, linearizer, vertexingOptions);

    // Cool the system down, i.e., reduce the temperature parameter. At lower
    // temperatures, outlying tracks are downweighted more.
    if (!state.annealingState.equilibriumReached) {
      m_cfg.annealingTool.anneal(state.annealingState);
    }

    isSmallShift = checkSmallShift(state);
    ++nIter;
  }
  // Multivertex fit is finished

  // Check if smoothing is required
  if (m_cfg.doSmoothing) {
    doVertexSmoothing(state);
  }

  return {};
}

template <typename input_track_t, typename linearizer_t>
Acts::Result<void>
Acts::AdaptiveMultiVertexFitter<input_track_t, linearizer_t>::addVtxToFit(
    State& state, Vertex<input_track_t>& newVertex,
    const linearizer_t& linearizer,
    const VertexingOptions<input_track_t>& vertexingOptions) const {
  if (state.vtxInfoMap[&newVertex].trackLinks.empty()) {
    return VertexingError::EmptyInput;
  }

  std::vector<Vertex<input_track_t>*> verticesToFit;

  // List of vertices added in last iteration
  std::vector<Vertex<input_track_t>*> lastIterAddedVertices = {&newVertex};
  // List of vertices added in current iteration
  std::vector<Vertex<input_track_t>*> currentIterAddedVertices;

  // Fill verticesToFit with vertices that are connected to newVertex (via
  // tracks and/or other vertices).
  while (!lastIterAddedVertices.empty()) {
    for (auto& lastIterAddedVertex : lastIterAddedVertices) {
      // Loop over all tracks at lastIterAddedVertex
      const std::vector<const input_track_t*>& trks =
          state.vtxInfoMap[lastIterAddedVertex].trackLinks;
      for (const auto& trk : trks) {
        // Range of vertices that are associated with trk. The range is
        // represented via its bounds: begin refers to the first iterator of the
        // range; end refers to the iterator after the last iterator of the
        // range.
        auto [begin, end] = state.trackToVerticesMultiMap.equal_range(trk);

        for (auto it = begin; it != end; ++it) {
          // it->first corresponds to trk, it->second to one of its associated
          // vertices
          auto vtxToFit = it->second;
          // Add vertex to the fit if it is not already included
          if (!isAlreadyInList(vtxToFit, verticesToFit)) {
            verticesToFit.push_back(vtxToFit);

            // Collect vertices that were added this iteration
            if (vtxToFit != lastIterAddedVertex) {
              currentIterAddedVertices.push_back(vtxToFit);
            }
          }
        }  // End for loop over range of associated vertices
      }    // End loop over trackLinks
    }      // End loop over lastIterAddedVertices

    lastIterAddedVertices = currentIterAddedVertices;
    currentIterAddedVertices.clear();
  }  // End while loop

  state.vertexCollection = verticesToFit;

  // Perform fit on all added vertices
  auto fitRes = fitImpl(state, linearizer, vertexingOptions);
  if (!fitRes.ok()) {
    return fitRes.error();
  }

  return {};
}

template <typename input_track_t, typename linearizer_t>
bool Acts::AdaptiveMultiVertexFitter<input_track_t, linearizer_t>::
    isAlreadyInList(Vertex<input_track_t>* vtx,
                    const std::vector<Vertex<input_track_t>*>& vertices) const {
  return std::find(vertices.begin(), vertices.end(), vtx) != vertices.end();
}

template <typename input_track_t, typename linearizer_t>
Acts::Result<void>
Acts::AdaptiveMultiVertexFitter<input_track_t, linearizer_t>::
    updateImpactParams3D(
        State& state, Vertex<input_track_t>* vtx,
        const VertexingOptions<input_track_t>& vertexingOptions) const {
  // Vertex info object
  auto& vtxInfo = state.vtxInfoMap[vtx];
  // Vertex position, i.e., point wrt which the impact parameters are estimated
  const Vector3& vtxPosition = vtxInfo.oldPosition.template head<3>();

  // Loop over all tracks at the vertex
  for (const auto& trk : vtxInfo.trackLinks) {
    // Track parameters
    auto trkParams = m_extractParameters(*trk);

    // Origin of the track reference surface
    Vector3 surfaceOrigin = trkParams.referenceSurface()
                                .transform(vertexingOptions.geoContext)
                                .translation();
    // Skip the impact point estimation if the impact parameters of trk wrt the
    // current vertex position were already calculated.
    if (surfaceOrigin == vtxPosition &&
        vtxInfo.impactParams3D.find(trk) != vtxInfo.impactParams3D.end()) {
      continue;
    }

    auto res = m_cfg.ipEst.estimate3DImpactParameters(
        vertexingOptions.geoContext, vertexingOptions.magFieldContext,
        trkParams, vtxPosition, state.ipState);
    if (!res.ok()) {
      return res.error();
    }

    // Try to create a new map entry. If "trk" already exists as key, the
    // corresponding value will not be overwritten and the boolean "inserted"
    // will be set to false ...
    auto [it, inserted] = vtxInfo.impactParams3D.emplace(trk, *res);
    // ... and we have to overwrite manually.
    if (not inserted) {
      it->second = *res;
    }
  }
  return {};
}

template <typename input_track_t, typename linearizer_t>
Acts::Result<void>
Acts::AdaptiveMultiVertexFitter<input_track_t, linearizer_t>::
    setAllVertexCompatibilities(
        State& state, Vertex<input_track_t>* vtx,
        const VertexingOptions<input_track_t>& vertexingOptions) const {
  // Update the 3D impact parameters of all tracks
  updateImpactParams3D(state, vtx, vertexingOptions);

  VertexInfo<input_track_t>& vtxInfo = state.vtxInfoMap[vtx];
  // Loop over the tracks that are associated with vtx and estimate their
  // compatibility
  for (const auto& trk : vtxInfo.trackLinks) {
    auto& trkAtVtx = state.tracksAtVerticesMap.at(std::make_pair(trk, vtx));
    Acts::Result<double> compatibilityResult(0.);
    if (m_cfg.useTime) {
      compatibilityResult = m_cfg.ipEst.template getVertexCompatibility<4>(
          vertexingOptions.geoContext, &(vtxInfo.impactParams3D.at(trk)),
          vtxInfo.oldPosition);
    } else {
      compatibilityResult = m_cfg.ipEst.template getVertexCompatibility<3>(
          vertexingOptions.geoContext, &(vtxInfo.impactParams3D.at(trk)),
          VectorHelpers::position(vtxInfo.oldPosition));
    }
    if (!compatibilityResult.ok()) {
      return compatibilityResult.error();
    }
    trkAtVtx.vertexCompatibility = *compatibilityResult;
  }
  return {};
}

template <typename input_track_t, typename linearizer_t>
Acts::Result<void> Acts::
    AdaptiveMultiVertexFitter<input_track_t, linearizer_t>::setWeightsAndUpdate(
        State& state, const linearizer_t& linearizer,
        const VertexingOptions<input_track_t>& vertexingOptions) const {
  for (auto vtx : state.vertexCollection) {
    VertexInfo<input_track_t>& vtxInfo = state.vtxInfoMap[vtx];

    if (vtxInfo.relinearize) {
      vtxInfo.linPoint = vtxInfo.oldPosition;
    }

    const std::shared_ptr<PerigeeSurface> vtxPerigeeSurface =
        Surface::makeShared<PerigeeSurface>(
            VectorHelpers::position(vtxInfo.linPoint));

    for (const auto& trk : vtxInfo.trackLinks) {
      auto& trkAtVtx = state.tracksAtVerticesMap.at(std::make_pair(trk, vtx));

      // Set trackWeight for current track
      trkAtVtx.trackWeight = m_cfg.annealingTool.getWeight(
          state.annealingState, trkAtVtx.vertexCompatibility,
          collectTrackToVertexCompatibilities(state, trk));

      if (trkAtVtx.trackWeight > m_cfg.minWeight) {
        // Check if track is already linearized and whether we need to
        // relinearize
        if (!trkAtVtx.isLinearized || vtxInfo.relinearize) {
          auto result = linearizer.linearizeTrack(
              m_extractParameters(*trk), vtxInfo.linPoint[3],
              *vtxPerigeeSurface, vertexingOptions.geoContext,
              vertexingOptions.magFieldContext, state.linearizerState);
          if (!result.ok()) {
            return result.error();
          }

          trkAtVtx.linearizedState = *result;
          trkAtVtx.isLinearized = true;
        }
        // Update the vertex with the new track
        KalmanVertexUpdater::updateVertexWithTrack<input_track_t>(*vtx,
                                                                  trkAtVtx);
      } else {
        ACTS_VERBOSE("Track weight too low. Skip track.");
      }
    }  // End loop over tracks at vertex
    ACTS_VERBOSE("New vertex position: " << vtx->fullPosition().transpose());
  }  // End loop over vertex collection

  return {};
}

template <typename input_track_t, typename linearizer_t>
std::vector<double>
Acts::AdaptiveMultiVertexFitter<input_track_t, linearizer_t>::
    collectTrackToVertexCompatibilities(State& state,
                                        const input_track_t* trk) const {
  // Compatibilities of trk wrt all of its associated vertices
  std::vector<double> trkToVtxCompatibilities;

  // Range of vertices that are associated with trk. The range is
  // represented via its bounds: begin refers to the first iterator of the
  // range; end refers to the iterator after the last iterator of the range.
  auto [begin, end] = state.trackToVerticesMultiMap.equal_range(trk);

  // Allocate space in memory for the vector of compatibilities
  trkToVtxCompatibilities.reserve(std::distance(begin, end));

  for (auto it = begin; it != end; ++it) {
    // it->first corresponds to trk, it->second to one of its associated
    // vertices
    trkToVtxCompatibilities.push_back(
        state.tracksAtVerticesMap.at(std::make_pair(trk, it->second))
            .vertexCompatibility);
  }

  return trkToVtxCompatibilities;
}

template <typename input_track_t, typename linearizer_t>
bool Acts::AdaptiveMultiVertexFitter<
    input_track_t, linearizer_t>::checkSmallShift(State& state) const {
  for (auto vtx : state.vertexCollection) {
    Vector3 diff =
        state.vtxInfoMap[vtx].oldPosition.template head<3>() - vtx->position();
    const SquareMatrix3& vtxCov = vtx->covariance();
    double relativeShift = diff.dot(vtxCov.inverse() * diff);
    if (relativeShift > m_cfg.maxRelativeShift) {
      return false;
    }
  }
  return true;
}

template <typename input_track_t, typename linearizer_t>
void Acts::AdaptiveMultiVertexFitter<
    input_track_t, linearizer_t>::doVertexSmoothing(State& state) const {
  for (const auto vtx : state.vertexCollection) {
    for (const auto trk : state.vtxInfoMap[vtx].trackLinks) {
      auto& trkAtVtx = state.tracksAtVerticesMap.at(std::make_pair(trk, vtx));
      if (trkAtVtx.trackWeight > m_cfg.minWeight) {
        KalmanVertexTrackUpdater::update<input_track_t>(trkAtVtx, *vtx);
      }
    }
  }
}
