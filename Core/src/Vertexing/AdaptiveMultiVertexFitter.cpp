// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/AdaptiveMultiVertexFitter.hpp"

#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Vertexing/KalmanVertexUpdater.hpp"
#include "Acts/Vertexing/VertexingError.hpp"

namespace Acts {

AdaptiveMultiVertexFitter::AdaptiveMultiVertexFitter(
    Config cfg, std::unique_ptr<const Logger> logger)
    : m_cfg(std::move(cfg)), m_logger(std::move(logger)) {
  if (!m_cfg.extractParameters.connected()) {
    throw std::invalid_argument(
        "AdaptiveMultiVertexFitter: No function to extract parameters "
        "from InputTrack_t provided.");
  }

  if (!m_cfg.trackLinearizer.connected()) {
    throw std::invalid_argument(
        "AdaptiveMultiVertexFitter: No track linearizer provided.");
  }
}

VertexScratch& AdaptiveMultiVertexFitter::scratchFor(
    const VertexFitProblem& problem, Vertex* vtx, Cache& cache) {
  auto [it, inserted] = cache.vertexScratch.try_emplace(vtx);
  if (inserted) {
    // Seed the linearization point and the previous position with the seed
    // position of the vertex, so that the first iteration measures the
    // distance to the seed rather than to the origin. A vertex without a
    // candidate entry keeps the zero-initialised default, matching the
    // previously default-constructed VertexInfo.
    if (auto candidateIt = problem.candidates.find(vtx);
        candidateIt != problem.candidates.end()) {
      it->second.linPoint = candidateIt->second.seedPosition;
      it->second.oldPosition = candidateIt->second.seedPosition;
    }
  }
  return it->second;
}

Result<void> AdaptiveMultiVertexFitter::fit(
    VertexFitProblem& problem, const VertexingOptions& vertexingOptions,
    Cache& cache) const {
  // Reset annealing tool
  cache.annealingState = AnnealingUtility::State();

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
         (!cache.annealingState.equilibriumReached || !isSmallShift)) {
    // Initial loop over all vertices in problem.vertices
    for (auto vtx : problem.vertices) {
      VertexScratch& scratch = scratchFor(problem, vtx, cache);
      scratch.relinearize = false;
      // Store old position of vertex, i.e. seed position
      // in case of first iteration or position determined
      // in previous iteration afterwards
      scratch.oldPosition = vtx->fullPosition();

      // Calculate the x-y-distance between the current vertex position
      // and the linearization point of the tracks. If it is too large,
      // we relinearize the tracks and recalculate their 3D impact
      // parameters.
      Vector2 xyDiff = scratch.oldPosition.template head<2>() -
                       scratch.linPoint.template head<2>();
      if (xyDiff.norm() > m_cfg.maxDistToLinPoint) {
        // Set flag for relinearization
        scratch.relinearize = true;
        // Recalculate the track impact parameters at the current vertex
        // position
        auto prepareVertexResult =
            prepareVertexForFit(problem, vtx, vertexingOptions, cache);
        if (!prepareVertexResult.ok()) {
          // Print vertices and associated tracks if logger is in debug mode
          if (logger().doPrint(Logging::DEBUG)) {
            logDebugData(problem, vertexingOptions.geoContext, cache);
          }
          return prepareVertexResult.error();
        }
      }

      // Check if we use the constraint during the vertex fit
      if (problem.candidates[vtx].constraint.fullCovariance() !=
          SquareMatrix4::Zero()) {
        const Vertex& constraint = problem.candidates[vtx].constraint;
        vtx->setFullPosition(constraint.fullPosition());
        vtx->setFitQuality(constraint.fitQuality());
        vtx->setFullCovariance(constraint.fullCovariance());
      } else if (vtx->fullCovariance() == SquareMatrix4::Zero()) {
        return VertexingError::NoCovariance;
      }

      // Set vertexCompatibility for all TrackAtVertex objects
      // at the current vertex
      auto setCompatibilitiesResult =
          setAllVertexCompatibilities(problem, vtx, vertexingOptions, cache);
      if (!setCompatibilitiesResult.ok()) {
        // Print vertices and associated tracks if logger is in debug mode
        if (logger().doPrint(Logging::DEBUG)) {
          logDebugData(problem, vertexingOptions.geoContext, cache);
        }
        return setCompatibilitiesResult.error();
      }
    }  // End loop over vertex collection

    // Recalculate all track weights and update vertices
    auto setWeightsResult =
        setWeightsAndUpdate(problem, vertexingOptions, cache);
    if (!setWeightsResult.ok()) {
      // Print vertices and associated tracks if logger is in debug mode
      if (logger().doPrint(Logging::DEBUG)) {
        logDebugData(problem, vertexingOptions.geoContext, cache);
      }
      return setWeightsResult.error();
    }

    // Cool the system down, i.e., reduce the temperature parameter. At lower
    // temperatures, outlying tracks are downweighted more.
    if (!cache.annealingState.equilibriumReached) {
      m_cfg.annealingTool.anneal(cache.annealingState);
    }

    isSmallShift = checkSmallShift(problem, cache);
    ++nIter;
  }
  // Multivertex fit is finished

  // Check if smoothing is required
  if (m_cfg.doSmoothing) {
    doVertexSmoothing(problem);
  }

  return {};
}

Result<void> AdaptiveMultiVertexFitter::addVtxToFit(
    VertexFitProblem& problem, const std::vector<Vertex*>& newVertices,
    const VertexingOptions& vertexingOptions, Cache& cache) const {
  for (const auto& newVertex : newVertices) {
    if (problem.candidates[newVertex].trackLinks.empty()) {
      ACTS_ERROR(
          "newVertex does not have any associated tracks (i.e., its trackLinks "
          "are empty).");
      return VertexingError::EmptyInput;
    }
  }

  std::vector<Vertex*> verticesToFit = newVertices;

  // List of vertices added in last iteration
  std::vector<Vertex*> lastIterAddedVertices = newVertices;
  // List of vertices added in current iteration
  std::vector<Vertex*> currentIterAddedVertices;

  // Fill verticesToFit with vertices that are connected to newVertex (via
  // tracks and/or other vertices).
  while (!lastIterAddedVertices.empty()) {
    for (auto& lastIterAddedVertex : lastIterAddedVertices) {
      // Loop over all tracks at lastIterAddedVertex
      const std::vector<InputTrack>& trks =
          problem.candidates[lastIterAddedVertex].trackLinks;
      for (const auto& trk : trks) {
        // Range of vertices that are associated with trk. The range is
        // represented via its bounds: begin refers to the first iterator of the
        // range; end refers to the iterator after the last iterator of the
        // range.
        auto [begin, end] = problem.trackToVertices.equal_range(trk);

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
      }  // End loop over trackLinks
    }  // End loop over lastIterAddedVertices

    lastIterAddedVertices = currentIterAddedVertices;
    currentIterAddedVertices.clear();
  }  // End while loop

  problem.vertices = verticesToFit;

  for (Vertex* newVertexPtr : newVertices) {
    // Save the 3D impact parameters of all tracks associated with newVertex.
    auto res =
        prepareVertexForFit(problem, newVertexPtr, vertexingOptions, cache);
    if (!res.ok()) {
      // Print vertices and associated tracks if logger is in debug mode
      if (logger().doPrint(Logging::DEBUG)) {
        logDebugData(problem, vertexingOptions.geoContext, cache);
      }
      return res.error();
    }
  }

  // Perform fit on all added vertices
  auto fitRes = fit(problem, vertexingOptions, cache);
  if (!fitRes.ok()) {
    return fitRes.error();
  }

  return {};
}

bool AdaptiveMultiVertexFitter::isAlreadyInList(
    Vertex* vtx, const std::vector<Vertex*>& vertices) const {
  return rangeContainsValue(vertices, vtx);
}

Result<void> AdaptiveMultiVertexFitter::prepareVertexForFit(
    const VertexFitProblem& problem, Vertex* vtx,
    const VertexingOptions& vertexingOptions, Cache& cache) const {
  // Vertex seed position
  const Vector3& seedPos =
      problem.candidates.at(vtx).seedPosition.template head<3>();
  VertexScratch& scratch = scratchFor(problem, vtx, cache);

  // Loop over all tracks at the vertex
  for (const auto& trk : problem.candidates.at(vtx).trackLinks) {
    auto res = m_cfg.ipEst.estimate3DImpactParameters(
        vertexingOptions.geoContext, vertexingOptions.magFieldContext,
        m_cfg.extractParameters(trk), seedPos, cache.ipState);
    if (!res.ok()) {
      return res.error();
    }
    // Save 3D impact parameters of the track
    scratch.impactParams3D.emplace(trk, res.value());
  }
  return {};
}

Result<void> AdaptiveMultiVertexFitter::setAllVertexCompatibilities(
    VertexFitProblem& problem, Vertex* vtx,
    const VertexingOptions& vertexingOptions, Cache& cache) const {
  VertexFitCandidate& candidate = problem.candidates[vtx];
  VertexScratch& scratch = scratchFor(problem, vtx, cache);

  // Loop over all tracks that are associated with vtx and estimate their
  // compatibility
  for (const auto& trk : candidate.trackLinks) {
    auto& trkAtVtx = problem.tracksAtVertices.at(std::make_pair(trk, vtx));
    // Recover from cases where linearization point != 0 but
    // more tracks were added later on
    if (!scratch.impactParams3D.contains(trk)) {
      auto res = m_cfg.ipEst.estimate3DImpactParameters(
          vertexingOptions.geoContext, vertexingOptions.magFieldContext,
          m_cfg.extractParameters(trk),
          VectorHelpers::position(scratch.linPoint), cache.ipState);
      if (!res.ok()) {
        return res.error();
      }
      // Set impactParams3D for current trackAtVertex
      scratch.impactParams3D.emplace(trk, res.value());
    }
    // Set compatibility with current vertex
    Result<double> compatibilityResult(0.);
    if (m_cfg.useTime) {
      compatibilityResult = m_cfg.ipEst.getVertexCompatibility(
          vertexingOptions.geoContext, &(scratch.impactParams3D.at(trk)),
          scratch.oldPosition);
    } else {
      Vector3 vertexPosOnly = VectorHelpers::position(scratch.oldPosition);
      compatibilityResult = m_cfg.ipEst.getVertexCompatibility(
          vertexingOptions.geoContext, &(scratch.impactParams3D.at(trk)),
          vertexPosOnly);
    }

    if (!compatibilityResult.ok()) {
      return compatibilityResult.error();
    }
    trkAtVtx.vertexCompatibility = *compatibilityResult;
  }
  return {};
}

Result<void> AdaptiveMultiVertexFitter::setWeightsAndUpdate(
    VertexFitProblem& problem, const VertexingOptions& vertexingOptions,
    Cache& cache) const {
  for (auto vtx : problem.vertices) {
    VertexScratch& scratch = scratchFor(problem, vtx, cache);

    if (scratch.relinearize) {
      scratch.linPoint = scratch.oldPosition;
    }

    const std::shared_ptr<PerigeeSurface> vtxPerigeeSurface =
        Surface::makeShared<PerigeeSurface>(
            VectorHelpers::position(scratch.linPoint));

    for (const auto& trk : problem.candidates[vtx].trackLinks) {
      auto& trkAtVtx = problem.tracksAtVertices.at(std::make_pair(trk, vtx));

      // Set trackWeight for current track
      trkAtVtx.trackWeight = m_cfg.annealingTool.getWeight(
          cache.annealingState, trkAtVtx.vertexCompatibility,
          collectTrackToVertexCompatibilities(problem, trk));

      if (trkAtVtx.trackWeight > m_cfg.minWeight) {
        // Check if track is already linearized and whether we need to
        // relinearize
        if (!trkAtVtx.isLinearized || scratch.relinearize) {
          auto result = m_cfg.trackLinearizer(
              m_cfg.extractParameters(trk), scratch.linPoint[3],
              *vtxPerigeeSurface, vertexingOptions.geoContext,
              vertexingOptions.magFieldContext, cache.fieldCache);
          if (!result.ok()) {
            return result.error();
          }

          trkAtVtx.linearizedState = *result;
          trkAtVtx.isLinearized = true;
        }
        // Update the vertex with the new track. The second template
        // argument corresponds to the number of fitted vertex dimensions
        // (i.e., 3 if we only fit spatial coordinates and 4 if we also fit
        // time).
        KalmanVertexUpdater::updateVertexWithTrack(*vtx, trkAtVtx,
                                                   m_cfg.useTime ? 4 : 3);
      } else {
        ACTS_VERBOSE("Track weight too low. Skip track.");
      }
    }  // End loop over tracks at vertex
    ACTS_VERBOSE("New vertex position: " << vtx->fullPosition().transpose());
  }  // End loop over vertex collection

  return {};
}

std::vector<double>
AdaptiveMultiVertexFitter::collectTrackToVertexCompatibilities(
    const VertexFitProblem& problem, const InputTrack& trk) const {
  // Compatibilities of trk wrt all of its associated vertices
  std::vector<double> trkToVtxCompatibilities;

  // Range of vertices that are associated with trk. The range is
  // represented via its bounds: begin refers to the first iterator of the
  // range; end refers to the iterator after the last iterator of the range.
  auto [begin, end] = problem.trackToVertices.equal_range(trk);
  // Allocate space in memory for the vector of compatibilities
  trkToVtxCompatibilities.reserve(std::distance(begin, end));

  for (auto it = begin; it != end; ++it) {
    // it->first corresponds to trk, it->second to one of its associated
    // vertices
    trkToVtxCompatibilities.push_back(
        problem.tracksAtVertices.at(std::make_pair(trk, it->second))
            .vertexCompatibility);
  }

  return trkToVtxCompatibilities;
}

bool AdaptiveMultiVertexFitter::checkSmallShift(const VertexFitProblem& problem,
                                                Cache& cache) const {
  for (auto* vtx : problem.vertices) {
    Vector3 diff =
        scratchFor(problem, vtx, cache).oldPosition.template head<3>() -
        vtx->position();
    const SquareMatrix3& vtxCov = vtx->covariance();
    double relativeShift = diff.dot(vtxCov.inverse() * diff);
    if (relativeShift > m_cfg.maxRelativeShift) {
      return false;
    }
  }
  return true;
}

void AdaptiveMultiVertexFitter::doVertexSmoothing(
    VertexFitProblem& problem) const {
  for (const auto vtx : problem.vertices) {
    for (const auto& trk : problem.candidates[vtx].trackLinks) {
      auto& trkAtVtx = problem.tracksAtVertices.at(std::make_pair(trk, vtx));
      if (trkAtVtx.trackWeight > m_cfg.minWeight) {
        // Update the new track under the assumption that it originates at the
        // vertex. The second template argument corresponds to the number of
        // fitted vertex dimensions (i.e., 3 if we only fit spatial coordinates
        // and 4 if we also fit time).
        KalmanVertexUpdater::updateTrackWithVertex(trkAtVtx, *vtx,
                                                   m_cfg.useTime ? 4 : 3);
      }
    }
  }
}

void AdaptiveMultiVertexFitter::logDebugData(const VertexFitProblem& problem,
                                             const GeometryContext& geoContext,
                                             Cache& cache) const {
  ACTS_DEBUG("Encountered an error when fitting the following "
             << problem.vertices.size() << " vertices:");
  for (std::size_t vtxInd = 0; vtxInd < problem.vertices.size(); ++vtxInd) {
    auto vtx = problem.vertices[vtxInd];
    ACTS_DEBUG("Position of " << vtxInd << ". vertex seed:\n"
                              << problem.candidates.at(vtx).seedPosition);
    ACTS_DEBUG("Position of said vertex after the last fitting step:\n"
               << scratchFor(problem, vtx, cache).oldPosition);
    ACTS_DEBUG("Associated tracks:");
    const auto& trks = problem.candidates.at(vtx).trackLinks;
    for (std::size_t trkInd = 0; trkInd < trks.size(); ++trkInd) {
      const auto& trkAtVtx =
          problem.tracksAtVertices.at(std::make_pair(trks[trkInd], vtx));
      const auto& trkParams = m_cfg.extractParameters(trkAtVtx.originalParams);
      ACTS_DEBUG(trkInd << ". track parameters:\n" << trkParams.parameters());
      ACTS_DEBUG(trkInd << ". track covariance matrix:\n"
                        << trkParams.covariance().value());
      ACTS_DEBUG("Origin of corresponding reference surface:\n"
                 << trkParams.referenceSurface().center(geoContext));
    }
  }
}

}  // namespace Acts
