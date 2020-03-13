// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/VertexingError.hpp"

template <typename vfitter_t, typename sfinder_t>
auto Acts::AdaptiveMultiVertexFinder<vfitter_t, sfinder_t>::find(
    const std::vector<const InputTrack_t*>& allTracks,
    const VertexFinderOptions<InputTrack_t>& vFinderOptions) const
    -> Result<std::vector<Vertex<InputTrack_t>>> {
  if (allTracks.empty()) {
    return VertexingError::EmptyInput;
  }
  // Original tracks
  const std::vector<const InputTrack_t*>& origTracks = allTracks;

  // Seed tracks
  std::vector<const InputTrack_t*> seedTracks = allTracks;

  // Construct the vertex fitter options from vertex finder options
  VertexFitterOptions<InputTrack_t> vFitterOptions(
      vFinderOptions.geoContext, vFinderOptions.magFieldContext,
      vFinderOptions.vertexConstraint);

  FitterState_t fitterState;

  std::vector<std::unique_ptr<Vertex<InputTrack_t>>> allVertices;

  std::vector<Vertex<InputTrack_t>*> allVerticesPtr;

  int iteration = 0;
  while (((m_cfg.addSingleTrackVertices && seedTracks.size() > 0) ||
          ((!m_cfg.addSingleTrackVertices) && seedTracks.size() > 1)) &&
         iteration < m_cfg.maxIterations) {
    FitterState_t oldFitterState = fitterState;

    // Tracks that are used for searching compatible tracks
    // near a vertex candidate
    std::vector<const InputTrack_t*> allTracks;
    if (m_cfg.doRealMultiVertex) {
      allTracks = origTracks;
    } else {
      allTracks = seedTracks;
    }
    Vertex<InputTrack_t> currentConstraint = vFinderOptions.vertexConstraint;
    // Retrieve seed vertex from all remaining seedTracks
    auto seedResult = doSeeding(seedTracks, currentConstraint, vFinderOptions);
    if (!seedResult.ok()) {
      return seedResult.error();
    }
    allVertices.push_back(std::make_unique<Vertex<InputTrack_t>>(*seedResult));

    Vertex<InputTrack_t>& vtxCandidate = *allVertices.back();
    allVerticesPtr.push_back(&vtxCandidate);

    ACTS_DEBUG("Position of current vertex candidate after seeding: "
               << vtxCandidate.fullPosition());
    if (vtxCandidate.position().z() == 0.) {
      ACTS_DEBUG(
          "No seed found anymore. Break and stop primary vertex finding.");
      allVertices.pop_back();
      allVerticesPtr.pop_back();
      break;
    }
    auto prepResult = canPrepareVertexForFit(
        allTracks, seedTracks, vtxCandidate, currentConstraint, fitterState);

    if (!prepResult.ok()) {
      return prepResult.error();
    }
    if (!(*prepResult)) {
      ACTS_DEBUG("Could not prepare for fit anymore. Break.");
      allVertices.pop_back();
      allVerticesPtr.pop_back();
      break;
    }
    // Update fitter state with all vertices
    fitterState.addVertexToMultiMap(vtxCandidate);

    // Perform the fit
    auto fitResult = m_cfg.vertexFitter.addVtxToFit(
        fitterState, vtxCandidate, m_cfg.linearizer, vFitterOptions);
    if (!fitResult.ok()) {
      return fitResult.error();
    }
    ACTS_DEBUG("New position of current vertex candidate after fit: "
               << vtxCandidate.fullPosition());
    // Check if vertex is good vertex
    auto [nCompatibleTracks, isGoodVertex] =
        checkVertexAndCompatibleTracks(vtxCandidate, seedTracks, fitterState);

    ACTS_DEBUG("Vertex is good vertex: " << isGoodVertex);
    if (nCompatibleTracks > 0) {
      removeCompatibleTracksFromSeedTracks(vtxCandidate, seedTracks,
                                           fitterState);
    } else {
      bool removedNonCompatibleTrack =
          canRemoveNonCompatibleTrackFromSeedTracks(vtxCandidate, seedTracks,
                                                    fitterState);
      if (!removedNonCompatibleTrack) {
        ACTS_DEBUG(
            "Could not remove any further track from seed tracks. Break.");
        allVertices.pop_back();
        allVerticesPtr.pop_back();
        break;
      }
    }
    bool keepVertex = isGoodVertex &&
                      keepNewVertex(vtxCandidate, allVerticesPtr, fitterState);
    ACTS_DEBUG("New vertex will be saved: " << keepVertex);

    // Delete vertex from allVertices list again if it's not kept
    if (not keepVertex) {
      auto deleteVertexResult =
          deleteLastVertex(vtxCandidate, allVertices, allVerticesPtr,
                           fitterState, oldFitterState, vFitterOptions);
      if (not deleteVertexResult.ok()) {
        return deleteVertexResult.error();
      }
    }
    iteration++;
  }  // end while loop

  return getVertexOutputList(allVerticesPtr, fitterState);
}

template <typename vfitter_t, typename sfinder_t>
auto Acts::AdaptiveMultiVertexFinder<vfitter_t, sfinder_t>::doSeeding(
    const std::vector<const InputTrack_t*>& trackVector,
    Vertex<InputTrack_t>& currentConstraint,
    const VertexFinderOptions<InputTrack_t>& vFinderOptions) const
    -> Result<Vertex<InputTrack_t>> {
  VertexFinderOptions<InputTrack_t> seedOptions = vFinderOptions;
  seedOptions.vertexConstraint = currentConstraint;
  // Run seed finder
  auto seedResult = m_cfg.seedFinder.find(trackVector, seedOptions);

  if (!seedResult.ok()) {
    return seedResult.error();
  }

  Vertex<InputTrack_t> seedVertex = (*seedResult).back();
  // Update constraints according to seed vertex
  if (m_cfg.useBeamSpotConstraint) {
    if (currentConstraint.fullCovariance() == SpacePointSymMatrix::Zero()) {
      ACTS_WARNING(
          "No constraint provided, but useBeamSpotConstraint set to true.");
    }
    if (m_cfg.useSeedConstraint) {
      currentConstraint.setFullPosition(seedVertex.fullPosition());
      currentConstraint.setFullCovariance(seedVertex.fullCovariance());
    }
  } else {
    currentConstraint.setFullPosition(seedVertex.fullPosition());
    currentConstraint.setFullCovariance(SpacePointSymMatrix::Identity() *
                                        m_cfg.looseConstrValue);
    currentConstraint.setFitQuality(m_cfg.defaultConstrFitQuality);
  }

  return seedVertex;
}

template <typename vfitter_t, typename sfinder_t>
auto Acts::AdaptiveMultiVertexFinder<vfitter_t, sfinder_t>::estimateDeltaZ(
    const BoundParameters& track, const Vector3D& vtxPos) const -> double {
  Vector3D trackPos = track.position();

  double phi = track.parameters()[ParID_t::ePHI];
  double th = track.parameters()[ParID_t::eTHETA];

  double X = trackPos[eX] - vtxPos.x();
  double Y = trackPos[eY] - vtxPos.y();

  double deltaZ = trackPos[eZ] - vtxPos.z() -
                  1. / std::tan(th) * (X * std::cos(phi) + Y * std::sin(phi));

  return deltaZ;
}

template <typename vfitter_t, typename sfinder_t>
auto Acts::AdaptiveMultiVertexFinder<vfitter_t, sfinder_t>::getIPSignificance(
    const InputTrack_t* track, const Vertex<InputTrack_t>& vtx) const
    -> Result<double> {
  // TODO: In original implementation the covariance of the given vertex is set
  // to zero. I did the same here now, but consider removing this and just
  // passing the vtx object to the estimator without changing its covariance.
  // After all, the vertex seed does have a non-zero convariance in general and
  // it probably should be used.
  Vertex<InputTrack_t> newVtx = vtx;
  if (not m_cfg.useVertexCovForIPEstimation) {
    newVtx.setFullCovariance(SpacePointSymMatrix::Zero());
  }

  auto estRes = m_cfg.ipEstimator.estimate(*track, newVtx);
  if (!estRes.ok()) {
    return estRes.error();
  }

  ImpactParametersAndSigma ipas = *estRes;

  double significance = 0.;
  if (ipas.sigmad0 > 0 && ipas.sigmaz0 > 0) {
    significance = std::sqrt(std::pow(ipas.IPd0 / ipas.sigmad0, 2) +
                             std::pow(ipas.IPz0 / ipas.sigmaz0, 2));
  }

  return significance;
}

template <typename vfitter_t, typename sfinder_t>
auto Acts::AdaptiveMultiVertexFinder<vfitter_t, sfinder_t>::
    addCompatibleTracksToVertex(const std::vector<const InputTrack_t*>& tracks,
                                Vertex<InputTrack_t>& vtx,
                                FitterState_t& fitterState) const
    -> Result<void> {
  for (const auto& trk : tracks) {
    auto sigRes = getIPSignificance(trk, vtx);
    if (!sigRes.ok()) {
      return sigRes.error();
    }
    double ipSig = *sigRes;
    auto params = m_extractParameters(*trk);
    if ((std::abs(estimateDeltaZ(params, vtx.position())) <
         m_cfg.tracksMaxZinterval) &&
        (ipSig < m_cfg.tracksMaxSignificance)) {
      // Create TrackAtVertex objects, unique for each (track, vertex) pair
      // fitterState.tracksAtVerticesMap.clear();
      fitterState.tracksAtVerticesMap.emplace(std::make_pair(trk, &vtx),
                                              TrackAtVertex(params, trk));

      // Add the original track parameters to the list for vtx
      fitterState.vtxInfoMap[&vtx].trackLinks.push_back(trk);
    }
  }
  return {};
}

template <typename vfitter_t, typename sfinder_t>
auto Acts::AdaptiveMultiVertexFinder<vfitter_t, sfinder_t>::
    canRecoverFromNoCompatibleTracks(
        const std::vector<const InputTrack_t*>& allTracks,
        const std::vector<const InputTrack_t*>& seedTracks,
        Vertex<InputTrack_t>& vtx,
        const Vertex<InputTrack_t>& currentConstraint,
        FitterState_t& fitterState) const -> Result<bool> {
  // Recover from cases where no compatible tracks to vertex
  // candidate were found
  // TODO: This is for now how it's done in athena... this look a bit
  // nasty to me
  if (fitterState.vtxInfoMap[&vtx].trackLinks.empty()) {
    // Find nearest track to vertex candidate
    double smallestDeltaZ = std::numeric_limits<double>::max();
    double newZ = 0;
    bool nearTrackFound = false;
    for (const auto& trk : seedTracks) {
      double zDistance = std::abs(m_extractParameters(*trk).position()[eZ] -
                                  vtx.position()[eZ]);
      if (zDistance < smallestDeltaZ) {
        smallestDeltaZ = zDistance;
        nearTrackFound = true;

        newZ = m_extractParameters(*trk).position()[eZ];
      }
    }
    if (nearTrackFound) {
      vtx.setFullPosition(SpacePointVector(0., 0., newZ, 0.));

      // Update vertex info for current vertex
      fitterState.vtxInfoMap[&vtx] =
          VertexInfo<InputTrack_t>(currentConstraint, vtx.fullPosition());

      // Try to add compatible track with adapted vertex position
      auto res = addCompatibleTracksToVertex(allTracks, vtx, fitterState);
      if (!res.ok()) {
        return Result<bool>::failure(res.error());
      }

      if (fitterState.vtxInfoMap[&vtx].trackLinks.empty()) {
        ACTS_DEBUG(
            "No tracks near seed were found, while at least one was "
            "expected. Break.");
        return Result<bool>::success(false);
      }

    } else {
      ACTS_DEBUG("No nearest track to seed found. Break.");
      return Result<bool>::success(false);
    }
  }

  return Result<bool>::success(true);
}

template <typename vfitter_t, typename sfinder_t>
auto Acts::AdaptiveMultiVertexFinder<vfitter_t, sfinder_t>::
    canPrepareVertexForFit(const std::vector<const InputTrack_t*>& allTracks,
                           const std::vector<const InputTrack_t*>& seedTracks,
                           Vertex<InputTrack_t>& vtx,
                           const Vertex<InputTrack_t>& currentConstraint,
                           FitterState_t& fitterState) const -> Result<bool> {
  // Add vertex info to fitter state
  fitterState.vtxInfoMap[&vtx] =
      VertexInfo<InputTrack_t>(currentConstraint, vtx.fullPosition());

  // Add all compatible tracks to vertex
  auto resComp = addCompatibleTracksToVertex(allTracks, vtx, fitterState);
  if (!resComp.ok()) {
    return Result<bool>::failure(resComp.error());
  }

  // Try to recover from cases where adding compatible track was not possible
  auto resRec = canRecoverFromNoCompatibleTracks(
      allTracks, seedTracks, vtx, currentConstraint, fitterState);
  if (!resRec.ok()) {
    return Result<bool>::failure(resRec.error());
  }

  return Result<bool>::success(*resRec);
}

template <typename vfitter_t, typename sfinder_t>
auto Acts::AdaptiveMultiVertexFinder<vfitter_t, sfinder_t>::
    checkVertexAndCompatibleTracks(
        Vertex<InputTrack_t>& vtx,
        const std::vector<const InputTrack_t*>& seedTracks,
        FitterState_t& fitterState) const -> std::pair<int, bool> {
  bool isGoodVertex = false;
  int nCompatibleTracks = 0;
  for (const auto& trk : fitterState.vtxInfoMap[&vtx].trackLinks) {
    const auto& trkAtVtx =
        fitterState.tracksAtVerticesMap.at(std::make_pair(trk, &vtx));
    if ((trkAtVtx.vertexCompatibility < m_cfg.maxVertexChi2 &&
         m_cfg.useFastCompatibility) ||
        (trkAtVtx.trackWeight > m_cfg.minWeight &&
         trkAtVtx.chi2Track < m_cfg.maxVertexChi2 &&
         !m_cfg.useFastCompatibility)) {
      // TODO: Understand why looking for compatible tracks only in seed tracks
      // and not also in all tracks
      auto foundIter =
          std::find_if(seedTracks.begin(), seedTracks.end(),
                       [&trk, this](auto seedTrk) { return trk == seedTrk; });
      if (foundIter != seedTracks.end()) {
        nCompatibleTracks++;
        ACTS_DEBUG("Compatible track found.");

        if (m_cfg.addSingleTrackVertices && m_cfg.useBeamSpotConstraint) {
          isGoodVertex = true;
          break;
        }
        if (nCompatibleTracks > 1) {
          isGoodVertex = true;
          break;
        }
      }
    }
  }  // end loop over all tracks at vertex

  return {nCompatibleTracks, isGoodVertex};
}

template <typename vfitter_t, typename sfinder_t>
auto Acts::AdaptiveMultiVertexFinder<vfitter_t, sfinder_t>::
    removeCompatibleTracksFromSeedTracks(
        Vertex<InputTrack_t>& vtx, std::vector<const InputTrack_t*>& seedTracks,
        FitterState_t& fitterState) const -> void {
  for (const auto& trk : fitterState.vtxInfoMap[&vtx].trackLinks) {
    const auto& trkAtVtx =
        fitterState.tracksAtVerticesMap.at(std::make_pair(trk, &vtx));
    if ((trkAtVtx.vertexCompatibility < m_cfg.maxVertexChi2 &&
         m_cfg.useFastCompatibility) ||
        (trkAtVtx.trackWeight > m_cfg.minWeight &&
         trkAtVtx.chi2Track < m_cfg.maxVertexChi2 &&
         !m_cfg.useFastCompatibility)) {
      // Find and remove track from seedTracks
      auto foundSeedIter =
          std::find_if(seedTracks.begin(), seedTracks.end(),
                       [&trk, this](auto seedTrk) { return trk == seedTrk; });
      if (foundSeedIter != seedTracks.end()) {
        seedTracks.erase(foundSeedIter);
      }
    }
  }
}

template <typename vfitter_t, typename sfinder_t>
auto Acts::AdaptiveMultiVertexFinder<vfitter_t, sfinder_t>::
    canRemoveNonCompatibleTrackFromSeedTracks(
        Vertex<InputTrack_t>& vtx, std::vector<const InputTrack_t*>& seedTracks,
        FitterState_t& fitterState) const -> bool {
  // Try to find the track with highest compatibility
  double maxCompatibility = 0;

  auto maxCompSeedIt = seedTracks.end();
  for (const auto& trk : fitterState.vtxInfoMap[&vtx].trackLinks) {
    const auto& trkAtVtx =
        fitterState.tracksAtVerticesMap.at(std::make_pair(trk, &vtx));
    double compatibility = trkAtVtx.vertexCompatibility;
    if (compatibility > maxCompatibility) {
      // Try to find track in seed tracks
      auto foundSeedIter =
          std::find_if(seedTracks.begin(), seedTracks.end(),
                       [&trk, this](auto seedTrk) { return trk == seedTrk; });
      if (foundSeedIter != seedTracks.end()) {
        maxCompatibility = compatibility;
        maxCompSeedIt = foundSeedIter;
      }
    }
  }
  if (maxCompSeedIt != seedTracks.end()) {
    // Remove track with highest compatibility from seed tracks
    seedTracks.erase(maxCompSeedIt);
  } else {
    // Could not find any seed with compatibility > 0, use alternative
    // method to remove a track from seed tracks: Closest track in z to
    // vtx candidate
    double smallestDeltaZ = std::numeric_limits<double>::max();
    auto smallestDzSeedIter = seedTracks.end();
    for (auto trkIter = seedTracks.begin(); trkIter != seedTracks.end();
         trkIter++) {
      double zDistance = std::abs(
          m_extractParameters(**trkIter).position()[eZ] - vtx.position()[eZ]);
      if (zDistance < smallestDeltaZ) {
        smallestDeltaZ = zDistance;
        smallestDzSeedIter = trkIter;
      }
    }
    if (smallestDzSeedIter != seedTracks.end()) {
      seedTracks.erase(smallestDzSeedIter);
    } else {
      ACTS_DEBUG("No track found to remove. Stop vertex finding now.");
      return false;
    }
  }
  return true;
}

template <typename vfitter_t, typename sfinder_t>
auto Acts::AdaptiveMultiVertexFinder<vfitter_t, sfinder_t>::keepNewVertex(
    Vertex<InputTrack_t>& vtx,
    const std::vector<Vertex<InputTrack_t>*>& allVertices,
    FitterState_t& fitterState) const -> bool {
  double contamination = 0.;
  double contaminationNum = 0;
  double contaminationDeNom = 0;
  for (const auto& trk : fitterState.vtxInfoMap[&vtx].trackLinks) {
    const auto& trkAtVtx =
        fitterState.tracksAtVerticesMap.at(std::make_pair(trk, &vtx));
    double trackWeight = trkAtVtx.trackWeight;
    contaminationNum += trackWeight * (1. - trackWeight);
    contaminationDeNom += trackWeight * trackWeight;
  }
  if (contaminationDeNom != 0) {
    contamination = contaminationNum / contaminationDeNom;
  }
  if (contamination > m_cfg.maximumVertexContamination) {
    return false;
  }

  if (isMergedVertex(vtx, allVertices)) {
    return false;
  }

  return true;
}

template <typename vfitter_t, typename sfinder_t>
auto Acts::AdaptiveMultiVertexFinder<vfitter_t, sfinder_t>::isMergedVertex(
    const Vertex<InputTrack_t>& vtx,
    const std::vector<Vertex<InputTrack_t>*>& allVertices) const -> bool {
  const SpacePointVector& candidatePos = vtx.fullPosition();
  const SpacePointSymMatrix& candidateCov = vtx.fullCovariance();
  const double candidateZPos = candidatePos[eZ];
  const double candidateZCov = candidateCov(eZ, eZ);

  for (const auto otherVtx : allVertices) {
    if (&vtx == otherVtx) {
      continue;
    }
    const SpacePointVector& otherPos = otherVtx->fullPosition();
    const SpacePointSymMatrix& otherCov = otherVtx->fullCovariance();
    const double otherZPos = otherPos[eZ];
    const double otherZCov = otherCov(eZ, eZ);

    const auto deltaPos = otherPos - candidatePos;
    const auto deltaZPos = otherZPos - candidateZPos;
    const auto sumCovZ = otherZCov + candidateZCov;

    double significance;
    if (not m_cfg.do3dSplitting) {
      // Use only z significance
      if (sumCovZ > 0.) {
        significance = std::abs(deltaZPos) / std::sqrt(sumCovZ);
      } else {
        return true;
      }
    } else {
      // Use full 3d information for significance
      auto sumCov = candidateCov + otherCov;
      significance =
          std::sqrt(deltaPos.dot((sumCov.inverse().eval()) * deltaPos));
    }
    if (significance < m_cfg.maxMergeVertexSignificance) {
      return true;
    }
  }
  return false;
}

template <typename vfitter_t, typename sfinder_t>
auto Acts::AdaptiveMultiVertexFinder<vfitter_t, sfinder_t>::deleteLastVertex(
    Vertex<InputTrack_t>& vtx,
    std::vector<std::unique_ptr<Vertex<InputTrack_t>>>& allVertices,
    std::vector<Vertex<InputTrack_t>*>& allVerticesPtr,
    FitterState_t& fitterState, FitterState_t& oldFitterState,
    const VertexFitterOptions<InputTrack_t>& vFitterOptions) const
    -> Result<void> {
  allVertices.pop_back();
  allVerticesPtr.pop_back();

  if (!m_cfg.refitAfterBadVertex) {
    fitterState.vertexCollection = oldFitterState.vertexCollection;
    fitterState.annealingState = oldFitterState.annealingState;
    fitterState.vtxInfoMap.clear();
    for (const auto& vtx : allVerticesPtr) {
      fitterState.vtxInfoMap.emplace(vtx, oldFitterState.vtxInfoMap[vtx]);
    }
    fitterState.trackToVerticesMultiMap =
        oldFitterState.trackToVerticesMultiMap;
    fitterState.tracksAtVerticesMap = oldFitterState.tracksAtVerticesMap;

  } else {
    // Update fitter state with removed vertex candidate
    fitterState.removeVertexFromMultiMap(vtx);

    // TODO: clean tracksAtVerticesMap maybe here? i.e. remove all entries
    // with old vertex?

    // Do the fit with removed vertex
    auto fitResult = m_cfg.vertexFitter.fit(fitterState, allVerticesPtr,
                                            m_cfg.linearizer, vFitterOptions);
    if (!fitResult.ok()) {
      return fitResult.error();
    }
  }
  return {};
}

template <typename vfitter_t, typename sfinder_t>
auto Acts::AdaptiveMultiVertexFinder<vfitter_t, sfinder_t>::getVertexOutputList(
    const std::vector<Vertex<InputTrack_t>*>& allVerticesPtr,
    FitterState_t& fitterState) const -> std::vector<Vertex<InputTrack_t>> {
  std::vector<Vertex<InputTrack_t>> outputVec;
  for (auto vtx : allVerticesPtr) {
    auto& outVtx = *vtx;
    std::vector<TrackAtVertex<InputTrack_t>> tracksAtVtx;
    for (const auto& trk : fitterState.vtxInfoMap[vtx].trackLinks) {
      tracksAtVtx.push_back(
          fitterState.tracksAtVerticesMap.at(std::make_pair(trk, vtx)));
    }
    outVtx.setTracksAtVertex(tracksAtVtx);
    outputVec.push_back(outVtx);
  }
  return outputVec;
}
