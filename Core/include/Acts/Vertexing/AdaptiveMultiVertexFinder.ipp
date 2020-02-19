// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Vertexing/VertexFitterOptions.hpp"
#include "Acts/Vertexing/VertexingError.hpp"

template <typename vfitter_t, typename sfinder_t>
auto Acts::AdaptiveMultiVertexFinder<vfitter_t, sfinder_t>::find(
    const std::vector<InputTrack_t>& allTracks,
    const VertexFinderOptions<InputTrack_t>& vFinderOptions) const
    -> Result<std::vector<Vertex<InputTrack_t>>> {
  if (allTracks.empty()) {
    return VertexingError::EmptyInput;
  }

  // Create copy of finder options, will be modified after seeding
  // to set the correct vertex constraint
  VertexFinderOptions<InputTrack_t> finderOptions = vFinderOptions;

  // Original tracks
  const std::vector<InputTrack_t>& origTracks = allTracks;
  // Tracks for seeding
  // Note: Remains to be investigated if another container (e.g. std::list)
  // or also std::vector<InputTrack_t*> is a faster option since erasures
  // of tracks is quite expensive with std::vector.
  // std::vector<InputTrack_t*> would however also come with an overhead
  // since m_cfg.vertexFitter.fit and m_cfg.seedFinder.find take
  // vector<InputTrack_t> and hence a lot of copying would be required.
  // Maybe use std::vector<InputTrack_t*> and adapt fit accordingly to
  // also take pointers to tracks instead of the track object.
  std::vector<InputTrack_t> seedTracks = allTracks;

  // Construct the vertex fitter options from vertex finder options
  VertexFitterOptions<InputTrack_t> vFitterOptions(
      finderOptions.geoContext, finderOptions.magFieldContext,
      finderOptions.vertexConstraint);

  FitterState_t fitterState;

  std::vector<Vertex<InputTrack_t>> allVertices;

  int iteration = 0;
  while (((m_cfg.addSingleTrackVertices && seedTracks.size() > 0) ||
          ((!m_cfg.addSingleTrackVertices) && seedTracks.size() > 1)) &&
         iteration < m_cfg.maxIterations) {
    // Tracks that are used for searching compatible tracks
    // near a vertex candidate
    // TODO: This involves a lot of copying. Change the way of accessing tracks
    std::vector<InputTrack_t> myTracks;
    if (m_cfg.realMultiVertex == true) {
      myTracks = origTracks;
    } else {
      myTracks = seedTracks;
    }

    // Retrieve seed vertex from all remaining seedTracks
    auto seedRes = doSeeding(seedTracks, finderOptions);
    if (!seedRes.ok()) {
      return seedRes.error();
    }
    Vertex<InputTrack_t> vtxCandidate = *seedRes;

    if (vtxCandidate.position().z() == 0.) {
      ACTS_DEBUG(
          "No seed found anymore. Break and stop primary vertex finding.");
      break;
    }

    auto resPrep = canPrepareVertexForFit(myTracks, seedTracks, vtxCandidate,
                                          finderOptions, fitterState);
    if (!resPrep.ok()) {
      return resPrep.error();
    }
    if (!(*resPrep)) {
      // Could not prepare the vertex for the fit anymore, break.
      break;
    }

    // state.vtxInfoMap[&vtxCandidate] = vtxInfo1;
  }

  return allVertices;
}

template <typename vfitter_t, typename sfinder_t>
auto Acts::AdaptiveMultiVertexFinder<vfitter_t, sfinder_t>::doSeeding(
    const std::vector<InputTrack_t>& trackVector,
    VertexFinderOptions<InputTrack_t>& vFinderOptions) const
    -> Result<Vertex<InputTrack_t>> {
  // Run seed finder
  auto seedRes = m_cfg.seedFinder.find(trackVector, vFinderOptions);

  if (!seedRes.ok()) {
    return seedRes.error();
  }

  Vertex<InputTrack_t> seedVertex = (*seedRes).back();

  if (m_cfg.useBeamSpotConstraint) {
    if (m_cfg.useSeedConstraint) {
      vFinderOptions.vertexConstraint.setFullPosition(
          seedVertex.fullPosition());
      vFinderOptions.vertexConstraint.setFullCovariance(
          seedVertex.fullCovariance());
    }
  } else {
    vFinderOptions.vertexConstraint.setFullPosition(seedVertex.fullPosition());
    vFinderOptions.vertexConstraint.setFullCovariance(
        SpacePointSymMatrix::Identity() * m_cfg.looseConstrValue);
    vFinderOptions.vertexConstraint.setFitQuality(
        m_cfg.defaultConstrFitQuality);
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
    const BoundParameters& track, const Vertex<InputTrack_t>& vtx) const
    -> Result<double> {
  // TODO: In original implementation the covariance of the given vertex is set
  // to zero. I did the same here now, but consider removing this and just
  // passing the vtx object to the estimator without changing its covariance.
  // After all, the vertex seed does have a non-zero convariance in general and
  // it probably should be used.
  Vertex<InputTrack_t> newVtx = vtx;
  newVtx.setFullCovariance(SpacePointSymMatrix::Zero());

  double significance = 0.;

  auto estRes = m_cfg.ipEstimator.estimate(track, newVtx);
  if (!estRes.ok()) {
    return estRes.error();
  }

  auto ipas = std::move(*estRes);

  if (ipas->sigmad0 > 0 && ipas->sigmaz0 > 0) {
    significance = std::sqrt(std::pow(ipas->IPd0 / ipas->sigmad0, 2) +
                             std::pow(ipas->IPz0 / ipas->sigmaz0, 2));
  }

  return significance;
}

template <typename vfitter_t, typename sfinder_t>
auto Acts::AdaptiveMultiVertexFinder<vfitter_t, sfinder_t>::
    addCompatibleTracksToVertex(const std::vector<InputTrack_t>& tracks,
                                Vertex<InputTrack_t>& vtx) const
    -> Result<void> {
  std::vector<TrackAtVertex<InputTrack_t>> tracksAtVtx;

  for (const auto& trk : tracks) {
    BoundParameters params = m_extractParameters(trk);
    auto sigRes = getIPSignificance(params, vtx);
    if (!sigRes.ok()) {
      return sigRes.error();
    }
    double ipSig = *sigRes;
    if ((std::abs(estimateDeltaZ(params, vtx.position())) <
         m_cfg.tracksMaxZinterval) &&
        (ipSig < m_cfg.tracksMaxSignificance)) {
      tracksAtVtx.push_back(TrackAtVertex(params, trk));
    }
  }

  vtx.setTracksAtVertex(tracksAtVtx);

  return {};
}

template <typename vfitter_t, typename sfinder_t>
auto Acts::AdaptiveMultiVertexFinder<vfitter_t, sfinder_t>::
    canRecoverFromNoCompatibleTracks(
        const std::vector<InputTrack_t>& myTracks,
        const std::vector<InputTrack_t>& seedTracks, Vertex<InputTrack_t>& vtx,
        const VertexFinderOptions<InputTrack_t>& vFinderOptions,
        FitterState_t& fitterState) const -> Result<bool> {
  // Recover from cases where no compatible tracks to vertex
  // candidate were found
  // TODO: This is for now how it's done in athena... this look a bit
  // nasty to me
  if (vtx.tracks().empty()) {
    double smallestDeltaZ = std::numeric_limits<double>::max();
    double newZ = 0;
    bool nearTrackFound = false;
    for (const auto& trk : seedTracks) {
      double zDistance = std::abs(m_extractParameters(trk).position()[eZ] -
                                  vtx.position()[eZ]);
      if (zDistance < smallestDeltaZ) {
        smallestDeltaZ = zDistance;
        nearTrackFound = true;
        newZ = m_extractParameters(trk).position()[eZ];
      }
    }
    if (nearTrackFound) {
      // TODO: check athena actualcandidate position here (has not changed?)
      // TODO: so do I want to change the vtx position here?
      vtx.setFullPosition(SpacePointVector(0., 0., newZ, 0.));

      // TODO: do sth with fitterstate here
      // vtxInfo = VertexInfo<InputTrack_t>(vFinderOptions.vertexConstraint,
      //                                   vtx.fullPosition());

      // Try to add compatible track with adapted vertex position
      auto res = addCompatibleTracksToVertex(myTracks, vtx);
      if (!res.ok()) {
        return Result<bool>::failure(res.error());
      }

      if (vtx.tracks().empty()) {
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
    canPrepareVertexForFit(
        const std::vector<InputTrack_t>& myTracks,
        const std::vector<InputTrack_t>& seedTracks, Vertex<InputTrack_t>& vtx,
        const VertexFinderOptions<InputTrack_t>& vFinderOptions,
        FitterState_t& fitterState) const -> Result<bool> {
  // TODO: do something with state here
  VertexInfo<InputTrack_t> vtxCandidateInfo(vFinderOptions.vertexConstraint,
                                            vtx.fullPosition());

  auto resComp = addCompatibleTracksToVertex(myTracks, vtx);
  if (!resComp.ok()) {
    return Result<bool>::failure(resComp.error());
  }

  auto resRec = canRecoverFromNoCompatibleTracks(myTracks, seedTracks, vtx,
                                                 vFinderOptions, fitterState);
  if (!resRec.ok()) {
    return Result<bool>::failure(resRec.error());
  }

  return Result<bool>::success(*resRec);
}
