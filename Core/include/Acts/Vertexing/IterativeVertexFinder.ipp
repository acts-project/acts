// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

template <typename vfitter_t, typename sfinder_t>
auto Acts::IterativeVertexFinder<vfitter_t, sfinder_t>::find(
    const std::vector<const InputTrack_t*>& trackVector,
    const VertexingOptions<InputTrack_t>& vertexingOptions, State& state) const
    -> Result<std::vector<Vertex<InputTrack_t>>> {
  // Original tracks
  const std::vector<const InputTrack_t*>& origTracks = trackVector;
  // Tracks for seeding
  std::vector<const InputTrack_t*> seedTracks = trackVector;

  // List of vertices to be filled below
  std::vector<Vertex<InputTrack_t>> vertexCollection;

  int nInterations = 0;
  // begin iterating
  while (seedTracks.size() > 1 && nInterations < m_cfg.maxVertices) {
    /// Do seeding
    auto seedRes = getVertexSeed(seedTracks, vertexingOptions);

    if (!seedRes.ok()) {
      return seedRes.error();
    }

    const auto& seedVertex = *seedRes;

    if (seedVertex.fullPosition()[eZ] ==
        vertexingOptions.constraint.position().z()) {
      ACTS_DEBUG("No more seed found. Break and stop primary vertex finding.");
      break;
    }

    /// End seeding
    /// Now take only tracks compatible with current seed
    // Tracks used for the fit in this iteration
    std::vector<const InputTrack_t*> tracksToFit;
    std::vector<const InputTrack_t*> tracksToFitSplitVertex;

    // Fill vector with tracks to fit, only compatible with seed:
    auto res = fillTracksToFit(seedTracks, seedVertex, tracksToFit,
                               tracksToFitSplitVertex, vertexingOptions, state);

    if (!res.ok()) {
      return res.error();
    }

    ACTS_DEBUG("Number of tracks used for fit: " << tracksToFit.size());

    /// Begin vertex fit
    Vertex<InputTrack_t> currentVertex;
    Vertex<InputTrack_t> currentSplitVertex;

    if (vertexingOptions.useConstraintInFit && !tracksToFit.empty()) {
      auto fitResult = m_cfg.vertexFitter.fit(
          tracksToFit, m_cfg.linearizer, vertexingOptions, state.fitterState);
      if (fitResult.ok()) {
        currentVertex = std::move(*fitResult);
      } else {
        return fitResult.error();
      }
    } else if (!vertexingOptions.useConstraintInFit && tracksToFit.size() > 1) {
      auto fitResult = m_cfg.vertexFitter.fit(
          tracksToFit, m_cfg.linearizer, vertexingOptions, state.fitterState);
      if (fitResult.ok()) {
        currentVertex = std::move(*fitResult);
      } else {
        return fitResult.error();
      }
    }
    if (m_cfg.createSplitVertices && tracksToFitSplitVertex.size() > 1) {
      auto fitResult =
          m_cfg.vertexFitter.fit(tracksToFitSplitVertex, m_cfg.linearizer,
                                 vertexingOptions, state.fitterState);
      if (fitResult.ok()) {
        currentSplitVertex = std::move(*fitResult);
      } else {
        return fitResult.error();
      }
    }
    /// End vertex fit
    ACTS_DEBUG("Vertex position after fit: "
               << currentVertex.fullPosition().transpose());

    // Number degrees of freedom
    double ndf = currentVertex.fitQuality().second;
    double ndfSplitVertex = currentSplitVertex.fitQuality().second;

    // Number of significant tracks
    int nTracksAtVertex = countSignificantTracks(currentVertex);
    int nTracksAtSplitVertex = countSignificantTracks(currentSplitVertex);

    bool isGoodVertex = ((!vertexingOptions.useConstraintInFit && ndf > 0 &&
                          nTracksAtVertex >= 2) ||
                         (vertexingOptions.useConstraintInFit && ndf > 3 &&
                          nTracksAtVertex >= 2));

    if (!isGoodVertex) {
      removeTracks(tracksToFit, seedTracks);
    } else {
      if (m_cfg.reassignTracksAfterFirstFit && (!m_cfg.createSplitVertices)) {
        // vertex is good vertex here
        // but add tracks which may have been missed

        auto result = reassignTracksToNewVertex(
            vertexCollection, currentVertex, tracksToFit, seedTracks,
            origTracks, vertexingOptions, state);
        if (!result.ok()) {
          return result.error();
        }
        isGoodVertex = *result;

      }  // end reassignTracksAfterFirstFit case
         // still good vertex? might have changed in the meanwhile
      if (isGoodVertex) {
        removeUsedCompatibleTracks(currentVertex, tracksToFit, seedTracks,
                                   vertexingOptions, state);

        ACTS_DEBUG(
            "Number of seed tracks after removal of compatible tracks "
            "and outliers: "
            << seedTracks.size());
      }
    }  // end case if good vertex

    // now splitvertex
    bool isGoodSplitVertex = false;
    if (m_cfg.createSplitVertices) {
      isGoodSplitVertex = (ndfSplitVertex > 0 && nTracksAtSplitVertex >= 2);

      if (!isGoodSplitVertex) {
        removeTracks(tracksToFitSplitVertex, seedTracks);
      } else {
        removeUsedCompatibleTracks(currentSplitVertex, tracksToFitSplitVertex,
                                   seedTracks, vertexingOptions, state);
      }
    }
    // Now fill vertex collection with vertex
    if (isGoodVertex) {
      vertexCollection.push_back(currentVertex);
    }
    if (isGoodSplitVertex && m_cfg.createSplitVertices) {
      vertexCollection.push_back(currentSplitVertex);
    }

    nInterations++;
  }  // end while loop

  return vertexCollection;
}

template <typename vfitter_t, typename sfinder_t>
auto Acts::IterativeVertexFinder<vfitter_t, sfinder_t>::getVertexSeed(
    const std::vector<const InputTrack_t*>& seedTracks,
    const VertexingOptions<InputTrack_t>& vertexingOptions) const
    -> Result<Vertex<InputTrack_t>> {
  typename sfinder_t::State finderState;
  auto res = m_cfg.seedFinder.find(seedTracks, vertexingOptions, finderState);

  if (!res.ok()) {
    ACTS_DEBUG("Seeding error: internal. Number of input tracks: "
               << seedTracks.size());
    return VertexingError::SeedingError;
  }

  const auto& vertexCollection = *res;

  if (vertexCollection.empty()) {
    ACTS_DEBUG("Seeding error: no seeds. Number of input tracks: "
               << seedTracks.size());
    return VertexingError::SeedingError;
  }

  ACTS_DEBUG("Found " << vertexCollection.size() << " seeds");

  // retrieve the seed vertex as the last element in
  // the seed vertexCollection
  Vertex<InputTrack_t> seedVertex = vertexCollection.back();

  ACTS_DEBUG("Considering seed at position: ("
             << seedVertex.fullPosition()[eX] << ", "
             << seedVertex.fullPosition()[eY] << ", "
             << seedVertex.fullPosition()[eZ] << ", " << seedVertex.time()
             << "). Number of input tracks: " << seedTracks.size());

  return seedVertex;
}

template <typename vfitter_t, typename sfinder_t>
void Acts::IterativeVertexFinder<vfitter_t, sfinder_t>::removeTracks(
    const std::vector<const InputTrack_t*>& tracksToRemove,
    std::vector<const InputTrack_t*>& seedTracks) const {
  for (const auto& trk : tracksToRemove) {
    const BoundTrackParameters& params = m_extractParameters(*trk);
    // Find track in seedTracks
    auto foundIter =
        std::find_if(seedTracks.begin(), seedTracks.end(),
                     [&params, this](const auto seedTrk) {
                       return params == m_extractParameters(*seedTrk);
                     });
    if (foundIter != seedTracks.end()) {
      // Remove track from seed tracks
      seedTracks.erase(foundIter);
    } else {
      ACTS_WARNING("Track to be removed not found in seed tracks.")
    }
  }
}

template <typename vfitter_t, typename sfinder_t>
Acts::Result<double>
Acts::IterativeVertexFinder<vfitter_t, sfinder_t>::getCompatibility(
    const BoundTrackParameters& params, const Vertex<InputTrack_t>& vertex,
    const Surface& perigeeSurface,
    const VertexingOptions<InputTrack_t>& vertexingOptions,
    State& state) const {
  // Linearize track
  auto result = m_cfg.linearizer.linearizeTrack(
      params, vertex.fullPosition()[3], perigeeSurface,
      vertexingOptions.geoContext, vertexingOptions.magFieldContext,
      state.linearizerState);
  if (!result.ok()) {
    return result.error();
  }

  auto linTrack = std::move(*result);

  // Calculate reduced weight
  SquareMatrix2 weightReduced =
      linTrack.covarianceAtPCA.template block<2, 2>(0, 0);

  SquareMatrix2 errorVertexReduced =
      (linTrack.positionJacobian *
       (vertex.fullCovariance() * linTrack.positionJacobian.transpose()))
          .template block<2, 2>(0, 0);
  weightReduced += errorVertexReduced;
  weightReduced = weightReduced.inverse().eval();

  // Calculate compatibility / chi2
  Vector2 trackParameters2D =
      linTrack.parametersAtPCA.template block<2, 1>(0, 0);
  double compatibility =
      trackParameters2D.dot(weightReduced * trackParameters2D);

  return compatibility;
}

template <typename vfitter_t, typename sfinder_t>
Acts::Result<void>
Acts::IterativeVertexFinder<vfitter_t, sfinder_t>::removeUsedCompatibleTracks(
    Vertex<InputTrack_t>& vertex, std::vector<const InputTrack_t*>& tracksToFit,
    std::vector<const InputTrack_t*>& seedTracks,
    const VertexingOptions<InputTrack_t>& vertexingOptions,
    State& state) const {
  std::vector<TrackAtVertex<InputTrack_t>> tracksAtVertex = vertex.tracks();

  for (const auto& trackAtVtx : tracksAtVertex) {
    // Check compatibility
    if (trackAtVtx.trackWeight < m_cfg.cutOffTrackWeight) {
      // Do not remove track here, since it is not compatible with the vertex
      continue;
    }
    // Find and remove track from seedTracks
    auto foundSeedIter =
        std::find_if(seedTracks.begin(), seedTracks.end(),
                     [&trackAtVtx](const auto seedTrk) {
                       return trackAtVtx.originalParams == seedTrk;
                     });
    if (foundSeedIter != seedTracks.end()) {
      seedTracks.erase(foundSeedIter);
    } else {
      ACTS_WARNING("Track trackAtVtx not found in seedTracks!");
    }

    // Find and remove track from tracksToFit
    auto foundFitIter =
        std::find_if(tracksToFit.begin(), tracksToFit.end(),
                     [&trackAtVtx](const auto fitTrk) {
                       return trackAtVtx.originalParams == fitTrk;
                     });
    if (foundFitIter != tracksToFit.end()) {
      tracksToFit.erase(foundFitIter);
    } else {
      ACTS_WARNING("Track trackAtVtx not found in tracksToFit!");
    }
  }  // end iteration over tracksAtVertex

  ACTS_DEBUG("After removal of tracks belonging to vertex, "
             << seedTracks.size() << " seed tracks left.");

  // Now start considering outliers
  // tracksToFit that are left here were below
  // m_cfg.cutOffTrackWeight threshold and are hence outliers
  ACTS_DEBUG("Number of outliers: " << tracksToFit.size());

  const std::shared_ptr<PerigeeSurface> vertexPerigeeSurface =
      Surface::makeShared<PerigeeSurface>(
          VectorHelpers::position(vertex.fullPosition()));

  for (const auto& trk : tracksToFit) {
    // calculate chi2 w.r.t. last fitted vertex
    auto result =
        getCompatibility(m_extractParameters(*trk), vertex,
                         *vertexPerigeeSurface, vertexingOptions, state);

    if (!result.ok()) {
      return result.error();
    }

    double chi2 = *result;

    // check if sufficiently compatible with last fitted vertex
    // (quite loose constraint)
    if (chi2 < m_cfg.maximumChi2cutForSeeding) {
      auto foundIter =
          std::find_if(seedTracks.begin(), seedTracks.end(),
                       [&trk](const auto seedTrk) { return trk == seedTrk; });
      if (foundIter != seedTracks.end()) {
        // Remove track from seed tracks
        seedTracks.erase(foundIter);
      }

    } else {
      // Track not compatible with vertex
      // Remove track from current vertex
      auto foundIter = std::find_if(
          tracksAtVertex.begin(), tracksAtVertex.end(),
          [&trk](auto trkAtVtx) { return trk == trkAtVtx.originalParams; });
      if (foundIter != tracksAtVertex.end()) {
        // Remove track from seed tracks
        tracksAtVertex.erase(foundIter);
      }
    }
  }

  // set updated (possibly with removed outliers) tracksAtVertex to vertex
  vertex.setTracksAtVertex(tracksAtVertex);

  return {};
}

template <typename vfitter_t, typename sfinder_t>
Acts::Result<void>
Acts::IterativeVertexFinder<vfitter_t, sfinder_t>::fillTracksToFit(
    const std::vector<const InputTrack_t*>& seedTracks,
    const Vertex<InputTrack_t>& seedVertex,
    std::vector<const InputTrack_t*>& tracksToFitOut,
    std::vector<const InputTrack_t*>& tracksToFitSplitVertexOut,
    const VertexingOptions<InputTrack_t>& vertexingOptions,
    State& state) const {
  int numberOfTracks = seedTracks.size();

  // Count how many tracks are used for fit
  int count = 0;
  // Fill tracksToFit vector with tracks compatible with seed
  for (const auto& sTrack : seedTracks) {
    // If there are only few tracks left, add them to fit regardless of their
    // position:
    if (numberOfTracks <= 2) {
      tracksToFitOut.push_back(sTrack);
      ++count;
    } else if (numberOfTracks <= 4 && !m_cfg.createSplitVertices) {
      tracksToFitOut.push_back(sTrack);
      ++count;
    } else if (numberOfTracks <= 4 * m_cfg.splitVerticesTrkInvFraction &&
               m_cfg.createSplitVertices) {
      if (count % m_cfg.splitVerticesTrkInvFraction != 0) {
        tracksToFitOut.push_back(sTrack);
        ++count;
      } else {
        tracksToFitSplitVertexOut.push_back(sTrack);
        ++count;
      }
    }
    // If a large amount of tracks is available, we check their compatibility
    // with the vertex before adding them to the fit:
    else {
      const BoundTrackParameters& sTrackParams = m_extractParameters(*sTrack);
      auto distanceRes = m_cfg.ipEst.calculateDistance(
          vertexingOptions.geoContext, sTrackParams, seedVertex.position(),
          state.ipState);
      if (!distanceRes.ok()) {
        return distanceRes.error();
      }

      if (!sTrackParams.covariance()) {
        return VertexingError::NoCovariance;
      }

      // sqrt(sigma(d0)^2+sigma(z0)^2), where sigma(d0)^2 is the variance of d0
      double hypotVariance =
          sqrt((*(sTrackParams.covariance()))(eBoundLoc0, eBoundLoc0) +
               (*(sTrackParams.covariance()))(eBoundLoc1, eBoundLoc1));

      if (hypotVariance == 0.) {
        ACTS_WARNING(
            "Track impact parameter covariances are zero. Track was not "
            "assigned to vertex.");
        continue;
      }

      if (*distanceRes / hypotVariance < m_cfg.significanceCutSeeding) {
        if (!m_cfg.createSplitVertices ||
            count % m_cfg.splitVerticesTrkInvFraction != 0) {
          tracksToFitOut.push_back(sTrack);
          ++count;
        } else {
          tracksToFitSplitVertexOut.push_back(sTrack);
          ++count;
        }
      }
    }
  }
  return {};
}

template <typename vfitter_t, typename sfinder_t>
Acts::Result<bool>
Acts::IterativeVertexFinder<vfitter_t, sfinder_t>::reassignTracksToNewVertex(
    std::vector<Vertex<InputTrack_t>>& vertexCollection,
    Vertex<InputTrack_t>& currentVertex,
    std::vector<const InputTrack_t*>& tracksToFit,
    std::vector<const InputTrack_t*>& seedTracks,
    const std::vector<const InputTrack_t*>& /* origTracks */,
    const VertexingOptions<InputTrack_t>& vertexingOptions,
    State& state) const {
  int numberOfAddedTracks = 0;

  const std::shared_ptr<PerigeeSurface> currentVertexPerigeeSurface =
      Surface::makeShared<PerigeeSurface>(
          VectorHelpers::position(currentVertex.fullPosition()));

  // iterate over all vertices and check if tracks need to be reassigned
  // to new (current) vertex
  for (auto& vertexIt : vertexCollection) {
    // tracks at vertexIt
    std::vector<TrackAtVertex<InputTrack_t>> tracksAtVertex = vertexIt.tracks();
    auto tracksBegin = tracksAtVertex.begin();
    auto tracksEnd = tracksAtVertex.end();

    const std::shared_ptr<PerigeeSurface> vertexItPerigeeSurface =
        Surface::makeShared<PerigeeSurface>(
            VectorHelpers::position(vertexIt.fullPosition()));

    for (auto tracksIter = tracksBegin; tracksIter != tracksEnd;) {
      // consider only tracks that are not too tightly assigned to other
      // vertex
      if (tracksIter->trackWeight > m_cfg.cutOffTrackWeightReassign) {
        tracksIter++;
        continue;
      }
      // use original perigee parameters
      BoundTrackParameters origParams =
          m_extractParameters(*(tracksIter->originalParams));

      // compute compatibility
      auto resultNew = getCompatibility(origParams, currentVertex,
                                        *currentVertexPerigeeSurface,
                                        vertexingOptions, state);
      if (!resultNew.ok()) {
        return Result<bool>::failure(resultNew.error());
      }
      double chi2NewVtx = *resultNew;

      auto resultOld =
          getCompatibility(origParams, vertexIt, *vertexItPerigeeSurface,
                           vertexingOptions, state);
      if (!resultOld.ok()) {
        return Result<bool>::failure(resultOld.error());
      }
      double chi2OldVtx = *resultOld;

      ACTS_DEBUG("Compatibility to new vs old vertex: " << chi2NewVtx << " vs "
                                                        << chi2OldVtx);

      if (chi2NewVtx < chi2OldVtx) {
        tracksToFit.push_back(tracksIter->originalParams);
        // origTrack was already deleted from seedTracks previously
        // (when assigned to old vertex)
        // add it now back to seedTracks to be able to consistently
        // delete it later
        // when all tracks used to fit current vertex are deleted
        seedTracks.push_back(tracksIter->originalParams);
        // seedTracks.push_back(*std::find_if(
        //     origTracks.begin(), origTracks.end(),
        //     [&origParams, this](auto origTrack) {
        //       return origParams == m_extractParameters(*origTrack);
        //     }));

        numberOfAddedTracks += 1;

        // remove track from old vertex
        tracksIter = tracksAtVertex.erase(tracksIter);
        tracksBegin = tracksAtVertex.begin();
        tracksEnd = tracksAtVertex.end();

      }  // end chi2NewVtx < chi2OldVtx

      else {
        // go and check next track
        ++tracksIter;
      }
    }  // end loop over tracks at old vertexIt

    vertexIt.setTracksAtVertex(tracksAtVertex);
  }  // end loop over all vertices

  ACTS_DEBUG("Added " << numberOfAddedTracks
                      << " tracks from old (other) vertices for new fit");

  // override current vertex with new fit
  // set first to default vertex to be able to check if still good vertex
  // later
  currentVertex = Vertex<InputTrack_t>();
  if (vertexingOptions.useConstraintInFit && !tracksToFit.empty()) {
    auto fitResult = m_cfg.vertexFitter.fit(
        tracksToFit, m_cfg.linearizer, vertexingOptions, state.fitterState);
    if (fitResult.ok()) {
      currentVertex = std::move(*fitResult);
    } else {
      return Result<bool>::success(false);
    }
  } else if (!vertexingOptions.useConstraintInFit && tracksToFit.size() > 1) {
    auto fitResult = m_cfg.vertexFitter.fit(
        tracksToFit, m_cfg.linearizer, vertexingOptions, state.fitterState);
    if (fitResult.ok()) {
      currentVertex = std::move(*fitResult);
    } else {
      return Result<bool>::success(false);
    }
  }

  // Number degrees of freedom
  double ndf = currentVertex.fitQuality().second;

  // Number of significant tracks
  int nTracksAtVertex = countSignificantTracks(currentVertex);

  bool isGoodVertex = ((!vertexingOptions.useConstraintInFit && ndf > 0 &&
                        nTracksAtVertex >= 2) ||
                       (vertexingOptions.useConstraintInFit && ndf > 3 &&
                        nTracksAtVertex >= 2));

  if (!isGoodVertex) {
    removeTracks(tracksToFit, seedTracks);

    ACTS_DEBUG("Going to new iteration with "
               << seedTracks.size() << "seed tracks after BAD vertex.");
  }

  return Result<bool>::success(isGoodVertex);
}

template <typename vfitter_t, typename sfinder_t>
int Acts::IterativeVertexFinder<vfitter_t, sfinder_t>::countSignificantTracks(
    const Vertex<InputTrack_t>& vtx) const {
  return std::count_if(vtx.tracks().begin(), vtx.tracks().end(),
                       [this](TrackAtVertex<InputTrack_t> trk) {
                         return trk.trackWeight > m_cfg.cutOffTrackWeight;
                       });
}
