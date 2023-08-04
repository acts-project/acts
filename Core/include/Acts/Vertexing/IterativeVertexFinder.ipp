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
        vertexingOptions.vertexConstraint.position().z()) {
      ACTS_DEBUG("No more seed found. Break and stop primary vertex finding.");
      break;
    }

    /// End seeding
    /// Now take only tracks compatible with current seed
    // Tracks used for the fit in this iteration
    std::vector<const InputTrack_t*> perigeesToFit;
    std::vector<const InputTrack_t*> perigeesToFitSplitVertex;

    // Fill vector with tracks to fit, only compatible with seed:
    auto res =
        fillPerigeesToFit(seedTracks, seedVertex, perigeesToFit,
                          perigeesToFitSplitVertex, vertexingOptions, state);

    if (!res.ok()) {
      return res.error();
    }

    ACTS_DEBUG("Perigees used for fit: " << perigeesToFit.size());

    /// Begin vertex fit
    Vertex<InputTrack_t> currentVertex;
    Vertex<InputTrack_t> currentSplitVertex;

    if (m_cfg.useBeamConstraint && !perigeesToFit.empty()) {
      auto fitResult = m_cfg.vertexFitter.fit(
          perigeesToFit, m_cfg.linearizer, vertexingOptions, state.fitterState);
      if (fitResult.ok()) {
        currentVertex = std::move(*fitResult);
      } else {
        return fitResult.error();
      }
    } else if (!m_cfg.useBeamConstraint && perigeesToFit.size() > 1) {
      auto fitResult = m_cfg.vertexFitter.fit(
          perigeesToFit, m_cfg.linearizer, vertexingOptions, state.fitterState);
      if (fitResult.ok()) {
        currentVertex = std::move(*fitResult);
      } else {
        return fitResult.error();
      }
    }
    if (m_cfg.createSplitVertices && perigeesToFitSplitVertex.size() > 1) {
      auto fitResult =
          m_cfg.vertexFitter.fit(perigeesToFitSplitVertex, m_cfg.linearizer,
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

    bool isGoodVertex =
        ((!m_cfg.useBeamConstraint && ndf > 0 && nTracksAtVertex >= 2) ||
         (m_cfg.useBeamConstraint && ndf > 3 && nTracksAtVertex >= 2));

    if (!isGoodVertex) {
      removeAllTracks(perigeesToFit, seedTracks);
    } else {
      if (m_cfg.reassignTracksAfterFirstFit && (!m_cfg.createSplitVertices)) {
        // vertex is good vertex here
        // but add tracks which may have been missed

        auto result = reassignTracksToNewVertex(
            vertexCollection, currentVertex, perigeesToFit, seedTracks,
            origTracks, vertexingOptions, state);
        if (!result.ok()) {
          return result.error();
        }
        isGoodVertex = *result;

      }  // end reassignTracksAfterFirstFit case
         // still good vertex? might have changed in the meanwhile
      if (isGoodVertex) {
        removeUsedCompatibleTracks(currentVertex, perigeesToFit, seedTracks,
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
        removeAllTracks(perigeesToFitSplitVertex, seedTracks);
      } else {
        removeUsedCompatibleTracks(currentSplitVertex, perigeesToFitSplitVertex,
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
void Acts::IterativeVertexFinder<vfitter_t, sfinder_t>::removeAllTracks(
    const std::vector<const InputTrack_t*>& perigeesToFit,
    std::vector<const InputTrack_t*>& seedTracks) const {
  for (const auto& fitPerigee : perigeesToFit) {
    const BoundTrackParameters& fitPerigeeParams =
        m_extractParameters(*fitPerigee);
    // Find track in seedTracks
    auto foundIter =
        std::find_if(seedTracks.begin(), seedTracks.end(),
                     [&fitPerigeeParams, this](const auto seedTrk) {
                       return fitPerigeeParams == m_extractParameters(*seedTrk);
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
    const VertexingOptions<InputTrack_t>& vertexingOptions,
    State& state) const {
  // Linearize track
  auto result = m_cfg.linearizer.linearizeTrack(
      params, vertex.fullPosition(), vertexingOptions.geoContext,
      vertexingOptions.magFieldContext, state.linearizerState);
  if (!result.ok()) {
    return result.error();
  }

  auto linTrack = std::move(*result);

  // Calculate reduced weight
  SymMatrix2 weightReduced =
      linTrack.covarianceAtPCA.template block<2, 2>(0, 0);

  SymMatrix2 errorVertexReduced =
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
    Vertex<InputTrack_t>& myVertex,
    std::vector<const InputTrack_t*>& perigeesToFit,
    std::vector<const InputTrack_t*>& seedTracks,
    const VertexingOptions<InputTrack_t>& vertexingOptions,
    State& state) const {
  std::vector<TrackAtVertex<InputTrack_t>> tracksAtVertex = myVertex.tracks();

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

    // Find and remove track from perigeesToFit
    auto foundFitIter =
        std::find_if(perigeesToFit.begin(), perigeesToFit.end(),
                     [&trackAtVtx](const auto fitTrk) {
                       return trackAtVtx.originalParams == fitTrk;
                     });
    if (foundFitIter != perigeesToFit.end()) {
      perigeesToFit.erase(foundFitIter);
    } else {
      ACTS_WARNING("Track trackAtVtx not found in perigeesToFit!");
    }
  }  // end iteration over tracksAtVertex

  ACTS_DEBUG("After removal of tracks belonging to vertex, "
             << seedTracks.size() << " seed tracks left.");

  // Now start considering outliers
  // perigeesToFit that are left here were below
  // m_cfg.cutOffTrackWeight threshold and are hence outliers
  ACTS_DEBUG("Number of outliers: " << perigeesToFit.size());

  for (const auto& myPerigeeToFit : perigeesToFit) {
    // calculate chi2 w.r.t. last fitted vertex
    auto result = getCompatibility(m_extractParameters(*myPerigeeToFit),
                                   myVertex, vertexingOptions, state);

    if (!result.ok()) {
      return result.error();
    }

    double chi2 = *result;

    // check if sufficiently compatible with last fitted vertex
    // (quite loose constraint)
    if (chi2 < m_cfg.maximumChi2cutForSeeding) {
      auto foundIter = std::find_if(seedTracks.begin(), seedTracks.end(),
                                    [&myPerigeeToFit](const auto seedTrk) {
                                      return myPerigeeToFit == seedTrk;
                                    });
      if (foundIter != seedTracks.end()) {
        // Remove track from seed tracks
        seedTracks.erase(foundIter);
      }

    } else {
      // Track not compatible with vertex
      // Remove track from current vertex
      auto foundIter =
          std::find_if(tracksAtVertex.begin(), tracksAtVertex.end(),
                       [&myPerigeeToFit](auto trk) {
                         return myPerigeeToFit == trk.originalParams;
                       });
      if (foundIter != tracksAtVertex.end()) {
        // Remove track from seed tracks
        tracksAtVertex.erase(foundIter);
      }
    }
  }

  // set updated (possibly with removed outliers) tracksAtVertex to vertex
  myVertex.setTracksAtVertex(tracksAtVertex);

  return {};
}

template <typename vfitter_t, typename sfinder_t>
Acts::Result<void>
Acts::IterativeVertexFinder<vfitter_t, sfinder_t>::fillPerigeesToFit(
    const std::vector<const InputTrack_t*>& perigeeList,
    const Vertex<InputTrack_t>& seedVertex,
    std::vector<const InputTrack_t*>& perigeesToFitOut,
    std::vector<const InputTrack_t*>& perigeesToFitSplitVertexOut,
    const VertexingOptions<InputTrack_t>& vertexingOptions,
    State& state) const {
  int numberOfTracks = perigeeList.size();

  // Count how many tracks are used for fit
  int count = 0;
  // Fill perigeesToFit vector with tracks compatible with seed
  for (const auto& sTrack : perigeeList) {
    if (numberOfTracks <= 2) {
      perigeesToFitOut.push_back(sTrack);
      ++count;
    } else if (numberOfTracks <= 4 && !m_cfg.createSplitVertices) {
      perigeesToFitOut.push_back(sTrack);
      ++count;
    } else if (numberOfTracks <= 4 * m_cfg.splitVerticesTrkInvFraction &&
               m_cfg.createSplitVertices) {
      // Only few tracks left, put them into fit regardless their position
      if (count % m_cfg.splitVerticesTrkInvFraction == 0) {
        perigeesToFitOut.push_back(sTrack);
        ++count;
      } else {
        perigeesToFitSplitVertexOut.push_back(sTrack);
        ++count;
      }
    }
    // still large amount of tracks available, check compatibility
    else {
      // check first that distance is not too large
      const BoundTrackParameters& sTrackParams = m_extractParameters(*sTrack);
      auto distanceRes = m_cfg.ipEst.calculate3dDistance(
          vertexingOptions.geoContext, sTrackParams, seedVertex.position(),
          state.ipState);
      if (!distanceRes.ok()) {
        return distanceRes.error();
      }

      if (!sTrackParams.covariance()) {
        return VertexingError::NoCovariance;
      }

      double error =
          sqrt((*(sTrackParams.covariance()))(eBoundLoc0, eBoundLoc0) +
               (*(sTrackParams.covariance()))(eBoundLoc1, eBoundLoc1));

      if (error == 0.) {
        ACTS_WARNING("Error is zero. Setting to 1.");
        error = 1.;
      }

      if (*distanceRes / error < m_cfg.significanceCutSeeding) {
        if (!m_cfg.createSplitVertices ||
            count % m_cfg.splitVerticesTrkInvFraction == 0) {
          perigeesToFitOut.push_back(sTrack);
          ++count;
        } else {
          perigeesToFitSplitVertexOut.push_back(sTrack);
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
    std::vector<const InputTrack_t*>& perigeesToFit,
    std::vector<const InputTrack_t*>& seedTracks,
    const std::vector<const InputTrack_t*>& /* origTracks */,
    const VertexingOptions<InputTrack_t>& vertexingOptions,
    State& state) const {
  int numberOfAddedTracks = 0;

  // iterate over all vertices and check if tracks need to be reassigned
  // to new (current) vertex
  for (auto& vertexIt : vertexCollection) {
    // tracks at vertexIt
    std::vector<TrackAtVertex<InputTrack_t>> tracksAtVertex = vertexIt.tracks();
    auto tracksBegin = tracksAtVertex.begin();
    auto tracksEnd = tracksAtVertex.end();

    for (auto tracksIter = tracksBegin; tracksIter != tracksEnd;) {
      // consider only tracks that are not too tightly assigned to other
      // vertex
      if (tracksIter->trackWeight > m_cfg.cutOffTrackWeightReassign) {
        tracksIter++;
        continue;
      }
      // use original perigee parameter of course
      BoundTrackParameters trackPerigee =
          m_extractParameters(*(tracksIter->originalParams));

      // compute compatibility
      auto resultNew = getCompatibility(trackPerigee, currentVertex,
                                        vertexingOptions, state);
      if (!resultNew.ok()) {
        return Result<bool>::failure(resultNew.error());
      }
      double chi2NewVtx = *resultNew;

      auto resultOld =
          getCompatibility(trackPerigee, vertexIt, vertexingOptions, state);
      if (!resultOld.ok()) {
        return Result<bool>::failure(resultOld.error());
      }
      double chi2OldVtx = *resultOld;

      ACTS_DEBUG("Compatibility to new vs old vertex: " << chi2NewVtx << " vs "
                                                        << chi2OldVtx);

      if (chi2NewVtx < chi2OldVtx) {
        perigeesToFit.push_back(tracksIter->originalParams);
        // origTrack was already deleted from seedTracks previously
        // (when assigned to old vertex)
        // add it now back to seedTracks to be able to consistently
        // delete it later
        // when all tracks used to fit current vertex are deleted
        seedTracks.push_back(tracksIter->originalParams);
        // seedTracks.push_back(*std::find_if(
        //     origTracks.begin(), origTracks.end(),
        //     [&trackPerigee, this](auto origTrack) {
        //       return trackPerigee == m_extractParameters(*origTrack);
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
  if (m_cfg.useBeamConstraint && !perigeesToFit.empty()) {
    auto fitResult = m_cfg.vertexFitter.fit(
        perigeesToFit, m_cfg.linearizer, vertexingOptions, state.fitterState);
    if (fitResult.ok()) {
      currentVertex = std::move(*fitResult);
    } else {
      return Result<bool>::success(false);
    }
  } else if (!m_cfg.useBeamConstraint && perigeesToFit.size() > 1) {
    auto fitResult = m_cfg.vertexFitter.fit(
        perigeesToFit, m_cfg.linearizer, vertexingOptions, state.fitterState);
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

  bool isGoodVertex =
      ((!m_cfg.useBeamConstraint && ndf > 0 && nTracksAtVertex >= 2) ||
       (m_cfg.useBeamConstraint && ndf > 3 && nTracksAtVertex >= 2));

  if (!isGoodVertex) {
    removeAllTracks(perigeesToFit, seedTracks);

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
