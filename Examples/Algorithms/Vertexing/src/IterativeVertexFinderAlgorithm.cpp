// This file is part of the Acts project.
//
// Copyright (C) 2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Vertexing/IterativeVertexFinderAlgorithm.hpp"

#include <Acts/Geometry/GeometryContext.hpp>
#include <Acts/MagneticField/MagneticFieldContext.hpp>
#include <iostream>

#include "ACTFW/Framework/RandomNumbers.hpp"
#include "ACTFW/TruthTracking/VertexAndTracks.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Vertexing/FullBilloirVertexFitter.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/IterativeVertexFinder.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexFinderConcept.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"
#include "Acts/Vertexing/ZScanVertexFinder.hpp"

FWE::IterativeVertexFinderAlgorithm::IterativeVertexFinderAlgorithm(
    const Config& cfg, Acts::Logging::Level level)
    : FW::BareAlgorithm("VertexFinding", level), m_cfg(cfg) {}

/// @brief Algorithm that receives all selected tracks from an event
/// and finds and fits its vertices
FW::ProcessCode FWE::IterativeVertexFinderAlgorithm::execute(
    const FW::AlgorithmContext& ctx) const {
  using MagneticField = Acts::ConstantBField;
  using Stepper = Acts::EigenStepper<MagneticField>;
  using Propagator = Acts::Propagator<Stepper>;
  using PropagatorOptions = Acts::PropagatorOptions<>;
  using TrackParameters = Acts::BoundParameters;
  using Linearizer = Acts::HelicalTrackLinearizer<Propagator>;
  using VertexFitter =
      Acts::FullBilloirVertexFitter<TrackParameters, Linearizer>;
  using ImpactPointEstimator =
      Acts::ImpactPointEstimator<TrackParameters, Propagator>;
  using VertexSeeder = Acts::ZScanVertexFinder<VertexFitter>;
  using VertexFinder = Acts::IterativeVertexFinder<VertexFitter, VertexSeeder>;
  using VertexFinderOptions = Acts::VertexingOptions<TrackParameters>;

  static_assert(Acts::VertexFinderConcept<VertexSeeder>,
                "VertexSeeder does not fulfill vertex finder concept.");
  static_assert(Acts::VertexFinderConcept<VertexFinder>,
                "VertexFinder does not fulfill vertex finder concept.");

  // Set up the magnetic field
  MagneticField bField(m_cfg.bField);
  // Set up propagator with void navigator
  auto propagator = std::make_shared<Propagator>(Stepper(bField));
  PropagatorOptions propagatorOpts(ctx.geoContext, ctx.magFieldContext);
  // Setup the vertex fitter
  VertexFitter::Config vertexFitterCfg;
  VertexFitter vertexFitter(std::move(vertexFitterCfg));
  // Setup the track linearizer
  Linearizer::Config linearizerCfg(bField, propagator);
  Linearizer linearizer(std::move(linearizerCfg));
  // Setup the seed finder
  ImpactPointEstimator::Config ipEstCfg(bField, propagator);
  ImpactPointEstimator ipEst(std::move(ipEstCfg));
  VertexSeeder::Config seederCfg(ipEst);
  VertexSeeder seeder(std::move(seederCfg));
  // Set up the actual vertex finder
  VertexFinder::Config finderCfg(std::move(vertexFitter), std::move(linearizer),
                                 std::move(seeder), ipEst);
  finderCfg.maxVertices = 200;
  finderCfg.reassignTracksAfterFirstFit = true;
  VertexFinder finder(finderCfg);
  VertexFinder::State state;
  VertexFinderOptions finderOpts(ctx.geoContext, ctx.magFieldContext);

  // Setup containers
  const auto& input = ctx.eventStore.get<std::vector<FW::VertexAndTracks>>(
      m_cfg.trackCollection);
  std::vector<Acts::BoundParameters> inputTrackCollection;

  int counte = 0;
  for (auto& bla : input) {
    counte += bla.tracks.size();
  }

  ACTS_INFO("Truth vertices in event: " << input.size());

  for (auto& vertexAndTracks : input) {
    ACTS_INFO("\t True vertex at ("
              << vertexAndTracks.vertex.position().x() << ","
              << vertexAndTracks.vertex.position().y() << ","
              << vertexAndTracks.vertex.position().z() << ") with "
              << vertexAndTracks.tracks.size() << " tracks.");
    inputTrackCollection.insert(inputTrackCollection.end(),
                                vertexAndTracks.tracks.begin(),
                                vertexAndTracks.tracks.end());
  }

  std::vector<const Acts::BoundParameters*> inputTrackPtrCollection;
  for (const auto& trk : inputTrackCollection) {
    inputTrackPtrCollection.push_back(&trk);
  }

  // Find vertices
  auto res = finder.find(inputTrackPtrCollection, finderOpts, state);

  if (res.ok()) {
    // Retrieve vertices found by vertex finder
    auto vertexCollection = *res;

    ACTS_INFO("Found " << vertexCollection.size() << " vertices in event.");

    unsigned int count = 0;
    for (const auto& vtx : vertexCollection) {
      ACTS_INFO("\t" << ++count << ". vertex at "
                     << "(" << vtx.position().x() << "," << vtx.position().y()
                     << "," << vtx.position().z() << ") with "
                     << vtx.tracks().size() << " tracks.");
    }
  } else {
    ACTS_ERROR("Error in vertex finder.");
  }

  return FW::ProcessCode::SUCCESS;
}
