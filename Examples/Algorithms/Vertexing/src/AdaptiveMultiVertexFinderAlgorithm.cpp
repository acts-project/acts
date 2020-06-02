// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Vertexing/AdaptiveMultiVertexFinderAlgorithm.hpp"

#include "ACTFW/Framework/RandomNumbers.hpp"
#include "ACTFW/TruthTracking/VertexAndTracks.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/Units.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFinder.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFitter.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"
#include "Acts/Vertexing/TrackDensityVertexFinder.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexFinderConcept.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"

FWE::AdaptiveMultiVertexFinderAlgorithm::AdaptiveMultiVertexFinderAlgorithm(
    const Config& cfg, Acts::Logging::Level level)
    : FW::BareAlgorithm("AMVF Algorithm", level), m_cfg(cfg) {}

/// @brief Algorithm that receives all selected tracks from an event
/// and finds and fits its vertices
FW::ProcessCode FWE::AdaptiveMultiVertexFinderAlgorithm::execute(
    const FW::AlgorithmContext& ctx) const {
  //////////////////////////////////////////////
  /* Full tutorial example code for reference */
  //////////////////////////////////////////////

  using namespace Acts::UnitLiterals;
  // Get the input track collection
  auto allTracks = getInputTrackCollection(ctx);
  // Create vector of track pointers for vertexing
  std::vector<const Acts::BoundParameters*> inputTrackPtrCollection;
  for (const auto& trk : allTracks) {
    inputTrackPtrCollection.push_back(&trk);
  }

  // Set up the magnetic field
  Acts::ConstantBField bField(Acts::Vector3D(0., 0., 2_T));

  // Set up EigenStepper
  Acts::EigenStepper<Acts::ConstantBField> stepper(bField);
  // Set up the propagator
  using Propagator = Acts::Propagator<Acts::EigenStepper<Acts::ConstantBField>>;
  auto propagator = std::make_shared<Propagator>(stepper);

  // Set up ImpactPointEstimator
  using IPEstimator =
      Acts::ImpactPointEstimator<Acts::BoundParameters, Propagator>;
  IPEstimator::Config ipEstimatorCfg(bField, propagator);
  IPEstimator ipEstimator(ipEstimatorCfg);

  // Set up the helical track linearizer
  using Linearizer = Acts::HelicalTrackLinearizer<Propagator>;
  Linearizer::Config ltConfig(bField, propagator);
  Linearizer linearizer(ltConfig);

  // Set up deterministic annealing with user-defined temperatures
  std::vector<double> temperatures{8.0, 4.0, 2.0, 1.4142136, 1.2247449, 1.0};
  Acts::AnnealingUtility::Config annealingConfig(temperatures);
  Acts::AnnealingUtility annealingUtility(annealingConfig);

  // Set up the vertex fitter with user-defined annealing
  using Fitter =
      Acts::AdaptiveMultiVertexFitter<Acts::BoundParameters, Linearizer>;
  Fitter::Config fitterCfg(ipEstimator);
  fitterCfg.annealingTool = annealingUtility;
  Fitter fitter(fitterCfg);

  // Set up the vertex seed finder
  using SeedFinder =
      Acts::TrackDensityVertexFinder<Fitter, Acts::GaussianTrackDensity>;
  SeedFinder seedFinder;

  // The vertex finder type
  using Finder = Acts::AdaptiveMultiVertexFinder<Fitter, SeedFinder>;

  Finder::Config finderConfig(std::move(fitter), seedFinder, ipEstimator,
                              linearizer);
  // We do not want to use a beamspot constraint here
  finderConfig.useBeamSpotConstraint = false;

  // Instantiate the finder
  Finder finder(finderConfig);
  // The vertex finder state
  Finder::State state;

  // Default vertexing options, this is where e.g. a constraint could be set
  using VertexingOptions = Acts::VertexingOptions<Acts::BoundParameters>;
  VertexingOptions finderOpts(ctx.geoContext, ctx.magFieldContext);

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
    ACTS_ERROR("Error in vertex finder: " << res.error().message());
  }

  return FW::ProcessCode::SUCCESS;
}

std::vector<Acts::BoundParameters>
FWE::AdaptiveMultiVertexFinderAlgorithm::getInputTrackCollection(
    const FW::AlgorithmContext& ctx) const {
  // Setup containers
  const auto& input = ctx.eventStore.get<std::vector<FW::VertexAndTracks>>(
      m_cfg.trackCollection);
  std::vector<Acts::BoundParameters> inputTrackCollection;

  for (auto& vertexAndTracks : input) {
    inputTrackCollection.insert(inputTrackCollection.end(),
                                vertexAndTracks.tracks.begin(),
                                vertexAndTracks.tracks.end());
  }

  return inputTrackCollection;
}
