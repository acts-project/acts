// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTFW/Vertexing/TutorialAMVFAlgorithm.hpp"

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

FWE::TutorialAMVFAlgorithm::TutorialAMVFAlgorithm(const Config& cfg,
                                                  Acts::Logging::Level level)
    : FW::BareAlgorithm("Tutorial AMVF Algorithm", level), m_cfg(cfg) {}

/// @brief Tutorial algorithm that receives all selected tracks from an event
/// and finds and fits its vertices using the AMVF
FW::ProcessCode FWE::TutorialAMVFAlgorithm::execute(
    const FW::AlgorithmContext& ctx) const {
  using namespace Acts::UnitLiterals;

  // Get the input track collection
  auto allTracks = getInputTrackCollection(ctx);
  // Create vector of track pointers for vertexing
  std::vector<const Acts::BoundParameters*> inputTrackPtrCollection;
  for (const auto& trk : allTracks) {
    inputTrackPtrCollection.push_back(&trk);
  }
  //* Do not change the code above this line *//

  //////////////////////////////////////////////////////////////////////////////
  /*****   Note: This is a skeleton file to be filled with tutorial code  *****/
  /*****   provided in the Acts Docs - Vertexing section under the link:  *****/
  /* https://acts.readthedocs.io/en/latest/howto/setup_and_run_vertexing.html */
  /*** or in the Acts repository in  docs/howto/setup_and_run_vertexing.md  ***/
  //////////////////////////////////////////////////////////////////////////////

  ///////////////////////////////////////////////
  /*     Add the tutorial example code here    */
  ///////////////////////////////////////////////
  ///////////////////////////////////////////////
  /*  For the full tutorial code please refer  */
  /* to AdaptiveMultiVertexFinderAlgorithm.cpp */
  ///////////////////////////////////////////////

  //* Do not change the code below this line *//
  return FW::ProcessCode::SUCCESS;
}

std::vector<Acts::BoundParameters>
FWE::TutorialAMVFAlgorithm::getInputTrackCollection(
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
