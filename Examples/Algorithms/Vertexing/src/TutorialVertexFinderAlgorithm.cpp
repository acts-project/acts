// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Vertexing/TutorialVertexFinderAlgorithm.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/Surfaces/PerigeeSurface.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFinder.hpp"
#include "Acts/Vertexing/AdaptiveMultiVertexFitter.hpp"
#include "Acts/Vertexing/HelicalTrackLinearizer.hpp"
#include "Acts/Vertexing/ImpactPointEstimator.hpp"
#include "Acts/Vertexing/LinearizedTrack.hpp"
#include "Acts/Vertexing/TrackDensityVertexFinder.hpp"
#include "Acts/Vertexing/Vertex.hpp"
#include "Acts/Vertexing/VertexFinderConcept.hpp"
#include "Acts/Vertexing/VertexingOptions.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include "VertexingHelpers.hpp"

ActsExamples::TutorialVertexFinderAlgorithm::TutorialVertexFinderAlgorithm(
    const Config& cfg, Acts::Logging::Level lvl)
    : ActsExamples::BareAlgorithm("TutorialVertexFinder", lvl), m_cfg(cfg) {
  if (m_cfg.inputTrackParameters.empty()) {
    throw std::invalid_argument("Missing input track parameters collection");
  }
  if (m_cfg.outputProtoVertices.empty()) {
    throw std::invalid_argument("Missing output proto vertices collection");
  }
}

ActsExamples::ProcessCode ActsExamples::TutorialVertexFinderAlgorithm::execute(
    const ActsExamples::AlgorithmContext& ctx) const {
  // retrieve input tracks and convert into the expected format
  const auto& inputTrackParameters =
      ctx.eventStore.get<TrackParametersContainer>(m_cfg.inputTrackParameters);
  const auto& inputTrackPointers =
      makeTrackParametersPointerContainer(inputTrackParameters);
  //* Do not change the code above this line *//

  //* Remove following 2 lines. Only here to suppress unused variable errors *//
  (void)(inputTrackParameters);
  (void)(inputTrackPointers);

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
  return ActsExamples::ProcessCode::SUCCESS;
}
