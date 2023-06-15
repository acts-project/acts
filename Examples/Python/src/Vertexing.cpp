// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Utilities/TypeTraits.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/Vertexing/AdaptiveMultiVertexFinderAlgorithm.hpp"
#include "ActsExamples/Vertexing/IterativeVertexFinderAlgorithm.hpp"
#include "ActsExamples/Vertexing/TutorialVertexFinderAlgorithm.hpp"
#include "ActsExamples/Vertexing/VertexFitterAlgorithm.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace Acts::Python {

void addVertexing(Context& ctx) {
  auto mex = ctx.get("examples");

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::AdaptiveMultiVertexFinderAlgorithm, mex,
      "AdaptiveMultiVertexFinderAlgorithm", inputTrackParameters,
      inputTrajectories, outputProtoVertices, outputVertices, bField);

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::IterativeVertexFinderAlgorithm,
                                mex, "IterativeVertexFinderAlgorithm",
                                inputTrackParameters, inputTrajectories,
                                outputProtoVertices, outputVertices, bField);

  ACTS_PYTHON_DECLARE_ALGORITHM(ActsExamples::TutorialVertexFinderAlgorithm,
                                mex, "TutorialVertexFinderAlgorithm",
                                inputTrackParameters, inputTrajectories,
                                outputProtoVertices, bField);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::VertexFitterAlgorithm, mex, "VertexFitterAlgorithm",
      inputTrackParameters, inputTrajectories, inputProtoVertices,
      outputVertices, bField, doConstrainedFit, constraintPos, constraintCov);
}

}  // namespace Acts::Python
