// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/Vertexing/AdaptiveMultiVertexFinderAlgorithm.hpp"
#include "ActsExamples/Vertexing/HoughVertexFinderAlgorithm.hpp"
#include "ActsExamples/Vertexing/IterativeVertexFinderAlgorithm.hpp"
#include "ActsExamples/Vertexing/VertexFitterAlgorithm.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace ActsExamples;
using namespace Acts;

namespace ActsPython {

void addVertexing(Context& ctx) {
  using Seeder = AdaptiveMultiVertexFinderAlgorithm::SeedFinder;
  auto mex = ctx.get("examples");
  auto& m = ctx.get("main");

  py::enum_<Seeder>(m, "VertexSeedFinder")
      .value("TruthSeeder", Seeder::TruthSeeder)
      .value("GaussianSeeder", Seeder::GaussianSeeder)
      .value("AdaptiveGridSeeder", Seeder::AdaptiveGridSeeder);

  ACTS_PYTHON_DECLARE_ALGORITHM(
      AdaptiveMultiVertexFinderAlgorithm, mex,
      "AdaptiveMultiVertexFinderAlgorithm", inputTrackParameters,
      inputTruthParticles, inputTruthVertices, outputProtoVertices,
      outputVertices, seedFinder, bField, minWeight, doSmoothing, maxIterations,
      useTime, tracksMaxZinterval, initialVariances, doFullSplitting,
      tracksMaxSignificance, maxMergeVertexSignificance, spatialBinExtent,
      temporalBinExtent);

  ACTS_PYTHON_DECLARE_ALGORITHM(IterativeVertexFinderAlgorithm, mex,
                                "IterativeVertexFinderAlgorithm",
                                inputTrackParameters, outputProtoVertices,
                                outputVertices, bField, maxIterations);

  ACTS_PYTHON_DECLARE_ALGORITHM(VertexFitterAlgorithm, mex,
                                "VertexFitterAlgorithm", inputTrackParameters,
                                inputProtoVertices, outputVertices, bField,
                                doConstrainedFit, constraintPos, constraintCov);

  ACTS_PYTHON_DECLARE_ALGORITHM(HoughVertexFinderAlgorithm, mex,
                                "HoughVertexFinderAlgorithm", inputSpacepoints,
                                outputVertices, targetSPs, minAbsEta, maxAbsEta,
                                minHits, defVtxPosition);
}

}  // namespace ActsPython
