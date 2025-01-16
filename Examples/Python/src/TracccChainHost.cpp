// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/EventData/Index.hpp"
#include "ActsExamples/Traccc/Host/TracccChainAlgorithm.hpp"

#include <memory>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

namespace Acts::Python {

void addTracccChainHost(Context& ctx) {
  auto m = ctx.get("traccc");

  ACTS_PYTHON_DECLARE_ALGORITHM(
      ActsExamples::Traccc::Host::TracccChainAlgorithm, m,
      "ReconstructionChainHostAlgorithm", inputCells, inputMeasurements,
      inputSpacePoints, inputSeeds, outputSpacePoints, outputSeeds,
      outputTracks, enableAmbiguityResolution, detrayStore, trackingGeometry,
      field, digitizationConfigs, seedfinderConfig, seedfilterConfig,
      findingConfig, fittingConfig, ambiguityResolutionConfig);
}

}  // namespace Acts::Python
