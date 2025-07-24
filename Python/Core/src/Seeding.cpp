// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/SeedConfirmationRangeConfig.hpp"
#include "Acts/Seeding/SeedFilterConfig.hpp"
#include "Acts/Seeding/SeedFinderConfig.hpp"
#include "Acts/Seeding/SeedFinderGbtsConfig.hpp"
#include "Acts/Seeding/SeedFinderOrthogonalConfig.hpp"
#include "Acts/Seeding/detail/CylindricalSpacePointGrid.hpp"
#include "ActsPython/Utilities/Macros.hpp"
#include "ActsPython/Utilities/Patchers.hpp"

#include <pybind11/eval.h>
#include <pybind11/pybind11.h>
#include <pybind11/pytypes.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Acts;

namespace ActsPython {

/// @brief This adds the classes from Core/Seeding to the python module
/// @param m is the pybind11 core module
void addSeeding(py::module_& m) {
  {
    using Config = SeedFilterConfig;
    auto c = py::class_<Config>(m, "SeedFilterConfig").def(py::init<>());
    ACTS_PYTHON_STRUCT(
        c, deltaInvHelixDiameter, impactWeightFactor, zOriginWeightFactor,
        compatSeedWeight, deltaRMin, maxSeedsPerSpM, compatSeedLimit,
        seedConfirmation, centralSeedConfirmationRange,
        forwardSeedConfirmationRange, useDeltaRorTopRadius, seedWeightIncrement,
        numSeedIncrement, maxSeedsPerSpMConf, maxQualitySeedsPerSpMConf);
    patchKwargsConstructor(c);
  }

  {
    using seedOptions = SeedFinderOptions;
    auto c = py::class_<seedOptions>(m, "SeedFinderOptions").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, beamPos, bFieldInZ);
    patchKwargsConstructor(c);
  }

  {
    using seedConf = SeedConfirmationRangeConfig;
    auto c = py::class_<seedConf>(m, "SeedConfirmationRangeConfig")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, zMinSeedConf, zMaxSeedConf, rMaxSeedConf,
                       nTopForLargeR, nTopForSmallR, seedConfMinBottomRadius,
                       seedConfMaxZOrigin, minImpactSeedConf);
    patchKwargsConstructor(c);
  }

  {
    using Config = CylindricalSpacePointGridConfig;
    auto c = py::class_<Config>(m, "SpacePointGridConfig").def(py::init<>());

    ACTS_PYTHON_STRUCT(c, minPt, rMax, zMax, zMin, phiMin, phiMax, deltaRMax,
                       cotThetaMax, phiBinDeflectionCoverage, maxPhiBins,
                       impactMax, zBinEdges);
    patchKwargsConstructor(c);
  }
  {
    using Options = CylindricalSpacePointGridOptions;
    auto c = py::class_<Options>(m, "SpacePointGridOptions").def(py::init<>());

    ACTS_PYTHON_STRUCT(c, bFieldInZ);
    patchKwargsConstructor(c);
  }
}
}  // namespace ActsPython
