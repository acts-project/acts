// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Seeding/SeedConfirmationRangeConfig.hpp"
#include "Acts/Seeding/SeedFilterConfig.hpp"
#include "Acts/Seeding/detail/CylindricalSpacePointGrid.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;

using namespace Acts;

namespace ActsPython {

/// Add the track finding related bindings to the module
/// @param m The module to add the bindings to
void addSeeding(py::module_& m) {
  {
    using Config = Acts::SeedFilterConfig;
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
    auto c = py::class_<Acts::SeedFinderOptions>(m, "SeedFinderOptions")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, beamPos, bFieldInZ);
    patchKwargsConstructor(c);
  }

  {
    auto c = py::class_<Acts::SeedConfirmationRangeConfig>(
                 m, "SeedConfirmationRangeConfig")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, zMinSeedConf, zMaxSeedConf, rMaxSeedConf,
                       nTopForLargeR, nTopForSmallR, seedConfMinBottomRadius,
                       seedConfMaxZOrigin, minImpactSeedConf);
    patchKwargsConstructor(c);
  }

  {
    auto c = py::class_<Acts::CylindricalSpacePointGridConfig>(
                 m, "SpacePointGridConfig")
                 .def(py::init<>());

    ACTS_PYTHON_STRUCT(c, minPt, rMax, zMax, zMin, phiMin, phiMax, deltaRMax,
                       cotThetaMax, phiBinDeflectionCoverage, maxPhiBins,
                       impactMax, zBinEdges);
    patchKwargsConstructor(c);
  }
  {
    auto c = py::class_<Acts::CylindricalSpacePointGridOptions>(
                 m, "SpacePointGridOptions")
                 .def(py::init<>());

    ACTS_PYTHON_STRUCT(c, bFieldInZ);
    patchKwargsConstructor(c);
  }
}
}  // namespace ActsPython
