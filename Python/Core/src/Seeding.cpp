// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Seeding/EstimateTrackParamsFromSeed.hpp"
#include "Acts/Seeding/SeedConfirmationRangeConfig.hpp"
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
    auto c = py::class_<Acts::SeedConfirmationRangeConfig>(
                 m, "SeedConfirmationRangeConfig")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, zMinSeedConf, zMaxSeedConf, rMaxSeedConf,
                       nTopForLargeR, nTopForSmallR, seedConfMinBottomRadius,
                       seedConfMaxZOrigin, minImpactSeedConf);
    patchKwargsConstructor(c);
  }

  m.def(
      "estimateTrackParamsFromSeed",
      [](const Vector3& sp0, double t0, const Vector3& sp1, const Vector3& sp2,
         const Vector3& bField) {
        return estimateTrackParamsFromSeed(sp0, t0, sp1, sp2, bField);
      },
      py::arg("sp0"), py::arg("t0"), py::arg("sp1"), py::arg("sp2"),
      py::arg("bField"));

  m.def(
      "estimateTrackTangentsFromSeed",
      [](const Vector3& sp0, const Vector3& sp1, const Vector3& sp2,
         const Vector3& bField) -> std::tuple<Vector3, Vector3, Vector3> {
        Vector3 tangent0, tangent1, tangent2;
        estimateTrackTangentsFromSeed(sp0, sp1, sp2, bField, tangent0, tangent1,
                                      tangent2);
        return {tangent0, tangent1, tangent2};
      },
      py::arg("sp0"), py::arg("sp1"), py::arg("sp2"), py::arg("bField"));
}

}  // namespace ActsPython
