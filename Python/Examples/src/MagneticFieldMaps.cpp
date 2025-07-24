// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/MagneticField/MagneticField.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/MagneticField/BFieldMapUtils.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/MagneticField/MultiRangeBField.hpp"
#include "Acts/MagneticField/NullBField.hpp"
#include "Acts/MagneticField/SolenoidBField.hpp"
#include "ActsExamples/MagneticField/FieldMapRootIo.hpp"
#include "ActsExamples/MagneticField/FieldMapTextIo.hpp"
#include "ActsPython/Utilities/Context.hpp"

#include <array>
#include <cstddef>
#include <filesystem>
#include <memory>
#include <stdexcept>
#include <string>
#include <tuple>
#include <utility>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ActsExamples;

namespace ActsPython {

void addMagneticFieldMaps(Context& ctx) {
  auto [mex, prop] = ctx.get("examples", "propagation");

  py::class_<detail::InterpolatedMagneticField2,
             Acts::InterpolatedMagneticField, Acts::MagneticFieldProvider,
             std::shared_ptr<detail::InterpolatedMagneticField2>>(
      mex, "InterpolatedMagneticField2");

  py::class_<detail::InterpolatedMagneticField3,
             Acts::InterpolatedMagneticField, Acts::MagneticFieldProvider,
             std::shared_ptr<detail::InterpolatedMagneticField3>>(
      mex, "InterpolatedMagneticField3");

  mex.def(
      "MagneticFieldMapXyz",
      [](const std::string& filename, const std::string& tree,
         double lengthUnit, double BFieldUnit, bool firstOctant) {
        const std::filesystem::path file = filename;

        auto mapBins = [](std::array<std::size_t, 3> bins,
                          std::array<std::size_t, 3> sizes) {
          return (bins[0] * (sizes[1] * sizes[2]) + bins[1] * sizes[2] +
                  bins[2]);
        };

        if (file.extension() == ".root") {
          auto map = makeMagneticFieldMapXyzFromRoot(
              std::move(mapBins), file.native(), tree, lengthUnit, BFieldUnit,
              firstOctant);
          return std::make_shared<
              detail::InterpolatedMagneticField3>(std::move(map));
        } else if (file.extension() == ".txt") {
          auto map = makeMagneticFieldMapXyzFromText(
              std::move(mapBins), file.native(), lengthUnit, BFieldUnit,
              firstOctant);
          return std::make_shared<
              detail::InterpolatedMagneticField3>(std::move(map));
        } else {
          throw std::runtime_error("Unsupported magnetic field map file type");
        }
      },
      py::arg("file"), py::arg("tree") = "bField",
      py::arg("lengthUnit") = Acts::UnitConstants::mm,
      py::arg("BFieldUnit") = Acts::UnitConstants::T,
      py::arg("firstOctant") = false);

  mex.def(
      "MagneticFieldMapRz",
      [](const std::string& filename, const std::string& tree,
         double lengthUnit, double BFieldUnit, bool firstQuadrant) {
        const std::filesystem::path file = filename;

        auto mapBins = [](std::array<std::size_t, 2> bins,
                          std::array<std::size_t, 2> sizes) {
          return (bins[1] * sizes[0] + bins[0]);
        };

        if (file.extension() == ".root") {
          auto map = makeMagneticFieldMapRzFromRoot(
              std::move(mapBins), file.native(), tree, lengthUnit, BFieldUnit,
              firstQuadrant);
          return std::make_shared<
              detail::InterpolatedMagneticField2>(std::move(map));
        } else if (file.extension() == ".txt") {
          auto map = makeMagneticFieldMapRzFromText(
              std::move(mapBins), file.native(), lengthUnit, BFieldUnit,
              firstQuadrant);
          return std::make_shared<
              detail::InterpolatedMagneticField2>(std::move(map));
        } else {
          throw std::runtime_error("Unsupported magnetic field map file type");
        }
      },
      py::arg("file"), py::arg("tree") = "bField",
      py::arg("lengthUnit") = Acts::UnitConstants::mm,
      py::arg("BFieldUnit") = Acts::UnitConstants::T,
      py::arg("firstQuadrant") = false);
}

}  // namespace ActsPython
