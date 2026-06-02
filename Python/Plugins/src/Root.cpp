// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Units.hpp"
#include "Acts/Material/IMaterialDecorator.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Root/RootMagneticFieldIo.hpp"
#include "ActsPlugins/Root/RootMaterialDecorator.hpp"
#include "ActsPlugins/Root/RootMaterialMapIo.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>

#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace py = pybind11;
using namespace pybind11::literals;
using namespace Acts;
using namespace ActsPlugins;

PYBIND11_MODULE(ActsPluginsPythonBindingsRoot, root) {
  {
    auto ac = py::class_<RootMaterialMapIo::Config>(root, "AccessorConfig")
                  .def(py::init<>());

    ACTS_PYTHON_STRUCT(ac, volumePrefix, portalPrefix, layerPrefix,
                       passivePrefix, sensitivePrefix, nBinsHistName,
                       axisDirHistName, axisBoundaryTypeHistName, indexHistName,
                       minRangeHistName, maxRangeHistName, thicknessHistName,
                       x0HistName, l0HistName, aHistName, zHistName,
                       rhoHistName);

    auto ao = py::class_<RootMaterialMapIo::Options>(root, "AccessorOptions")
                  .def(py::init<>());

    ACTS_PYTHON_STRUCT(ao, homogeneousMaterialTreeName, indexedMaterialTreeName,
                       folderSurfaceNameBase, folderVolumeNameBase,
                       indexedMaterial);

    auto rmd =
        py::class_<RootMaterialDecorator, IMaterialDecorator,
                   std::shared_ptr<RootMaterialDecorator>>(
            root, "RootMaterialDecorator")
            .def(py::init<RootMaterialDecorator::Config, Logging::Level>(),
                 py::arg("config"), py::arg("level"));

    using Config = RootMaterialDecorator::Config;
    auto c = py::class_<Config>(rmd, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT(c, accessorConfig, accessorOptions, fileName);
  }

  {
    root.def(
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
            return std::make_shared<decltype(map)>(std::move(map));
          } else {
            throw std::runtime_error(
                "Unsupported magnetic field map file type");
          }
        },
        py::arg("file"), py::arg("tree") = "bField",
        py::arg("lengthUnit") = UnitConstants::mm,
        py::arg("BFieldUnit") = UnitConstants::T,
        py::arg("firstOctant") = false);

    root.def(
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
            return std::make_shared<decltype(map)>(std::move(map));
          } else {
            throw std::runtime_error(
                "Unsupported magnetic field map file type");
          }
        },
        py::arg("file"), py::arg("tree") = "bField",
        py::arg("lengthUnit") = UnitConstants::mm,
        py::arg("BFieldUnit") = UnitConstants::T,
        py::arg("firstQuadrant") = false);
  }
}
