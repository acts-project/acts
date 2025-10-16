// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Material/IMaterialDecorator.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPlugins/Root/RootMaterialDecorator.hpp"
#include "ActsPlugins/Root/RootMaterialMapIo.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>
#include <vector>

#include <TGeoManager.h>
#include <TGeoVolume.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(ActsPluginsPythonBindingsRoot, root) {
  using namespace Acts;
  using namespace ActsPlugins;

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
}