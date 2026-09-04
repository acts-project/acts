// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPlugins/Geant4/Geant4DetectorElement.hpp"
#include "ActsPlugins/Geant4/Geant4DetectorSurfaceFactory.hpp"
#include "ActsPlugins/Geant4/Geant4PhysicalVolumeSelectors.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>
#include <string>

#include <pybind11/pybind11.h>

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(ActsPluginsPythonBindingsGeant4, geant4) {
  using namespace Acts;
  using namespace ActsPlugins;

  // Detector Element
  py::class_<Geant4DetectorElement, std::shared_ptr<Geant4DetectorElement>,
             SurfacePlacementBase>(geant4, "Geant4DetectorElement")
      .def("surface", [](const Geant4DetectorElement& self) {
        return self.surface().getSharedPtr();
      });

  // Surface creation
  {
    using ISelector = IGeant4PhysicalVolumeSelector;
    auto is = py::class_<ISelector, std::shared_ptr<ISelector>>(
        geant4, "IVolumeSelector");

    using NameSelector = Geant4PhysicalVolumeSelectors::NameSelector;
    auto ns = py::class_<NameSelector, std::shared_ptr<NameSelector>>(
                  geant4, "VolumeNameSelector", is)
                  .def(py::init<const std::vector<std::string>&, bool>());

    using Factory = Geant4DetectorSurfaceFactory;
    auto o = py::class_<Factory::Options>(geant4, "SurfaceFactoryOptions")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(o, scaleConversion, convertMaterial,
                       convertedMaterialThickness, sensitiveSurfaceSelector,
                       passiveSurfaceSelector);
  }
}
