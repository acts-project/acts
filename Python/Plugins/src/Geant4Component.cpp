// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Detector/Detector.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Plugins/Geant4/Geant4DetectorElement.hpp"
#include "Acts/Plugins/Geant4/Geant4DetectorSurfaceFactory.hpp"
#include "Acts/Plugins/Geant4/Geant4PhysicalVolumeSelectors.hpp"
#include "Acts/Surfaces/SurfaceVisitorConcept.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPython/Utilities/Context.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <algorithm>
#include <memory>
#include <ranges>
#include <string>
#include <tuple>
#include <unordered_map>
#include <utility>
#include <vector>

#include <G4RunManager.hh>
#include <G4Transform3D.hh>
#include <G4UserEventAction.hh>
#include <G4UserRunAction.hh>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsPython;


PYBIND11_MODULE(ActsPluginsPythonBindingsGeant4, mod) {

  {
    using ISelector = IGeant4PhysicalVolumeSelector;
    auto is = py::class_<ISelector, std::shared_ptr<ISelector>>(
        mod, "IVolumeSelector");

    using NameSelector = Geant4PhysicalVolumeSelectors::NameSelector;
    auto ns = py::class_<NameSelector, std::shared_ptr<NameSelector>>(
                  mod, "VolumeNameSelector", is)
                  .def(py::init<const std::vector<std::string>&, bool>());

    using Factory = Geant4DetectorSurfaceFactory;
    auto o = py::class_<Factory::Options>(mod, "SurfaceFactoryOptions")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(o, scaleConversion, convertMaterial,
                       convertedMaterialThickness, sensitiveSurfaceSelector,
                       passiveSurfaceSelector);
  }

  ActsPython::Context ctx;
  ctx.modules["geant4"] = mod;
}
