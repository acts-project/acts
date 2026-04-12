// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Material/BinnedSurfaceMaterialAccumulator.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/IMaterialDecorator.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Material/IntersectionMaterialAssigner.hpp"
#include "Acts/Material/MaterialMapper.hpp"
#include "Acts/Material/MaterialValidator.hpp"
#include "Acts/Material/PropagatorMaterialAssigner.hpp"
#include "Acts/Material/ProtoSurfaceMaterial.hpp"
#include "Acts/Material/SurfaceMaterialMapper.hpp"
#include "Acts/Material/VolumeMaterialMapper.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <array>
#include <map>
#include <memory>
#include <tuple>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace Acts {
class TrackingGeometry;
}  // namespace Acts

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;

namespace ActsPython {

/// @brief Add the material bindings to a module.
/// @param m the module to add the bindings to
void addMaterial(py::module_& m) {
  {
    py::class_<ISurfaceMaterial, std::shared_ptr<ISurfaceMaterial>>(
        m, "ISurfaceMaterial")
        .def("toString", &ISurfaceMaterial::toString);

    py::class_<ProtoSurfaceMaterial, ISurfaceMaterial,
               std::shared_ptr<ProtoSurfaceMaterial>>(m,
                                                      "ProtoSurfaceMaterial");

    py::class_<HomogeneousSurfaceMaterial, ISurfaceMaterial,
               std::shared_ptr<HomogeneousSurfaceMaterial>>(
        m, "HomogeneousSurfaceMaterial");

    py::class_<IVolumeMaterial, std::shared_ptr<IVolumeMaterial>>(
        m, "IVolumeMaterial");
  }

  {
    py::class_<IMaterialDecorator, std::shared_ptr<IMaterialDecorator>>(
        m, "IMaterialDecorator")
        .def("decorate", py::overload_cast<Surface&>(
                             &IMaterialDecorator::decorate, py::const_));
  }

  {
    py::class_<IAssignmentFinder, std::shared_ptr<IAssignmentFinder>>(
        m, "IAssignmentFinder");
  }

  {
    auto isma =
        py::class_<IntersectionMaterialAssigner, IAssignmentFinder,
                   std::shared_ptr<IntersectionMaterialAssigner>>(
            m, "IntersectionMaterialAssigner")
            .def(py::init([](const IntersectionMaterialAssigner::Config& config,
                             Logging::Level level) {
                   return std::make_shared<IntersectionMaterialAssigner>(
                       config,
                       getDefaultLogger("IntersectionMaterialAssigner", level));
                 }),
                 py::arg("config"), py::arg("level"))
            .def("assignmentCandidates",
                 &IntersectionMaterialAssigner::assignmentCandidates);

    auto c = py::class_<IntersectionMaterialAssigner::Config>(isma, "Config")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, surfaces, trackingVolumes);
  }

  {
    py::class_<ISurfaceMaterialAccumulator,
               std::shared_ptr<ISurfaceMaterialAccumulator>>(
        m, "ISurfaceMaterialAccumulator");
  }

  {
    auto bsma =
        py::class_<BinnedSurfaceMaterialAccumulator,
                   ISurfaceMaterialAccumulator,
                   std::shared_ptr<BinnedSurfaceMaterialAccumulator>>(
            m, "BinnedSurfaceMaterialAccumulator")
            .def(
                py::init(
                    [](const BinnedSurfaceMaterialAccumulator::Config& config,
                       Logging::Level level) {
                      return std::make_shared<BinnedSurfaceMaterialAccumulator>(
                          config,
                          getDefaultLogger("BinnedSurfaceMaterialAccumulator",
                                           level));
                    }),
                py::arg("config"), py::arg("level"))
            .def("createState", &BinnedSurfaceMaterialAccumulator::createState)
            .def("accumulate", &BinnedSurfaceMaterialAccumulator::accumulate)
            .def("finalizeMaterial",
                 &BinnedSurfaceMaterialAccumulator::finalizeMaterial);

    auto c =
        py::class_<BinnedSurfaceMaterialAccumulator::Config>(bsma, "Config")
            .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, emptyBinCorrection, materialSurfaces);
  }

  {
    auto mm = py::class_<MaterialMapper, std::shared_ptr<MaterialMapper>>(
                  m, "MaterialMapper")
                  .def(py::init([](const MaterialMapper::Config& config,
                                   Logging::Level level) {
                         return std::make_shared<MaterialMapper>(
                             config, getDefaultLogger("MaterialMapper", level));
                       }),
                       py::arg("config"), py::arg("level"));

    auto c = py::class_<MaterialMapper::Config>(mm, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, assignmentFinder, surfaceMaterialAccumulator);
  }

  {
    auto mvc =
        py::class_<MaterialValidator, std::shared_ptr<MaterialValidator>>(
            m, "MaterialValidator")
            .def(py::init([](const MaterialValidator::Config& config,
                             Logging::Level level) {
                   return std::make_shared<MaterialValidator>(
                       config, getDefaultLogger("MaterialValidator", level));
                 }),
                 py::arg("config"), py::arg("level"))
            .def("recordMaterial", &MaterialValidator::recordMaterial);

    auto c =
        py::class_<MaterialValidator::Config>(mvc, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, materialAssigner);
  }
}

}  // namespace ActsPython
