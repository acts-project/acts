// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Material/BinnedSurfaceMaterialAccumulater.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/IMaterialDecorator.hpp"
#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/IVolumeMaterial.hpp"
#include "Acts/Material/IntersectionMaterialAssigner.hpp"
#include "Acts/Material/MaterialMapper.hpp"
#include "Acts/Material/MaterialValidater.hpp"
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

    py::class_<ProtoGridSurfaceMaterial, ISurfaceMaterial,
               std::shared_ptr<ProtoGridSurfaceMaterial>>(
        m, "ProtoGridSurfaceMaterial");

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
    auto cls =
        py::class_<SurfaceMaterialMapper,
                   std::shared_ptr<SurfaceMaterialMapper>>(
            m, "SurfaceMaterialMapper")
            .def(py::init([](const SurfaceMaterialMapper::Config& config,
                             SurfaceMaterialMapper::StraightLinePropagator prop,
                             Logging::Level level) {
                   return std::make_shared<SurfaceMaterialMapper>(
                       config, std::move(prop),
                       getDefaultLogger("SurfaceMaterialMapper", level));
                 }),
                 py::arg("config"), py::arg("propagator"), py::arg("level"));

    auto c = py::class_<SurfaceMaterialMapper::Config>(cls, "Config")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, etaRange, emptyBinCorrection, mapperDebugOutput,
                       computeVariance);
  }

  {
    auto cls =
        py::class_<VolumeMaterialMapper, std::shared_ptr<VolumeMaterialMapper>>(
            m, "VolumeMaterialMapper")
            .def(py::init([](const VolumeMaterialMapper::Config& config,
                             VolumeMaterialMapper::StraightLinePropagator prop,
                             Logging::Level level) {
                   return std::make_shared<VolumeMaterialMapper>(
                       config, std::move(prop),
                       getDefaultLogger("VolumeMaterialMapper", level));
                 }),
                 py::arg("config"), py::arg("propagator"), py::arg("level"));

    auto c = py::class_<VolumeMaterialMapper::Config>(cls, "Config")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, mappingStep);
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
    py::class_<ISurfaceMaterialAccumulater,
               std::shared_ptr<ISurfaceMaterialAccumulater>>(
        m, "ISurfaceMaterialAccumulater");
  }

  {
    auto bsma =
        py::class_<BinnedSurfaceMaterialAccumulater,
                   ISurfaceMaterialAccumulater,
                   std::shared_ptr<BinnedSurfaceMaterialAccumulater>>(
            m, "BinnedSurfaceMaterialAccumulater")
            .def(
                py::init(
                    [](const BinnedSurfaceMaterialAccumulater::Config& config,
                       Logging::Level level) {
                      return std::make_shared<BinnedSurfaceMaterialAccumulater>(
                          config,
                          getDefaultLogger("BinnedSurfaceMaterialAccumulater",
                                           level));
                    }),
                py::arg("config"), py::arg("level"))
            .def("createState", &BinnedSurfaceMaterialAccumulater::createState)
            .def("accumulate", &BinnedSurfaceMaterialAccumulater::accumulate)
            .def("finalizeMaterial",
                 &BinnedSurfaceMaterialAccumulater::finalizeMaterial);

    auto c =
        py::class_<BinnedSurfaceMaterialAccumulater::Config>(bsma, "Config")
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
    ACTS_PYTHON_STRUCT(c, assignmentFinder, surfaceMaterialAccumulater);
  }

  {
    auto mvc =
        py::class_<MaterialValidater, std::shared_ptr<MaterialValidater>>(
            m, "MaterialValidater")
            .def(py::init([](const MaterialValidater::Config& config,
                             Logging::Level level) {
                   return std::make_shared<MaterialValidater>(
                       config, getDefaultLogger("MaterialValidater", level));
                 }),
                 py::arg("config"), py::arg("level"))
            .def("recordMaterial", &MaterialValidater::recordMaterial);

    auto c =
        py::class_<MaterialValidater::Config>(mvc, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, materialAssigner);
  }
}

}  // namespace ActsPython
