// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
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
#include "Acts/Plugins/Json/ActsJson.hpp"
#include "Acts/Plugins/Json/MaterialMapJsonConverter.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/MaterialMapping/CoreMaterialMapping.hpp"
#include "ActsExamples/MaterialMapping/MappingMaterialDecorator.hpp"
#include "ActsExamples/MaterialMapping/MaterialMapping.hpp"
#include "ActsExamples/MaterialMapping/MaterialValidation.hpp"
#include "ActsPython/Utilities/Context.hpp"
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
namespace ActsExamples {
class IAlgorithm;
}  // namespace ActsExamples

namespace py = pybind11;
using namespace pybind11::literals;

using namespace ActsExamples;

namespace ActsPython {
void addMaterial(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  {
    py::class_<Acts::ISurfaceMaterial, std::shared_ptr<Acts::ISurfaceMaterial>>(
        m, "ISurfaceMaterial")
        .def("toString", &Acts::ISurfaceMaterial::toString);

    py::class_<Acts::ProtoGridSurfaceMaterial, Acts::ISurfaceMaterial,
               std::shared_ptr<Acts::ProtoGridSurfaceMaterial>>(
        m, "ProtoGridSurfaceMaterial");

    py::class_<Acts::ProtoSurfaceMaterial, Acts::ISurfaceMaterial,
               std::shared_ptr<Acts::ProtoSurfaceMaterial>>(
        m, "ProtoSurfaceMaterial");

    py::class_<Acts::HomogeneousSurfaceMaterial, Acts::ISurfaceMaterial,
               std::shared_ptr<Acts::HomogeneousSurfaceMaterial>>(
        m, "HomogeneousSurfaceMaterial");

    py::class_<Acts::IVolumeMaterial, std::shared_ptr<Acts::IVolumeMaterial>>(
        m, "IVolumeMaterial");
  }

  {
    py::class_<Acts::IMaterialDecorator,
               std::shared_ptr<Acts::IMaterialDecorator>>(m,
                                                          "IMaterialDecorator")
        .def("decorate", py::overload_cast<Acts::Surface&>(
                             &Acts::IMaterialDecorator::decorate, py::const_));
  }

  {
    py::class_<Acts::MappingMaterialDecorator, Acts::IMaterialDecorator,
               std::shared_ptr<Acts::MappingMaterialDecorator>>(
        m, "MappingMaterialDecorator")
        .def(py::init<const Acts::TrackingGeometry&, Acts::Logging::Level, bool,
                      bool>(),
             py::arg("tGeometry"), py::arg("level"),
             py::arg("clearSurfaceMaterial") = true,
             py::arg("clearVolumeMaterial") = true)
        .def("binningMap", &Acts::MappingMaterialDecorator::binningMap)
        .def("setBinningMap", &Acts::MappingMaterialDecorator::setBinningMap);
  }

  {
    using Alg = ActsExamples::MaterialMapping;

    auto alg = py::class_<Alg, ActsExamples::IAlgorithm, std::shared_ptr<Alg>>(
                   mex, "MaterialMapping")
                   .def(py::init<const Alg::Config&, Acts::Logging::Level>(),
                        py::arg("config"), py::arg("level"))
                   .def("scoringParameters", &Alg::scoringParameters)
                   .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Alg::Config>(alg, "Config")
                 .def(py::init<const Acts::GeometryContext&,
                               const Acts::MagneticFieldContext&>());

    ACTS_PYTHON_STRUCT(c, inputMaterialTracks, mappingMaterialCollection,
                       materialSurfaceMapper, materialVolumeMapper,
                       materialWriters, trackingGeometry, geoContext,
                       magFieldContext);
  }

  {
    auto cls =
        py::class_<Acts::SurfaceMaterialMapper,
                   std::shared_ptr<Acts::SurfaceMaterialMapper>>(
            m, "SurfaceMaterialMapper")
            .def(
                py::init(
                    [](const Acts::SurfaceMaterialMapper::Config& config,
                       Acts::SurfaceMaterialMapper::StraightLinePropagator prop,
                       Acts::Logging::Level level) {
                      return std::make_shared<Acts::SurfaceMaterialMapper>(
                          config, std::move(prop),
                          Acts::getDefaultLogger("SurfaceMaterialMapper",
                                                 level));
                    }),
                py::arg("config"), py::arg("propagator"), py::arg("level"));

    auto c = py::class_<Acts::SurfaceMaterialMapper::Config>(cls, "Config")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, etaRange, emptyBinCorrection, mapperDebugOutput,
                       computeVariance);
  }

  {
    auto cls =
        py::class_<Acts::VolumeMaterialMapper,
                   std::shared_ptr<Acts::VolumeMaterialMapper>>(
            m, "VolumeMaterialMapper")
            .def(py::init(
                     [](const Acts::VolumeMaterialMapper::Config& config,
                        Acts::VolumeMaterialMapper::StraightLinePropagator prop,
                        Acts::Logging::Level level) {
                       return std::make_shared<Acts::VolumeMaterialMapper>(
                           config, std::move(prop),
                           Acts::getDefaultLogger("VolumeMaterialMapper",
                                                  level));
                     }),
                 py::arg("config"), py::arg("propagator"), py::arg("level"));

    auto c = py::class_<Acts::VolumeMaterialMapper::Config>(cls, "Config")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, mappingStep);
  }

  {
    py::class_<Acts::IAssignmentFinder,
               std::shared_ptr<Acts::IAssignmentFinder>>(m,
                                                         "IAssignmentFinder");
  }

  {
    auto isma =
        py::class_<Acts::IntersectionMaterialAssigner, Acts::IAssignmentFinder,
                   std::shared_ptr<Acts::IntersectionMaterialAssigner>>(
            m, "IntersectionMaterialAssigner")
            .def(py::init([](const Acts::IntersectionMaterialAssigner::Config&
                                 config,
                             Acts::Logging::Level level) {
                   return std::make_shared<Acts::IntersectionMaterialAssigner>(
                       config, Acts::getDefaultLogger(
                                   "IntersectionMaterialAssigner", level));
                 }),
                 py::arg("config"), py::arg("level"))
            .def("assignmentCandidates",
                 &Acts::IntersectionMaterialAssigner::assignmentCandidates);

    auto c =
        py::class_<Acts::IntersectionMaterialAssigner::Config>(isma, "Config")
            .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, surfaces, trackingVolumes, detectorVolumes);
  }

  {
    py::class_<Acts::ISurfaceMaterialAccumulater,
               std::shared_ptr<Acts::ISurfaceMaterialAccumulater>>(
        m, "ISurfaceMaterialAccumulater");
  }

  {
    auto bsma =
        py::class_<Acts::BinnedSurfaceMaterialAccumulater,
                   Acts::ISurfaceMaterialAccumulater,
                   std::shared_ptr<Acts::BinnedSurfaceMaterialAccumulater>>(
            m, "BinnedSurfaceMaterialAccumulater")
            .def(py::init(
                     [](const Acts::BinnedSurfaceMaterialAccumulater::Config&
                            config,
                        Acts::Logging::Level level) {
                       return std::make_shared<
                           Acts::BinnedSurfaceMaterialAccumulater>(
                           config,
                           Acts::getDefaultLogger(
                               "BinnedSurfaceMaterialAccumulater", level));
                     }),
                 py::arg("config"), py::arg("level"))
            .def("createState",
                 &Acts::BinnedSurfaceMaterialAccumulater::createState)
            .def("accumulate",
                 &Acts::BinnedSurfaceMaterialAccumulater::accumulate)
            .def("finalizeMaterial",
                 &Acts::BinnedSurfaceMaterialAccumulater::finalizeMaterial);

    auto c = py::class_<Acts::BinnedSurfaceMaterialAccumulater::Config>(
                 bsma, "Config")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, emptyBinCorrection, materialSurfaces);
  }

  {
    auto mm =
        py::class_<Acts::MaterialMapper, std::shared_ptr<Acts::MaterialMapper>>(
            m, "MaterialMapper")
            .def(py::init([](const Acts::MaterialMapper::Config& config,
                             Acts::Logging::Level level) {
                   return std::make_shared<Acts::MaterialMapper>(
                       config, Acts::getDefaultLogger("MaterialMapper", level));
                 }),
                 py::arg("config"), py::arg("level"));

    auto c = py::class_<Acts::MaterialMapper::Config>(mm, "Config")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, assignmentFinder, surfaceMaterialAccumulater);
  }

  {
    auto mmca = py::class_<CoreMaterialMapping, IAlgorithm,
                           std::shared_ptr<CoreMaterialMapping>>(
                    mex, "CoreMaterialMapping")
                    .def(py::init<const CoreMaterialMapping::Config&,
                                  Acts::Logging::Level>(),
                         py::arg("config"), py::arg("level"));

    auto c = py::class_<CoreMaterialMapping::Config>(mmca, "Config")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, inputMaterialTracks, mappedMaterialTracks,
                       unmappedMaterialTracks, materialMapper,
                       materiaMaplWriters);
  }

  {
    auto mvc =
        py::class_<Acts::MaterialValidater,
                   std::shared_ptr<Acts::MaterialValidater>>(
            m, "MaterialValidater")
            .def(py::init([](const Acts::MaterialValidater::Config& config,
                             Acts::Logging::Level level) {
                   return std::make_shared<Acts::MaterialValidater>(
                       config,
                       Acts::getDefaultLogger("MaterialValidater", level));
                 }),
                 py::arg("config"), py::arg("level"))
            .def("recordMaterial", &Acts::MaterialValidater::recordMaterial);

    auto c = py::class_<Acts::MaterialValidater::Config>(mvc, "Config")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, materialAssigner);
  }

  {
    auto mv = py::class_<MaterialValidation, IAlgorithm,
                         std::shared_ptr<MaterialValidation>>(
                  mex, "MaterialValidation")
                  .def(py::init<const MaterialValidation::Config&,
                                Acts::Logging::Level>(),
                       py::arg("config"), py::arg("level"))
                  .def("execute", &MaterialValidation::execute)
                  .def_property_readonly("config", &MaterialValidation::config);

    auto c =
        py::class_<MaterialValidation::Config>(mv, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, ntracks, startPosition, phiRange, etaRange,
                       randomNumberSvc, materialValidater,
                       outputMaterialTracks);
  }
}

}  // namespace ActsPython
