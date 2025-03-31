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
#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Io/Root/RootMaterialDecorator.hpp"
#include "ActsExamples/MaterialMapping/CoreMaterialMapping.hpp"
#include "ActsExamples/MaterialMapping/MappingMaterialDecorator.hpp"
#include "ActsExamples/MaterialMapping/MaterialMapping.hpp"
#include "ActsExamples/MaterialMapping/MaterialValidation.hpp"

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

namespace Acts::Python {
void addMaterial(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  {
    py::class_<Acts::ISurfaceMaterial, std::shared_ptr<ISurfaceMaterial>>(
        m, "ISurfaceMaterial")
        .def("toString", &Acts::ISurfaceMaterial::toString);

    py::class_<Acts::ProtoGridSurfaceMaterial, Acts::ISurfaceMaterial,
               std::shared_ptr<ProtoGridSurfaceMaterial>>(
        m, "ProtoGridSurfaceMaterial");

    py::class_<Acts::ProtoSurfaceMaterial, Acts::ISurfaceMaterial,
               std::shared_ptr<ProtoSurfaceMaterial>>(m,
                                                      "ProtoSurfaceMaterial");

    py::class_<Acts::HomogeneousSurfaceMaterial, Acts::ISurfaceMaterial,
               std::shared_ptr<HomogeneousSurfaceMaterial>>(
        m, "HomogeneousSurfaceMaterial");

    py::class_<Acts::IVolumeMaterial, std::shared_ptr<IVolumeMaterial>>(
        m, "IVolumeMaterial");
  }

  {
    py::class_<Acts::IMaterialDecorator,
               std::shared_ptr<Acts::IMaterialDecorator>>(m,
                                                          "IMaterialDecorator")
        .def("decorate", py::overload_cast<Surface&>(
                             &Acts::IMaterialDecorator::decorate, py::const_));
  }

  {
    auto rmd =
        py::class_<RootMaterialDecorator, Acts::IMaterialDecorator,
                   std::shared_ptr<RootMaterialDecorator>>(
            mex, "RootMaterialDecorator")
            .def(
                py::init<RootMaterialDecorator::Config, Acts::Logging::Level>(),
                py::arg("config"), py::arg("level"));

    using Config = RootMaterialDecorator::Config;
    auto c = py::class_<Config>(rmd, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT(c, voltag, boutag, laytag, apptag, sentag, ntag, vtag,
                       otag, mintag, maxtag, ttag, x0tag, l0tag, atag, ztag,
                       rhotag, fileName);
  }

  {
    py::class_<MappingMaterialDecorator, Acts::IMaterialDecorator,
               std::shared_ptr<MappingMaterialDecorator>>(
        m, "MappingMaterialDecorator")
        .def(py::init<const Acts::TrackingGeometry&, Acts::Logging::Level, bool,
                      bool>(),
             py::arg("tGeometry"), py::arg("level"),
             py::arg("clearSurfaceMaterial") = true,
             py::arg("clearVolumeMaterial") = true)
        .def("binningMap", &MappingMaterialDecorator::binningMap)
        .def("setBinningMap", &MappingMaterialDecorator::setBinningMap);
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
        py::class_<SurfaceMaterialMapper,
                   std::shared_ptr<SurfaceMaterialMapper>>(
            m, "SurfaceMaterialMapper")
            .def(py::init([](const SurfaceMaterialMapper::Config& config,
                             SurfaceMaterialMapper::StraightLinePropagator prop,
                             Acts::Logging::Level level) {
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
                             Acts::Logging::Level level) {
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
                       config,
                       getDefaultLogger("IntersectionMaterialAssigner", level));
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
        py::class_<BinnedSurfaceMaterialAccumulater,
                   ISurfaceMaterialAccumulater,
                   std::shared_ptr<BinnedSurfaceMaterialAccumulater>>(
            m, "BinnedSurfaceMaterialAccumulater")
            .def(
                py::init(
                    [](const BinnedSurfaceMaterialAccumulater::Config& config,
                       Acts::Logging::Level level) {
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
                                   Acts::Logging::Level level) {
                         return std::make_shared<MaterialMapper>(
                             config, getDefaultLogger("MaterialMapper", level));
                       }),
                       py::arg("config"), py::arg("level"));

    auto c = py::class_<MaterialMapper::Config>(mm, "Config").def(py::init<>());
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
        py::class_<MaterialValidater, std::shared_ptr<MaterialValidater>>(
            m, "MaterialValidater")
            .def(py::init([](const MaterialValidater::Config& config,
                             Acts::Logging::Level level) {
                   return std::make_shared<MaterialValidater>(
                       config, getDefaultLogger("MaterialValidater", level));
                 }),
                 py::arg("config"), py::arg("level"))
            .def("recordMaterial", &MaterialValidater::recordMaterial);

    auto c =
        py::class_<MaterialValidater::Config>(mvc, "Config").def(py::init<>());
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

}  // namespace Acts::Python
