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
void addMaterialLegacy(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");



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
