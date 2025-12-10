// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/MaterialMapping/MaterialMapping.hpp"

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/MaterialMapping/CoreMaterialMapping.hpp"
#include "ActsExamples/MaterialMapping/MappingMaterialDecorator.hpp"
#include "ActsExamples/MaterialMapping/MaterialValidation.hpp"
#include "ActsPlugins/Json/ActsJson.hpp"
#include "ActsPlugins/Json/MaterialMapJsonConverter.hpp"
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
namespace ActsExamples {
class IAlgorithm;
}  // namespace ActsExamples

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsExamples;

namespace ActsPython {
void addMaterialMapping(py::module& mex) {
  {
    using Alg = MaterialMapping;

    auto alg = py::class_<Alg, IAlgorithm, std::shared_ptr<Alg>>(
                   mex, "MaterialMapping")
                   .def(py::init<const Alg::Config&, Logging::Level>(),
                        py::arg("config"), py::arg("level"))
                   .def("scoringParameters", &Alg::scoringParameters)
                   .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Alg::Config>(alg, "Config")
                 .def(py::init<const GeometryContext&,
                               const MagneticFieldContext&>());

    ACTS_PYTHON_STRUCT(c, inputMaterialTracks, mappingMaterialCollection,
                       materialSurfaceMapper, materialVolumeMapper,
                       materialWriters, trackingGeometry, geoContext,
                       magFieldContext);
  }

  {
    py::class_<MappingMaterialDecorator, IMaterialDecorator,
               std::shared_ptr<MappingMaterialDecorator>>(
        mex, "MappingMaterialDecorator")
        .def(py::init<const TrackingGeometry&, Logging::Level, bool, bool>(),
             py::arg("tGeometry"), py::arg("level"),
             py::arg("clearSurfaceMaterial") = true,
             py::arg("clearVolumeMaterial") = true)
        .def("binningMap", &MappingMaterialDecorator::binningMap)
        .def("setBinningMap", &MappingMaterialDecorator::setBinningMap);
  }

  {
    auto mmca =
        py::class_<CoreMaterialMapping, IAlgorithm,
                   std::shared_ptr<CoreMaterialMapping>>(mex,
                                                         "CoreMaterialMapping")
            .def(py::init<const CoreMaterialMapping::Config&, Logging::Level>(),
                 py::arg("config"), py::arg("level"));

    auto c = py::class_<CoreMaterialMapping::Config>(mmca, "Config")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, inputMaterialTracks, mappedMaterialTracks,
                       unmappedMaterialTracks, materialMapper,
                       materiaMaplWriters);
  }

  {
    auto mv =
        py::class_<MaterialValidation, IAlgorithm,
                   std::shared_ptr<MaterialValidation>>(mex,
                                                        "MaterialValidation")
            .def(py::init<const MaterialValidation::Config&, Logging::Level>(),
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
