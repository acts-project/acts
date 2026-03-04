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
    auto [alg, c] =
        declareAlgorithm<MaterialMapping, IAlgorithm>(mex, "MaterialMapping");
    alg.def("scoringParameters", &MaterialMapping::scoringParameters);
    // Config also needs the geometry+field constructor because geoContext and
    // magFieldContext are reference_wrapper fields with no default constructor.
    c.def(py::init([](const Acts::GeometryContext& gc,
                      const Acts::MagneticFieldContext& mfc) {
            MaterialMapping::Config cfg{gc, mfc};
            return cfg;
          }),
          py::arg("geoContext"), py::arg("magFieldContext"));
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
    auto [mmca, c] = declareAlgorithm<CoreMaterialMapping, IAlgorithm>(
        mex, "CoreMaterialMapping");
    ACTS_PYTHON_STRUCT(c, inputMaterialTracks, mappedMaterialTracks,
                       unmappedMaterialTracks, materialMapper,
                       materiaMaplWriters);
  }

  {
    auto [mv, c] = declareAlgorithm<MaterialValidation, IAlgorithm>(
        mex, "MaterialValidation");
    mv.def("execute", &MaterialValidation::execute);
    ACTS_PYTHON_STRUCT(c, ntracks, startPosition, phiRange, etaRange,
                       randomNumberSvc, materialValidater,
                       outputMaterialTracks);
  }
}

}  // namespace ActsPython
