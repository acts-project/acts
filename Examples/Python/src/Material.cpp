// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/Detector.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"
#include "Acts/Material/IMaterialDecorator.hpp"
#include "Acts/Material/SurfaceMaterialMapper.hpp"
#include "Acts/Material/VolumeMaterialMapper.hpp"
#include "Acts/Material/interface/IMaterialMapper.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"
#include "ActsExamples/Io/Root/RootMaterialDecorator.hpp"
#include "ActsExamples/MaterialMapping/MappingMaterialDecorator.hpp"
#include "ActsExamples/MaterialMapping/MaterialMapping.hpp"

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
    py::class_<Acts::IMaterialDecorator,
               std::shared_ptr<Acts::IMaterialDecorator>>(m,
                                                          "IMaterialDecorator");
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

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(voltag);
    ACTS_PYTHON_MEMBER(boutag);
    ACTS_PYTHON_MEMBER(laytag);
    ACTS_PYTHON_MEMBER(apptag);
    ACTS_PYTHON_MEMBER(sentag);
    ACTS_PYTHON_MEMBER(ntag);
    ACTS_PYTHON_MEMBER(vtag);
    ACTS_PYTHON_MEMBER(otag);
    ACTS_PYTHON_MEMBER(mintag);
    ACTS_PYTHON_MEMBER(maxtag);
    ACTS_PYTHON_MEMBER(ttag);
    ACTS_PYTHON_MEMBER(x0tag);
    ACTS_PYTHON_MEMBER(l0tag);
    ACTS_PYTHON_MEMBER(atag);
    ACTS_PYTHON_MEMBER(ztag);
    ACTS_PYTHON_MEMBER(rhotag);
    ACTS_PYTHON_MEMBER(fileName);
    ACTS_PYTHON_STRUCT_END();
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
                   .def_property_readonly("config", &Alg::config);

    auto c = py::class_<Alg::Config>(alg, "Config")
                 .def(py::init<const Acts::GeometryContext&,
                               const Acts::MagneticFieldContext&>());

    ACTS_PYTHON_STRUCT_BEGIN(c, Alg::Config);
    ACTS_PYTHON_MEMBER(collection);
    ACTS_PYTHON_MEMBER(mappedCollection);
    ACTS_PYTHON_MEMBER(unmappedCollection);
    ACTS_PYTHON_MEMBER(materialMappers);
    ACTS_PYTHON_MEMBER(materialWriters);
    ACTS_PYTHON_MEMBER(trackingGeometry);
    ACTS_PYTHON_MEMBER(detector);
    ACTS_PYTHON_MEMBER(geoContext);
    ACTS_PYTHON_MEMBER(magFieldContext);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    py::class_<IMaterialMapper, std::shared_ptr<IMaterialMapper>>(
        m, "IMaterialMapper");
  }

  {
    auto sMapper =
        py::class_<SurfaceMaterialMapper, IMaterialMapper,
                   std::shared_ptr<SurfaceMaterialMapper>>(
            m, "SurfaceMaterialMapper")
            .def(py::init([](const SurfaceMaterialMapper::Config& config,
                             SurfaceMaterialMapper::StraightLineTGPropagator& prop,
                             Acts::Logging::Level level) {
                   return std::make_shared<SurfaceMaterialMapper>(
                       config, 
                       prop,
                       getDefaultLogger("SurfaceMaterialMapper", level));
                 }),
                 py::arg("config"), py::arg("propagator"), py::arg("level"))
            .def(py::init([](const SurfaceMaterialMapper::Config& config,
                             SurfaceMaterialMapper::StraightLineDetPropagator& prop,
                             Acts::Logging::Level level) {
                   return std::make_shared<SurfaceMaterialMapper>(
                       config, 
                       prop,
                       getDefaultLogger("SurfaceMaterialMapper", level));
                 }),
                 py::arg("config"), py::arg("propagator"), py::arg("level"));

    auto c = py::class_<SurfaceMaterialMapper::Config>(sMapper, "Config")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, SurfaceMaterialMapper::Config);
    ACTS_PYTHON_MEMBER(emptyBinCorrection);
    ACTS_PYTHON_MEMBER(mapperDebugOutput);
    ACTS_PYTHON_MEMBER(computeVariance);
    ACTS_PYTHON_MEMBER(veto);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    auto vMapper =
        py::class_<VolumeMaterialMapper, IMaterialMapper,
                   std::shared_ptr<VolumeMaterialMapper>>(
            m, "VolumeMaterialMapper")
            .def(py::init([](const VolumeMaterialMapper::Config& config,
                             VolumeMaterialMapper::StraightLineTGPropagator& prop,
                             Acts::Logging::Level level) {
                   return std::make_shared<VolumeMaterialMapper>(
                       config, 
                       prop,
                       getDefaultLogger("VolumeMaterialMapper", level));
                 }),
                 py::arg("config"), py::arg("propagator"), py::arg("level"))
            .def(py::init([](const VolumeMaterialMapper::Config& config,
                             VolumeMaterialMapper::StraightLineDetPropagator& prop,
                             Acts::Logging::Level level) {
                   return std::make_shared<VolumeMaterialMapper>(
                       config, 
                       prop,
                       getDefaultLogger("VolumeMaterialMapper", level));
                 }),
                 py::arg("config"), py::arg("propagator"), py::arg("level"));

    auto c = py::class_<VolumeMaterialMapper::Config>(vMapper, "Config")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, VolumeMaterialMapper::Config);
    ACTS_PYTHON_MEMBER(mappingStep);
    ACTS_PYTHON_MEMBER(veto);
    ACTS_PYTHON_STRUCT_END();
  }
}
}  // namespace Acts::Python
