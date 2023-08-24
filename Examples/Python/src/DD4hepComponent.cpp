// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/DD4hepDetectorElement.hpp"
#include "Acts/Plugins/DD4hep/DD4hepDetectorSurfaceFactory.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepGeometryService.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsExamples/Framework/ProcessCode.hpp"

#include <array>
#include <memory>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace Acts {
class IMaterialDecorator;
}  // namespace Acts

namespace py = pybind11;
using namespace ActsExamples;
using namespace Acts::Python;

PYBIND11_MODULE(ActsPythonBindingsDD4hep, m) {
  {
    using Config = ActsExamples::DD4hep::DD4hepGeometryService::Config;
    auto s = py::class_<DD4hep::DD4hepGeometryService,
                        std::shared_ptr<DD4hep::DD4hepGeometryService>>(
                 m, "DD4hepGeometryService")
                 .def(py::init<const Config&>());

    auto c = py::class_<Config>(s, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(logLevel);
    ACTS_PYTHON_MEMBER(dd4hepLogLevel);
    ACTS_PYTHON_MEMBER(xmlFileNames);
    ACTS_PYTHON_MEMBER(name);
    ACTS_PYTHON_MEMBER(bTypePhi);
    ACTS_PYTHON_MEMBER(bTypeR);
    ACTS_PYTHON_MEMBER(bTypeZ);
    ACTS_PYTHON_MEMBER(envelopeR);
    ACTS_PYTHON_MEMBER(envelopeZ);
    ACTS_PYTHON_MEMBER(defaultLayerThickness);
    ACTS_PYTHON_MEMBER(geometryIdentifierHook);
    ACTS_PYTHON_STRUCT_END();

    patchKwargsConstructor(c);
  }

  {
    py::class_<DD4hep::DD4hepDetector, std::shared_ptr<DD4hep::DD4hepDetector>>(
        m, "DD4hepDetector")
        .def(py::init<>())
        .def("finalize",
             py::overload_cast<DD4hep::DD4hepGeometryService::Config,
                               std::shared_ptr<const Acts::IMaterialDecorator>>(
                 &DD4hep::DD4hepDetector::finalize));
  }

  {
    py::class_<Acts::DD4hepDetectorElement,
               std::shared_ptr<Acts::DD4hepDetectorElement>>(
        m, "DD4hepDetectorElement");
  }

  {
    /// This method allows to test conversion of DD4hep surfaces to Acts
    /// surfaces as a first check for detector conversion
    m.def("convertSurfaces", [](const std::vector<std::string>& xmlFileNames,
                                Acts::Logging::Level logLevel) {
      ActsExamples::DD4hep::DD4hepGeometryService::Config cfg;
      cfg.dd4hepLogLevel = logLevel;
      cfg.xmlFileNames = xmlFileNames;
      cfg.name = "DD4hepConvertSurfacesTestService";

      ActsExamples::DD4hep::DD4hepGeometryService dd4HepService(cfg);
      auto dd4HepGeometry = dd4HepService.geometry();

      Acts::DD4hepDetectorSurfaceFactory::Cache cache;
      Acts::DD4hepDetectorSurfaceFactory dd4HepSurfaceFactory(
          Acts::getDefaultLogger("DD4hepDetectorSurfaceFactory", logLevel));
      dd4HepSurfaceFactory.construct(cache, dd4HepGeometry);

      // Capture the sensitive elements and the surfaces
      using Elements =
          std::vector<std::shared_ptr<Acts::DD4hepDetectorElement>>;
      Elements detectorElements;
      detectorElements.reserve(cache.sensitiveSurfaces.size());
      using Surfaces = std::vector<std::shared_ptr<Acts::Surface>>;
      Surfaces surfaces;
      surfaces.reserve(cache.sensitiveSurfaces.size());
      std::for_each(cache.sensitiveSurfaces.begin(),
                    cache.sensitiveSurfaces.end(), [&](const auto& sensitive) {
                      detectorElements.push_back(std::get<0>(sensitive));
                      surfaces.push_back(std::get<1>(sensitive));
                    });

      // Capture the passive surfaces
      Surfaces passiveSurfaces;
      passiveSurfaces.reserve(cache.passiveSurfaces.size());
      for (const auto& passive : cache.passiveSurfaces) {
        passiveSurfaces.push_back(passive);
      }

      // Return a convenient tuple for drawing
      return std::tuple<Elements, Surfaces, Surfaces>(
          std::move(detectorElements), std::move(surfaces),
          std::move(passiveSurfaces));
    });
  }
}
