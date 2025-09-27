// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Detector/GeometryIdGenerator.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/DD4hepDetector/AlignedDD4hepDetectorElement.hpp"
#include "ActsExamples/DD4hepDetector/DD4hepDetector.hpp"
#include "ActsPlugins/DD4hep/DD4hepDetectorElement.hpp"
#include "ActsPlugins/DD4hep/DD4hepDetectorStructure.hpp"
#include "ActsPlugins/DD4hep/DD4hepFieldAdapter.hpp"
#include "ActsPlugins/DD4hep/DD4hepIdentifierMapper.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <memory>
#include <utility>

#include <DD4hep/DetElement.h>
#include <DD4hep/Fields.h>
#include <pybind11/functional.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace Acts;
using namespace ActsPlugins;
using namespace ActsExamples;
using namespace ActsPython;

PYBIND11_MODULE(ActsPythonBindingsDD4hep, m) {
  {
    py::class_<DD4hepDetectorElement, DetectorElementBase,
               std::shared_ptr<DD4hepDetectorElement>>(m,
                                                       "DD4hepDetectorElement");

    py::class_<AlignedDD4hepDetectorElement, DD4hepDetectorElement,
               std::shared_ptr<AlignedDD4hepDetectorElement>>(
        m, "AlignedDD4hepDetectorElement");
  }

  {
    py::class_<dd4hep::DetElement, std::shared_ptr<dd4hep::DetElement>>(
        m, "DD4hepDetElement");
  }

  {
    auto f =
        py::class_<DD4hepDetector, Detector, std::shared_ptr<DD4hepDetector>>(
            m, "DD4hepDetector")
            .def(py::init<const DD4hepDetector::Config&>())
            .def_property_readonly("field", &DD4hepDetector::field);

    auto c = py::class_<DD4hepDetector::Config>(f, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, logLevel, dd4hepLogLevel, xmlFileNames, name,
                       bTypePhi, bTypeR, bTypeZ, envelopeR, envelopeZ,
                       defaultLayerThickness, materialDecorator,
                       geometryIdentifierHook, detectorElementFactory);
    patchKwargsConstructor(c);

    m.def("alignedDD4hepDetectorElementFactory",
          &alignedDD4hepDetectorElementFactory);
  }

  {
    py::class_<DD4hepFieldAdapter, MagneticFieldProvider,
               std::shared_ptr<DD4hepFieldAdapter>>(m, "DD4hepFieldAdapter");
  }

  {
    m.def(
        "createDD4hepIdGeoIdMap",
        [](const TrackingGeometry& tGeometry)
            -> std::map<DD4hepDetectorElement::DD4hepVolumeID,
                        GeometryIdentifier> {
          // The surface visitor
          struct DD4hepIdGrabber {
            std::map<DD4hepDetectorElement::DD4hepVolumeID, GeometryIdentifier>
                dd4hepIdGeoIdMap;

            void operator()(const Surface* surface) {
              const auto* dde = surface->associatedDetectorElement();
              const auto* dd4hepDetElement =
                  dynamic_cast<const DD4hepDetectorElement*>(dde);
              // Check if it is valid
              if (dd4hepDetElement != nullptr) {
                dd4hep::DDSegmentation::VolumeID dd4hepID =
                    dd4hepDetElement->sourceElement().volumeID();
                auto geoID = surface->geometryId();
                dd4hepIdGeoIdMap[dd4hepID] = geoID;
              }
            }
          };

          // Create an instance
          DD4hepIdGrabber dd4hepIdGrabber;
          // Visit the surfaces & return what you have
          tGeometry.visitSurfaces(dd4hepIdGrabber);
          return dd4hepIdGrabber.dd4hepIdGeoIdMap;
        });
  }

  {
    using Options = DD4hepDetectorStructure::Options;
    auto o = py::class_<Options>(m, "DD4hepDetectorOptions").def(py::init<>());
    ACTS_PYTHON_STRUCT(o, logLevel, emulateToGraph, geoIdGenerator,
                       materialDecorator);

    patchKwargsConstructor(o);

    m.def("attachDD4hepGeoIdMapper",
          [](DD4hepDetectorStructure::Options& options,
             const std::map<DD4hepDetectorElement::DD4hepVolumeID,
                            GeometryIdentifier>& dd4hepIdGeoIdMap) {
            // The Geo mapper
            auto geoIdMapper = std::make_shared<const DD4hepIdentifierMapper>(
                DD4hepIdentifierMapper::Config{dd4hepIdGeoIdMap},
                getDefaultLogger("GeometryIdMapper", options.logLevel));

            // A remaining recursive logger
            auto geoIdGenerator =
                std::make_shared<const Experimental::GeometryIdGenerator>(
                    Experimental::GeometryIdGenerator::Config{},
                    getDefaultLogger("GeometryIdGenerator", options.logLevel));

            std::tuple<std::shared_ptr<const Experimental::GeometryIdGenerator>,
                       std::shared_ptr<const DD4hepIdentifierMapper>>
                chainedGenerators = {geoIdGenerator, geoIdMapper};

            auto chainedGeoIdGenerator =
                std::make_shared<const Experimental::ChainedGeometryIdGenerator<
                    std::shared_ptr<const Experimental::GeometryIdGenerator>,
                    std::shared_ptr<const DD4hepIdentifierMapper>>>(
                    std::move(chainedGenerators),
                    getDefaultLogger("ChainedGeometryIdGenerator",
                                     options.logLevel));

            options.geoIdGenerator = chainedGeoIdGenerator;
          });
  }
}
