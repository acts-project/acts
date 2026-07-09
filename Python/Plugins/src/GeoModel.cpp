// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Must be on top to avoid some conflict between forward declare and typedef
// Needed until https://gitlab.cern.ch/GeoModelDev/GeoModel/-/merge_requests/351
// is deployed
// clang-format off
#include <GeoModelRead/ReadGeoModel.h>
// clang-format on

#include "Acts/Geometry/ITrackingGeometryBuilder.hpp"
#include "Acts/Surfaces/AnnulusBounds.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Surfaces/PlaneSurface.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "ActsPlugins/GeoModel/GeoModelConverters.hpp"
#include "ActsPlugins/GeoModel/GeoModelDetectorElement.hpp"
#include "ActsPlugins/GeoModel/GeoModelDetectorElementITk.hpp"
#include "ActsPlugins/GeoModel/GeoModelDetectorObjectFactory.hpp"
#include "ActsPlugins/GeoModel/GeoModelReader.hpp"
#include "ActsPlugins/GeoModel/GeoModelTree.hpp"
#include "ActsPlugins/GeoModel/IGeoShapeConverter.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <string>

#include <GeoModelKernel/GeoFullPhysVol.h>
#include <GeoModelKernel/GeoVPhysVol.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

PYBIND11_MODULE(ActsPluginsPythonBindingsGeoModel, gm) {
  using namespace Acts;
  using namespace ActsPython;
  using namespace ActsPlugins;

  // Basic bindings
  {
    py::class_<GeoModelTree::FpvConstLink>(gm, "GeoModelTree::FpvConstLink")
        .def(py::init<>())
        .def("get", &GeoModelTree::FpvConstLink::get,
             py::return_value_policy::reference);

    py::class_<GeoModelTree>(gm, "GeoModelTree").def(py::init<>());

    gm.def("readFromDb", &GeoModelReader::readFromDb);

    py::class_<GeoModelDetectorElement,
               std::shared_ptr<GeoModelDetectorElement>>(
        gm, "GeoModelDetectorElement")
        .def("logVolName", &GeoModelDetectorElement::logVolName)
        .def("databaseEntryName", &GeoModelDetectorElement::databaseEntryName)
        .def("surface", [](GeoModelDetectorElement self) {
          return self.surface().getSharedPtr();
        });

    py::class_<GeoModelDetectorElementITk,
               std::shared_ptr<GeoModelDetectorElementITk>>(
        gm, "GeoModelDetectorElementITk")
        .def("surface", [](GeoModelDetectorElementITk& self) {
          return self.surface().getSharedPtr();
        });
    gm.def("convertToItk", &GeoModelDetectorElementITk::convertFromGeomodel);
  }

  // Shape converters
  {
    py::class_<IGeoShapeConverter, std::shared_ptr<IGeoShapeConverter>>(
        gm, "IGeoShapeConverter");

    py::class_<GeoBoxConverter, IGeoShapeConverter,
               std::shared_ptr<GeoBoxConverter>>(gm, "GeoBoxConverter")
        .def(py::init<>())
        .def("toSensitiveSurface", &GeoBoxConverter::toSensitiveSurface)
        .def("toPassiveSurface", &GeoBoxConverter::toPassiveSurface);

    py::class_<GeoTrdConverter, IGeoShapeConverter,
               std::shared_ptr<GeoTrdConverter>>(gm, "GeoTrdConverter")
        .def(py::init<>())
        .def("toSensitiveSurface", &GeoTrdConverter::toSensitiveSurface)
        .def("toPassiveSurface", &GeoTrdConverter::toPassiveSurface);

    py::class_<GeoTubeConverter, IGeoShapeConverter,
               std::shared_ptr<GeoTubeConverter>>(gm, "GeoTubeConverter")
        .def(py::init<>())
        .def("toSensitiveSurface", &GeoTubeConverter::toSensitiveSurface)
        .def("toPassiveSurface", &GeoTubeConverter::toPassiveSurface);

    py::class_<GeoUnionDoubleTrdConverter, IGeoShapeConverter,
               std::shared_ptr<GeoUnionDoubleTrdConverter>>(
        gm, "GeoUnionDoubleTrdConverter")
        .def(py::init<>())
        .def("toSensitiveSurface",
             &GeoUnionDoubleTrdConverter::toSensitiveSurface)
        .def("toPassiveSurface", &GeoUnionDoubleTrdConverter::toPassiveSurface);

    py::class_<GeoIntersectionAnnulusConverter, IGeoShapeConverter,
               std::shared_ptr<GeoIntersectionAnnulusConverter>>(
        gm, "GeoIntersectionAnnulusConverter")
        .def(py::init<>())
        .def("toSensitiveSurface",
             &GeoIntersectionAnnulusConverter::toSensitiveSurface)
        .def("toPassiveSurface",
             &GeoIntersectionAnnulusConverter::toPassiveSurface);

    py::class_<GeoShiftConverter, IGeoShapeConverter,
               std::shared_ptr<GeoShiftConverter>>(gm, "GeoShiftConverter")
        .def(py::init<>())
        .def("toSensitiveSurface", &GeoShiftConverter::toSensitiveSurface)
        .def("toPassiveSurface", &GeoShiftConverter::toPassiveSurface);
  }

  // Volume factory
  {
    auto a =
        py::class_<GeoModelDetectorObjectFactory,
                   std::shared_ptr<GeoModelDetectorObjectFactory>>(
            gm, "GeoModelDetectorObjectFactory")
            .def(py::init([](const GeoModelDetectorObjectFactory::Config& cfg,
                             Logging::Level level) {
              return std::make_shared<GeoModelDetectorObjectFactory>(
                  cfg,
                  getDefaultLogger("GeoModelDetectorObjectFactory", level));
            }))
            .def("construct", &GeoModelDetectorObjectFactory::construct);

    py::class_<GeoModelDetectorObjectFactory::Config>(a, "Config")
        .def(py::init<>())
        .def_readwrite(
            "convertSubVolumes",
            &GeoModelDetectorObjectFactory::Config::convertSubVolumes)
        .def_readwrite("nameList",
                       &GeoModelDetectorObjectFactory::Config::nameList)
        .def_readwrite("convertBox",
                       &GeoModelDetectorObjectFactory::Config::convertBox)
        .def_readwrite("materialList",
                       &GeoModelDetectorObjectFactory::Config::materialList);

    auto convVol = py::class_<GeoModelDetectorObjectFactory::ConvertedGeoVol>(
        a, "ConvertedGeoVol");

    ACTS_PYTHON_STRUCT(convVol, volume, fullPhysVol, name, surfaces);
    py::class_<GeoModelDetectorObjectFactory::Cache>(a, "Cache")
        .def(py::init<>())
        .def_readwrite("sensitiveSurfaces",
                       &GeoModelDetectorObjectFactory::Cache::sensitiveSurfaces)
        .def_readwrite("boundingBoxes",
                       &GeoModelDetectorObjectFactory::Cache::volumeBoxFPVs);

    py::class_<GeoModelDetectorObjectFactory::Options>(a, "Options")
        .def(py::init<>())
        .def_readwrite("queries",
                       &GeoModelDetectorObjectFactory::Options::queries);
  }
}
