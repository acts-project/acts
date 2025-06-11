// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Detector/Detector.hpp"

#include "Acts/Geometry/DetectorElementBase.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Material/IMaterialDecorator.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "ActsExamples/ContextualDetector/AlignedDetector.hpp"
#include "ActsExamples/DetectorCommons/Detector.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsExamples/GenericDetector/GenericDetector.hpp"
#include "ActsExamples/TGeoDetector/TGeoDetector.hpp"
#include "ActsExamples/TelescopeDetector/TelescopeDetector.hpp"
#include "ActsExamples/Utilities/Options.hpp"

#include <memory>
#include <optional>
#include <string>
#include <utility>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/stl/filesystem.h>

namespace py = pybind11;
using namespace ActsExamples;

namespace Acts::Python {

void addDetector(Context& ctx) {
  auto [m, mex] = ctx.get("main", "examples");

  {
    py::class_<IContextDecorator, std::shared_ptr<IContextDecorator>>(
        mex, "IContextDecorator")
        .def("decorate", &IContextDecorator::decorate)
        .def("name", &IContextDecorator::name);
  }

  {
    py::class_<Detector, std::shared_ptr<Detector>>(mex, "DetectorBase")
        .def("nominalGeometryContext", &Detector::nominalGeometryContext)
        .def("trackingGeometry", &Detector::trackingGeometry)
        .def("gen2Geometry", &Detector::gen2Geometry)
        .def("contextDecorators", &Detector::contextDecorators)
        .def("__enter__",
             [](const std::shared_ptr<Detector>& self) { return self; })
        .def("__exit__",
             [](std::shared_ptr<Detector>& self,
                const std::optional<py::object>&,
                const std::optional<py::object>&,
                const std::optional<py::object>&) { self.reset(); });
  }

  {
    auto d =
        py::class_<GenericDetector, Detector, std::shared_ptr<GenericDetector>>(
            mex, "GenericDetector")
            .def(py::init<const GenericDetector::Config&>());

    auto c = py::class_<GenericDetector::Config>(d, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, buildLevel, logLevel, surfaceLogLevel, layerLogLevel,
                       volumeLogLevel, buildProto, materialDecorator, gen3,
                       graphvizFile);
  }

  {
    auto d =
        py::class_<TelescopeDetector, Detector,
                   std::shared_ptr<TelescopeDetector>>(mex, "TelescopeDetector")
            .def(py::init<const TelescopeDetector::Config&>());

    auto c =
        py::class_<TelescopeDetector::Config>(d, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, positions, stereos, offsets, bounds, thickness,
                       surfaceType, binValue, materialDecorator, logLevel);
  }

  {
    auto d =
        py::class_<AlignedDetector, Detector, std::shared_ptr<AlignedDetector>>(
            mex, "AlignedDetector")
            .def(py::init<const AlignedDetector::Config&>());

    auto c = py::class_<AlignedDetector::Config, GenericDetector::Config>(
                 d, "Config")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, seed, iovSize, flushSize, doGarbageCollection,
                       sigmaInPlane, sigmaOutPlane, sigmaInRot, sigmaOutRot,
                       firstIovNominal, decoratorLogLevel, mode);

    py::enum_<AlignedDetector::Config::Mode>(c, "Mode")
        .value("Internal", AlignedDetector::Config::Mode::Internal)
        .value("External", AlignedDetector::Config::Mode::External);
  }

  {
    auto d = py::class_<TGeoDetector, Detector, std::shared_ptr<TGeoDetector>>(
                 mex, "TGeoDetector")
                 .def(py::init<const TGeoDetector::Config&>());

    py::class_<Options::Interval>(mex, "Interval")
        .def(py::init<>())
        .def(py::init<std::optional<double>, std::optional<double>>())
        .def_readwrite("lower", &Options::Interval::lower)
        .def_readwrite("upper", &Options::Interval::upper);

    auto c = py::class_<TGeoDetector::Config>(d, "Config").def(py::init<>());

    c.def_property("jsonFile", nullptr,
                   [](TGeoDetector::Config& cfg, const std::string& file) {
                     cfg.readJson(file);
                   });

    py::enum_<TGeoDetector::Config::SubVolume>(c, "SubVolume")
        .value("Negative", TGeoDetector::Config::SubVolume::Negative)
        .value("Central", TGeoDetector::Config::SubVolume::Central)
        .value("Positive", TGeoDetector::Config::SubVolume::Positive);

    py::enum_<Acts::BinningType>(c, "BinningType")
        .value("equidistant", Acts::BinningType::equidistant)
        .value("arbitrary", Acts::BinningType::arbitrary);

    auto volume =
        py::class_<TGeoDetector::Config::Volume>(c, "Volume").def(py::init<>());
    ACTS_PYTHON_STRUCT(
        volume, name, binToleranceR, binTolerancePhi, binToleranceZ,
        cylinderDiscSplit, cylinderNZSegments, cylinderNPhiSegments,
        discNRSegments, discNPhiSegments, itkModuleSplit, barrelMap, discMap,
        splitPatterns, layers, subVolumeName, sensitiveNames, sensitiveAxes,
        rRange, zRange, splitTolR, splitTolZ, binning0, binning1);

    auto regTriplet = [&c](const std::string& name, auto v) {
      using type = decltype(v);
      py::class_<TGeoDetector::Config::LayerTriplet<type>>(c, name.c_str())
          .def(py::init<>())
          .def(py::init<type>())
          .def(py::init<type, type, type>())
          .def_readwrite("negative",
                         &TGeoDetector::Config::LayerTriplet<type>::negative)
          .def_readwrite("central",
                         &TGeoDetector::Config::LayerTriplet<type>::central)
          .def_readwrite("positive",
                         &TGeoDetector::Config::LayerTriplet<type>::positive)
          .def("at", py::overload_cast<TGeoDetector::Config::SubVolume>(
                         &TGeoDetector::Config::LayerTriplet<type>::at));
    };

    regTriplet("LayerTripletBool", true);
    regTriplet("LayerTripletString", std::string{""});
    regTriplet("LayerTripletVectorString", std::vector<std::string>{});
    regTriplet("LayerTripletInterval", Options::Interval{});
    regTriplet("LayerTripletDouble", double{5.5});
    regTriplet("LayerTripletVectorBinning",
               std::vector<std::pair<int, Acts::BinningType>>{});

    ACTS_PYTHON_STRUCT(c, surfaceLogLevel, layerLogLevel, volumeLogLevel,
                       fileName, buildBeamPipe, beamPipeRadius,
                       beamPipeHalflengthZ, beamPipeLayerThickness,
                       beamPipeEnvelopeR, layerEnvelopeR, unitScalor,
                       materialDecorator, volumes);

    patchKwargsConstructor(c);
  }

  {
    py::class_<Acts::DetectorElementBase,
               std::shared_ptr<Acts::DetectorElementBase>>(
        mex, "DetectorElementBase");
  }
}

}  // namespace Acts::Python
