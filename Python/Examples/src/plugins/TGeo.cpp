// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/TGeoDetector/TGeoDetector.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsExamples;
using namespace ActsPython;

PYBIND11_MODULE(ActsExamplesPythonBindingsTGeo, tgeo) {
  {
    auto d = py::class_<TGeoDetector, Detector, std::shared_ptr<TGeoDetector>>(
                 tgeo, "TGeoDetector")
                 .def(py::init<const TGeoDetector::Config&>());

    py::class_<Options::Interval>(tgeo, "Interval")
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

    py::enum_<BinningType>(c, "BinningType")
        .value("equidistant", BinningType::equidistant)
        .value("arbitrary", BinningType::arbitrary);

    auto volume =
        py::class_<TGeoDetector::Config::Volume>(c, "Volume").def(py::init<>());
    ACTS_PYTHON_STRUCT(
        volume, name, binToleranceR, binTolerancePhi, binToleranceZ,
        cylinderDiscSplit, cylinderNZSegments, cylinderNPhiSegments,
        discNRSegments, discNPhiSegments, itkModuleSplit, barrelMap, discMap,
        splitPatterns, layers, subVolumeName, sensitiveNames, sensitiveAxes,
        rRange, zRange, splitTolR, splitTolZ, binning0, binning1);

    auto regTriplet = [&c](const std::string& name, const auto& v) {
      using type = std::decay_t<decltype(v)>;
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
               std::vector<std::pair<int, BinningType>>{});

    ACTS_PYTHON_STRUCT(c, surfaceLogLevel, layerLogLevel, volumeLogLevel,
                       fileName, buildBeamPipe, beamPipeRadius,
                       beamPipeHalflengthZ, beamPipeLayerThickness,
                       beamPipeEnvelopeR, layerEnvelopeR, unitScalor,
                       materialDecorator, volumes);

    patchKwargsConstructor(c);
  }
}
