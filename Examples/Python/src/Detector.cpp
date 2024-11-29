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
#include "ActsExamples/DetectorCommons/DetectorBase.hpp"
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
    py::class_<DetectorBase, std::shared_ptr<DetectorBase>>(mex, "DetectorBase")
        .def("geometryContext", &DetectorBase::geometryContext)
        .def("gen1Geometry", &DetectorBase::gen1Geometry)
        .def("gen2Geometry", &DetectorBase::gen2Geometry)
        .def("contextDecorators", &DetectorBase::contextDecorators)
        .def("__enter__",
             [](const std::shared_ptr<DetectorBase>& self) { return self; })
        .def("__exit__",
             [](std::shared_ptr<DetectorBase>& self,
                const std::optional<py::object>&,
                const std::optional<py::object>&,
                const std::optional<py::object>&) { self.reset(); });

    py::class_<DetectorFactoryBase, std::shared_ptr<DetectorFactoryBase>>(
        mex, "DetectorFactoryBase")
        .def("buildDetector", &DetectorFactoryBase::buildDetector);
  }

  {
    using DetectorFactory = GenericDetectorFactory;
    using Config = DetectorFactory::Config;

    auto d = py::class_<DetectorFactory, DetectorFactoryBase,
                        std::shared_ptr<DetectorFactory>>(
                 mex, "GenericDetectorFactory")
                 .def(py::init<const Config&>());

    auto c = py::class_<Config>(d, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(buildLevel);
    ACTS_PYTHON_MEMBER(logLevel);
    ACTS_PYTHON_MEMBER(surfaceLogLevel);
    ACTS_PYTHON_MEMBER(layerLogLevel);
    ACTS_PYTHON_MEMBER(volumeLogLevel);
    ACTS_PYTHON_MEMBER(buildProto);
    ACTS_PYTHON_MEMBER(materialDecorator);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using DetectorFactory = TelescopeDetectorFactory;
    using Config = DetectorFactory::Config;

    auto d = py::class_<DetectorFactory, DetectorFactoryBase,
                        std::shared_ptr<DetectorFactory>>(
                 mex, "TelescopeDetectorFactory")
                 .def(py::init<const Config&>());

    auto c = py::class_<Config>(d, "Config").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(positions);
    ACTS_PYTHON_MEMBER(stereos);
    ACTS_PYTHON_MEMBER(offsets);
    ACTS_PYTHON_MEMBER(bounds);
    ACTS_PYTHON_MEMBER(thickness);
    ACTS_PYTHON_MEMBER(surfaceType);
    ACTS_PYTHON_MEMBER(binValue);
    ACTS_PYTHON_MEMBER(materialDecorator);
    ACTS_PYTHON_MEMBER(logLevel);
    ACTS_PYTHON_STRUCT_END();
  }

  {
    using DetectorFactory = AlignedDetectorFactory;
    using Config = DetectorFactory::Config;

    auto d = py::class_<DetectorFactory, DetectorFactoryBase,
                        std::shared_ptr<DetectorFactory>>(
                 mex, "AlignedDetectorFactory")
                 .def(py::init<const Config&>());

    auto c = py::class_<Config, GenericDetectorFactory::Config>(d, "Config")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(seed);
    ACTS_PYTHON_MEMBER(iovSize);
    ACTS_PYTHON_MEMBER(flushSize);
    ACTS_PYTHON_MEMBER(doGarbageCollection);
    ACTS_PYTHON_MEMBER(sigmaInPlane);
    ACTS_PYTHON_MEMBER(sigmaOutPlane);
    ACTS_PYTHON_MEMBER(sigmaInRot);
    ACTS_PYTHON_MEMBER(sigmaOutRot);
    ACTS_PYTHON_MEMBER(firstIovNominal);
    ACTS_PYTHON_MEMBER(decoratorLogLevel);
    ACTS_PYTHON_MEMBER(mode);
    ACTS_PYTHON_STRUCT_END();

    py::enum_<Config::Mode>(c, "Mode")
        .value("Internal", Config::Mode::Internal)
        .value("External", Config::Mode::External);
  }

  {
    using DetectorFactory = TGeoDetectorFactory;
    using Config = DetectorFactory::Config;

    auto d =
        py::class_<DetectorFactory, DetectorFactoryBase,
                   std::shared_ptr<DetectorFactory>>(mex, "TGeoDetectorFactory")
            .def(py::init<const Config&>());

    py::class_<Options::Interval>(mex, "Interval")
        .def(py::init<>())
        .def(py::init<std::optional<double>, std::optional<double>>())
        .def_readwrite("lower", &Options::Interval::lower)
        .def_readwrite("upper", &Options::Interval::upper);

    auto c = py::class_<Config>(d, "Config").def(py::init<>());

    c.def_property(
        "jsonFile", nullptr,
        [](Config& cfg, const std::string& file) { cfg.readJson(file); });

    py::enum_<Config::SubVolume>(c, "SubVolume")
        .value("Negative", Config::SubVolume::Negative)
        .value("Central", Config::SubVolume::Central)
        .value("Positive", Config::SubVolume::Positive);

    py::enum_<Acts::BinningType>(c, "BinningType")
        .value("equidistant", Acts::BinningType::equidistant)
        .value("arbitrary", Acts::BinningType::arbitrary);

    auto volume = py::class_<Config::Volume>(c, "Volume").def(py::init<>());
    ACTS_PYTHON_STRUCT_BEGIN(volume, Config::Volume);
    ACTS_PYTHON_MEMBER(name);
    ACTS_PYTHON_MEMBER(binToleranceR);
    ACTS_PYTHON_MEMBER(binTolerancePhi);
    ACTS_PYTHON_MEMBER(binToleranceZ);
    ACTS_PYTHON_MEMBER(cylinderDiscSplit);
    ACTS_PYTHON_MEMBER(cylinderNZSegments);
    ACTS_PYTHON_MEMBER(cylinderNPhiSegments);
    ACTS_PYTHON_MEMBER(discNRSegments);
    ACTS_PYTHON_MEMBER(discNPhiSegments);
    ACTS_PYTHON_MEMBER(itkModuleSplit);
    ACTS_PYTHON_MEMBER(barrelMap);
    ACTS_PYTHON_MEMBER(discMap);
    ACTS_PYTHON_MEMBER(splitPatterns);

    ACTS_PYTHON_MEMBER(layers);
    ACTS_PYTHON_MEMBER(subVolumeName);
    ACTS_PYTHON_MEMBER(sensitiveNames);
    ACTS_PYTHON_MEMBER(sensitiveAxes);
    ACTS_PYTHON_MEMBER(rRange);
    ACTS_PYTHON_MEMBER(zRange);
    ACTS_PYTHON_MEMBER(splitTolR);
    ACTS_PYTHON_MEMBER(splitTolZ);
    ACTS_PYTHON_MEMBER(binning0);
    ACTS_PYTHON_MEMBER(binning1);
    ACTS_PYTHON_STRUCT_END();

    auto regTriplet = [&c](const std::string& name, auto v) {
      using type = decltype(v);
      py::class_<Config::LayerTriplet<type>>(c, name.c_str())
          .def(py::init<>())
          .def(py::init<type>())
          .def(py::init<type, type, type>())
          .def_readwrite("negative", &Config::LayerTriplet<type>::negative)
          .def_readwrite("central", &Config::LayerTriplet<type>::central)
          .def_readwrite("positive", &Config::LayerTriplet<type>::positive)
          .def("at", py::overload_cast<Config::SubVolume>(
                         &Config::LayerTriplet<type>::at));
    };

    regTriplet("LayerTripletBool", true);
    regTriplet("LayerTripletString", std::string{""});
    regTriplet("LayerTripletVectorString", std::vector<std::string>{});
    regTriplet("LayerTripletInterval", Options::Interval{});
    regTriplet("LayerTripletDouble", double{5.5});
    regTriplet("LayerTripletVectorBinning",
               std::vector<std::pair<int, Acts::BinningType>>{});

    ACTS_PYTHON_STRUCT_BEGIN(c, Config);
    ACTS_PYTHON_MEMBER(surfaceLogLevel);
    ACTS_PYTHON_MEMBER(layerLogLevel);
    ACTS_PYTHON_MEMBER(volumeLogLevel);
    ACTS_PYTHON_MEMBER(fileName);
    ACTS_PYTHON_MEMBER(buildBeamPipe);
    ACTS_PYTHON_MEMBER(beamPipeRadius);
    ACTS_PYTHON_MEMBER(beamPipeHalflengthZ);
    ACTS_PYTHON_MEMBER(beamPipeLayerThickness);
    ACTS_PYTHON_MEMBER(beamPipeEnvelopeR);
    ACTS_PYTHON_MEMBER(layerEnvelopeR);
    ACTS_PYTHON_MEMBER(unitScalor);
    ACTS_PYTHON_MEMBER(volumes);
    ACTS_PYTHON_STRUCT_END();

    patchKwargsConstructor(c);
  }

  {
    py::class_<Acts::DetectorElementBase,
               std::shared_ptr<Acts::DetectorElementBase>>(
        mex, "DetectorElementBase");
  }
}

}  // namespace Acts::Python
