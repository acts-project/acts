// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/CylinderVolumeHelper.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/LayerArrayCreator.hpp"
#include "Acts/Geometry/LayerCreator.hpp"
#include "Acts/Geometry/SurfaceArrayCreator.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"

#include <memory>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

namespace Acts::Python {
void addGeometryBuildingGen1(Context &ctx) {
  auto m = ctx.get("main");

  using SurfacePtrVector = std::vector<std::shared_ptr<const Surface>>;

  py::class_<Acts::Layer, std::shared_ptr<Acts::Layer>>(m, "Layer");

  {
    auto creator =
        py::class_<Acts::LayerCreator>(m, "LayerCreator")
            .def(py::init([](const Acts::LayerCreator::Config &cfg,
                             Acts::Logging::Level level) {
              return Acts::LayerCreator(
                  cfg, Acts::getDefaultLogger("LayerCreator", level));
            }))
            .def(
                "cylinderLayer",
                [](const Acts::LayerCreator &self, const GeometryContext &gctx,
                   SurfacePtrVector surfaces, std::size_t binsPhi,
                   std::size_t binsZ) {
                  return self.cylinderLayer(gctx, std::move(surfaces), binsPhi,
                                            binsZ);
                },
                "gctx"_a, "surfaces"_a, "binsPhi"_a, "binsZ"_a)
            .def(
                "discLayer",
                [](const Acts::LayerCreator &self, const GeometryContext &gctx,
                   SurfacePtrVector surfaces, std::size_t binsR,
                   std::size_t binsPhi) {
                  return self.discLayer(gctx, std::move(surfaces), binsR,
                                        binsPhi);
                },
                "gctx"_a, "surfaces"_a, "binsR"_a, "binsPhi"_a);

    auto config =
        py::class_<LayerCreator::Config>(creator, "Config").def(py::init<>());

    ACTS_PYTHON_STRUCT(config, surfaceArrayCreator, cylinderZtolerance,
                       cylinderPhiTolerance, defaultEnvelopeR,
                       defaultEnvelopeZ);
  }

  {
    using Creator = Acts::SurfaceArrayCreator;
    using Config = typename Creator::Config;

    auto creator =
        py::class_<Creator, std::shared_ptr<Creator>>(m, "SurfaceArrayCreator")
            .def(py::init<Config>());

    py::class_<Config>(creator, "Config").def(py::init<>());
  }

  {
    using Base = Acts::ILayerArrayCreator;
    using Creator = Acts::LayerArrayCreator;
    using Config = typename Creator::Config;

    py::class_<Base, std::shared_ptr<Base>>(m, "ILayerArrayCreator");

    auto creator = py::class_<Creator, std::shared_ptr<Creator>, Base>(
                       m, "LayerArrayCreator")
                       .def(py::init<Config>());

    py::class_<Config>(creator, "Config").def(py::init<>());
  }

  {
    using Base = Acts::ITrackingVolumeArrayCreator;
    using Creator = Acts::TrackingVolumeArrayCreator;
    using Config = typename Creator::Config;

    py::class_<Base, std::shared_ptr<Base>>(m, "ITrackingVolumeArrayCreator");

    auto creator = py::class_<Creator, std::shared_ptr<Creator>, Base>(
                       m, "TrackingVolumeArrayCreator")
                       .def(py::init<Config>());

    py::class_<Config>(creator, "Config").def(py::init<>());
  }

  {
    auto helper =
        py::class_<Acts::CylinderVolumeHelper>(m, "CylinderVolumeHelper")
            .def(py::init([](const Acts::CylinderVolumeHelper::Config &cfg,
                             Acts::Logging::Level level) {
              return Acts::CylinderVolumeHelper(
                  cfg, Acts::getDefaultLogger("CylinderVolumeHelper", level));
            }))
            .def("createTrackingVolume",
                 [](const Acts::CylinderVolumeHelper &self,
                    GeometryContext gctx, const LayerVector &layers,
                    std::shared_ptr<VolumeBounds> volumeBounds,
                    const Transform3 &trafo, const std::string &name) {
                   return self.createTrackingVolume(gctx, layers, {},
                                                    std::move(volumeBounds), {},
                                                    trafo, name);
                 })
            .def("createContainerTrackingVolume",
                 &Acts::CylinderVolumeHelper::createContainerTrackingVolume);

    auto config = py::class_<CylinderVolumeHelper::Config>(helper, "Config")
                      .def(py::init<>());

    ACTS_PYTHON_STRUCT(config, layerArrayCreator, trackingVolumeArrayCreator,
                       passiveLayerThickness, passiveLayerPhiBins,
                       passiveLayerRzBins);
  }
}

}  // namespace Acts::Python
