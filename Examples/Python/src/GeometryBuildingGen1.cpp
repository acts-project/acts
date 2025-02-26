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
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Geometry/TrackingVolumeArrayCreator.hpp"
#include "Acts/Geometry/VolumeBounds.hpp"
#include "Acts/Geometry/detail/PolyhedronIntersection.hpp"
#include "Acts/Plugins/Python/Utilities.hpp"
#include "Acts/Surfaces/SurfaceBounds.hpp"

#include <algorithm>
#include <memory>
#include <vector>

#include <boost/container/flat_set.hpp>
#include <boost/geometry.hpp>
#include <boost/timer/progress_display.hpp>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

namespace Acts::detail {

namespace bg = boost::geometry;

std::vector<std::array<GeometryIdentifier, 2>> checkSensitiveOverlaps(
    std::shared_ptr<const TrackingGeometry> tg) {
  std::vector<std::pair<Polyhedron, GeometryIdentifier>> polyhedrons;
  tg->visitSurfaces(
      [&](const Surface *s) {
        auto nPoints = 0;
        if (s->type() == Surface::Plane &&
            s->bounds().type() == SurfaceBounds::BoundsType::eRectangle) {
          nPoints = 4;
        }
        auto polyhedron = s->polyhedronRepresentation({}, nPoints);

        polyhedrons.emplace_back(polyhedron, s->geometryId());
      },
      true);

  std::vector<std::array<GeometryIdentifier, 2>> overlaps;

  auto nComparisons = polyhedrons.size() * (polyhedrons.size() - 1) / 2;
  boost::timer::progress_display show_progress(nComparisons);

  auto &polygons = polyhedrons;
  std::vector<Vector3> set;
  std::size_t skipShared = 0;
  for (std::size_t i = 0; i < polygons.size(); ++i) {
    for (std::size_t j = i + 1; j < polygons.size(); ++j) {
      const auto &p1 = polygons[i].first;
      const auto &p2 = polygons[j].first;
      set.clear();
      std::ranges::copy(p1.vertices, std::back_inserter(set));
      std::ranges::copy(p2.vertices, std::back_inserter(set));
      std::sort(set.begin(), set.end(), [](const auto &a, const auto &b) {
        return a.x() < b.x() || (a.x() == b.x() && a.y() < b.y()) ||
               (a.x() == b.x() && a.y() == b.y() && a.z() < b.z());
      });
      set.erase(std::unique(set.begin(), set.end(),
                            [](const auto &a, const auto &b) {
                              return a.isApprox(b, 1e-4);
                            }),
                set.end());
      if (set.size() < p1.vertices.size() + p2.vertices.size()) {
        ++skipShared;
      } else if (polyhedronIntersection(p2, p1, 1e-4)) {
        overlaps.push_back({polygons[i].second, polygons[j].second});
      }
      ++show_progress;
    }
  }

  std::cout << "skipped because of shared vertices: " << skipShared
            << std::endl;
  std::cout << "found overlaps: " << overlaps.size() << std::endl;
  std::cout << "total checks: " << nComparisons << std::endl;

  std::set<GeometryIdentifier> uniqueIds;
  for (const auto &[a, b] : overlaps) {
    uniqueIds.insert(a);
    uniqueIds.insert(b);
  }

  std::cout << "involved surfaces: " << uniqueIds.size() << " / "
            << polyhedrons.size() << std::endl;

  return overlaps;
}

}  // namespace Acts::detail

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

    ACTS_PYTHON_STRUCT_BEGIN(config, LayerCreator::Config);
    ACTS_PYTHON_MEMBER(surfaceArrayCreator);
    ACTS_PYTHON_MEMBER(cylinderZtolerance);
    ACTS_PYTHON_MEMBER(cylinderPhiTolerance);
    ACTS_PYTHON_MEMBER(defaultEnvelopeR);
    ACTS_PYTHON_MEMBER(defaultEnvelopeZ);
    ACTS_PYTHON_STRUCT_END();
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

    ACTS_PYTHON_STRUCT_BEGIN(config, CylinderVolumeHelper::Config);
    ACTS_PYTHON_MEMBER(layerArrayCreator);
    ACTS_PYTHON_MEMBER(trackingVolumeArrayCreator);
    ACTS_PYTHON_MEMBER(passiveLayerThickness);
    ACTS_PYTHON_MEMBER(passiveLayerPhiBins);
    ACTS_PYTHON_MEMBER(passiveLayerRzBins);
    ACTS_PYTHON_STRUCT_END();
  }

  { m.def("checkSensitiveOverlaps", &Acts::detail::checkSensitiveOverlaps); }
}

}  // namespace Acts::Python
