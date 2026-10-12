// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Nodd/DD4hepNoddDetector.hpp"

#include "Acts/Definitions/Units.hpp"
#include "Acts/Geometry/Blueprint.hpp"
#include "Acts/Geometry/BlueprintOptions.hpp"
#include "Acts/Geometry/ContainerBlueprintNode.hpp"
#include "Acts/Geometry/NavigationPolicyFactory.hpp"
#include "Acts/Geometry/VolumeAttachmentStrategy.hpp"
#include "Acts/Geometry/VolumeResizeStrategy.hpp"
#include "Acts/Material/HomogeneousSurfaceMaterial.hpp"
#include "Acts/Material/Material.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Navigation/TryAllNavigationPolicy.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/RectangleBounds.hpp"
#include "ActsPlugins/DD4hep/BlueprintBuilder.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <map>
#include <set>
#include <stdexcept>
#include <utility>

#include <DD4hep/Alignments.h>
#include <TGeoMaterial.h>
#include <TGeoMedium.h>
#include <TGeoTube.h>
#include <nlohmann/json.hpp>

namespace ActsExamples {
namespace {
using Element = dd4hep::DetElement;
using Ids = std::map<std::string, int>;

Ids placementIds(Element element) {
  Ids ids;
  for (; element.isValid(); element = element.parent()) {
    if (element.placement().isValid()) {
      for (const auto& [key, value] : element.placement().volIDs()) {
        if (!ids.emplace(key, value).second && ids.at(key) != value) {
          throw std::runtime_error("Conflicting DD4hep placement ID: " + key);
        }
      }
    }
  }
  return ids;
}

std::shared_ptr<const Acts::ISurfaceMaterial> nativeMaterial(
    const TGeoMaterial& material, double thickness) {
  using namespace Acts::UnitConstants;
  return std::make_shared<Acts::HomogeneousSurfaceMaterial>(Acts::MaterialSlab(
      Acts::Material::fromMassDensity(
          material.GetRadLen() * cm, material.GetIntLen() * cm, material.GetA(),
          material.GetZ(), material.GetDensity() * g / (cm * cm * cm)),
      thickness));
}

std::unique_ptr<Acts::TrackingGeometry> buildGeometry(
    const dd4hep::Detector& detector, const Acts::GeometryContext& gctx,
    const Acts::Logger& logger, double envelope) {
  using namespace Acts;
  using enum AxisDirection;
  ActsPlugins::DD4hep::BlueprintBuilder builder{
      {.dd4hepDetector = &detector,
       .lengthScale = UnitConstants::cm,
       .gctx = gctx},
      logger.cloneWithSuffix("Gen3")};
  std::map<std::pair<int, int>, std::vector<Element>> groups;
  std::function<void(Element)> visit = [&](Element element) {
    if (element.volume().isSensitive()) {
      const auto ids = placementIds(element);
      const auto system = ids.at("system");
      if (system < 1 || system > 7) {
        throw std::runtime_error("Unexpected nODD tracker system ID");
      }
      groups[{system, ids.at("layer")}].push_back(element);
    }
    for (const auto& [name, child] : element.children()) {
      (void)name;
      visit(child);
    }
  };
  visit(detector.world());
  for (int system = 1; system <= 7; ++system) {
    if (std::none_of(groups.begin(), groups.end(), [system](const auto& group) {
          return group.first.first == system;
        })) {
      throw std::runtime_error("Gen3 requires complete nODD systems 1–7");
    }
  }
  Blueprint::Config cfg;
  cfg.envelope = ExtentEnvelope{Envelope{envelope, envelope}};
  Blueprint root{cfg};
  auto& outer = root.addCylinderContainer("NoddTracker", AxisR);
  outer.setAttachmentStrategy(VolumeAttachmentStrategy::Gap);
  outer.addChild(builder.backend().makeBeampipe());
  for (int barrelSystem : {1, 4, 6}) {
    auto& subsystem = outer.addCylinderContainer(
        "Subsystem" + std::to_string(barrelSystem), AxisZ);
    subsystem.setAttachmentStrategy(VolumeAttachmentStrategy::Gap);
    subsystem.setResizeStrategies(VolumeResizeStrategy::Gap,
                                  VolumeResizeStrategy::Gap);
    auto addLayer = [&](const auto& key, const auto& sensors,
                        BlueprintNode& parent) {
      auto layer =
          builder.layerFromSensors()
              .setSensorAxes("XYZ")
              .setSensors(sensors)
              .setLayerName("System" + std::to_string(key.first) + "Layer" +
                            std::to_string(key.second))
              .setEnvelope(ExtentEnvelope{Envelope{envelope, envelope}})
              .onLayer([](const std::optional<Element>&,
                          LayerBlueprintNode& node) {
                // Correctness-first policy supports staggered chips,
                // tilted strips and paired faces without bin assumptions.
                node.setUseCenterOfGravity(false, false, true);
                node.setNavigationPolicyFactory(
                    NavigationPolicyFactory{}
                        .add<TryAllNavigationPolicy>()
                        .asUniquePtr());
                for (const auto& surface : node.surfaces()) {
                  const auto* source =
                      dynamic_cast<const ActsPlugins::DD4hepDetectorElement*>(
                          surface->surfacePlacement());
                  if (source == nullptr) {
                    throw std::runtime_error("Missing DD4hep-backed sensor");
                  }
                  auto material = source->sourceElement()
                                      .volume()
                                      .material()
                                      .ptr()
                                      ->GetMaterial();
                  surface->assignSurfaceMaterial(
                      nativeMaterial(*material, source->thickness()));
                }
              });
      if (key.first == barrelSystem) {
        std::move(layer).barrel().addTo(parent);
      } else {
        std::move(layer).endcap().addTo(parent);
      }
    };
    std::vector<std::pair<std::pair<int, int>, std::vector<Element>>> endcaps;
    for (const auto& [key, sensors] : groups) {
      if ((barrelSystem == 1 ? (key.first == 2 || key.first == 3)
                             : key.first == barrelSystem + 1)) {
        endcaps.emplace_back(key, sensors);
      }
    }
    auto z = [](const auto& group) {
      return group.second.front()
          .nominal()
          .worldTransformation()
          .GetTranslation()[2];
    };
    std::ranges::sort(
        endcaps, [&](const auto& a, const auto& b) { return z(a) < z(b); });
    for (const auto& group : endcaps) {
      if (z(group) < 0) {
        addLayer(group.first, group.second, subsystem);
      }
    }
    auto& barrel = subsystem.addCylinderContainer(
        "Barrel" + std::to_string(barrelSystem), AxisR);
    barrel.setAttachmentStrategy(VolumeAttachmentStrategy::Gap);
    barrel.setResizeStrategies(VolumeResizeStrategy::Gap,
                               VolumeResizeStrategy::Gap);
    for (const auto& [key, sensors] : groups) {
      if (key.first == barrelSystem) {
        addLayer(key, sensors, barrel);
      }
    }
    for (const auto& group : endcaps) {
      if (z(group) > 0) {
        addLayer(group.first, group.second, subsystem);
      }
    }
  }
  auto geometry = root.construct(BlueprintOptions{}, gctx, logger);
  const auto pipe = detector.detector("BeamPipe");
  const auto* tube = dynamic_cast<const TGeoTube*>(pipe.volume().solid().ptr());
  if (tube == nullptr) {
    throw std::runtime_error("BeamPipe is not a Tube");
  }
  const double radius =
      (tube->GetRmin() + tube->GetRmax()) * UnitConstants::cm / 2;
  const double thickness =
      (tube->GetRmax() - tube->GetRmin()) * UnitConstants::cm;
  auto material =
      nativeMaterial(*pipe.volume().material().ptr()->GetMaterial(), thickness);
  std::set<const Surface*> assigned;
  geometry->highestTrackingVolume()->visitVolumes(
      [&](const TrackingVolume* volume) {
        for (auto& portal : const_cast<TrackingVolume*>(volume)->portals()) {
          auto& surface = portal.surface();
          auto bounds = dynamic_cast<const CylinderBounds*>(&surface.bounds());
          if (bounds != nullptr &&
              std::abs(bounds->get(CylinderBounds::eR) - radius) < 1e-8) {
            surface.assignSurfaceMaterial(material);
            assigned.insert(&surface);
          }
        }
      });
  if (assigned.size() != 1) {
    throw std::runtime_error(
        "Expected one shared finite cylindrical pipe portal");
  }
  return geometry;
}
}  // namespace

DD4hepNoddDetector::Config::Config() {
  name = "NoddDetector";
}

DD4hepNoddDetector::DD4hepNoddDetector(const Config& cfg)
    : DD4hepDetectorBase{cfg}, m_cfg{cfg} {
  if (cfg.gen3) {
    if (!std::isfinite(cfg.navigationEnvelopeMm) ||
        cfg.navigationEnvelopeMm <= 0) {
      throw std::invalid_argument(
          "navigationEnvelopeMm must be finite and positive");
    }
    m_trackingGeometry =
        buildGeometry(dd4hepDetector(), nominalGeometryContext(), logger(),
                      cfg.navigationEnvelopeMm);
  }
}

auto DD4hepNoddDetector::config() const -> const Config& {
  return m_cfg;
}

std::string DD4hepNoddDetector::geometryReport() const {
  if (m_trackingGeometry == nullptr) {
    throw std::runtime_error("Gen3 geometry was not requested");
  }
  using nlohmann::json;
  json report = {{"generation", 3},
                 {"sensors", json::array()},
                 {"volumes", json::array()},
                 {"beampipe_portals", json::array()}};
  const auto gctx = nominalGeometryContext();
  m_trackingGeometry->visitSurfaces([&](const Acts::Surface* surface) {
    const auto* element =
        dynamic_cast<const ActsPlugins::DD4hepDetectorElement*>(
            surface->surfacePlacement());
    if (element == nullptr) {
      return;
    }
    auto bounds =
        dynamic_cast<const Acts::RectangleBounds*>(&surface->bounds());
    if (bounds == nullptr) {
      throw std::runtime_error("Non-rectangular nODD sensor");
    }
    const auto transform = surface->localToGlobalTransform(gctx);
    const auto center = surface->center(gctx);
    json item = {{"ids", placementIds(element->sourceElement())},
                 {"geometry_id", surface->geometryId().value()},
                 {"source_path", element->sourceElement().path()},
                 {"center_mm", {center.x(), center.y(), center.z()}},
                 {"size_mm",
                  {2 * bounds->halfLengthX(), 2 * bounds->halfLengthY(),
                   element->thickness()}}};
    for (int i = 0; i < 3; ++i) {
      const auto axis = transform.linear().col(i);
      item[i == 0   ? "u"
           : i == 1 ? "v"
                    : "normal"] = {axis.x(), axis.y(), axis.z()};
    }
    report["sensors"].push_back(item);
  });
  m_trackingGeometry->visitVolumes([&](const Acts::TrackingVolume* volume) {
    if (volume->volumeName() == "BeamPipe") {
      for (const auto& portal : volume->portals()) {
        const auto& surface = portal.surface();
        const auto* bounds =
            dynamic_cast<const Acts::CylinderBounds*>(&surface.bounds());
        const auto* material = surface.surfaceMaterial();
        if (bounds != nullptr && material != nullptr) {
          const auto& slab = material->materialSlab(Acts::Vector2::Zero());
          report["beampipe_portals"].push_back(
              {{"geometry_id", surface.geometryId().value()},
               {"radius_mm", bounds->get(Acts::CylinderBounds::eR)},
               {"half_length_mm",
                bounds->get(Acts::CylinderBounds::eHalfLengthZ)},
               {"thickness_mm", slab.thickness()},
               {"radiation_length_mm", slab.material().X0()},
               {"elemental_Z", slab.material().Z()}});
        }
      }
    }
    report["volumes"].push_back(
        {{"name", volume->volumeName()},
         {"geometry_id", volume->geometryId().value()},
         {"portals", std::ranges::distance(volume->portals())},
         {"surfaces", std::ranges::distance(volume->surfaces())}});
  });
  return report.dump();
}
}  // namespace ActsExamples
