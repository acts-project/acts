// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/ActSVG/TrackingGeometrySvgConverter.hpp"

#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/GridPortalLink.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Plugins/ActSVG/DetectorVolumeSvgConverter.hpp"
#include "Acts/Plugins/ActSVG/PortalSvgConverter.hpp"
#include "Acts/Plugins/ActSVG/SurfaceSvgConverter.hpp"
#include "Acts/Utilities/Zip.hpp"

#include <sstream>
#include <stdexcept>

namespace Acts::Svg {

std::vector<actsvg::svg::object> TrackingGeometryConverter::convert(
    const GeometryContext& gctx, const TrackingGeometry& tGeometry,
    const TrackingGeometryConverter::Options& cOptions) {
  // Get the world volume
  const TrackingVolume* world = tGeometry.highestTrackingVolume();

  // Initiate the cache
  TrackingGeometryConverter::State cState;

  // Run the conversion recursively
  convert(gctx, *world, cOptions, cState);

  // Digest the views and globals
  std::vector<actsvg::svg::object> finalViews = cState.finalViews;
  if (!cState.xyCrossSection.empty()) {
    finalViews.push_back(
        Acts::Svg::group(cState.xyCrossSection, cOptions.prefix + "layers_xy"));
  }
  if (!cState.zrCrossSection.empty()) {
    finalViews.push_back(
        Acts::Svg::group(cState.zrCrossSection, cOptions.prefix + "layers_zr"));
  }
  // return all final Views
  return finalViews;
}

void TrackingGeometryConverter::convert(
    const GeometryContext& gctx, const TrackingVolume& tVolume,
    const TrackingGeometryConverter::Options& cOptions,
    TrackingGeometryConverter::State& cState) {
  // Process confined layers first
  if (tVolume.confinedLayers() != nullptr) {
    for (const auto& layer : tVolume.confinedLayers()->arrayObjects()) {
      if (layer->surfaceArray() != nullptr) {
        GeometryIdentifier geoID = layer->geometryId();
        std::string layerName = cOptions.prefix + "vol_" +
                                std::to_string(geoID.volume()) + "_layer_" +
                                std::to_string(geoID.layer());

        LayerConverter::Options lOptions;
        // Search for predefined layer options
        auto deflOptions = cOptions.layerOptions.find(geoID);
        if (deflOptions != cOptions.layerOptions.end()) {
          lOptions = (*deflOptions);
          lOptions.name = layerName;
        }
        // Retrieve the layer sheets
        auto layerSheets = LayerConverter::convert(gctx, *layer, lOptions);

        // Record the sheets
        for (const auto& lSheet : layerSheets) {
          if (lSheet.is_defined()) {
            cState.finalViews.push_back(lSheet);
          }
        }
        // Collect the xy views
        if (layerSheets[LayerConverter::eCrossSectionXY].is_defined() &&
            layer->surfaceRepresentation().type() == Acts::Surface::Cylinder) {
          cState.xyCrossSection.push_back(
              layerSheets[LayerConverter::eCrossSectionXY]);
        }
        // Collect the zr views
        if (layerSheets[LayerConverter::eCrossSectionZR].is_defined()) {
          cState.zrCrossSection.push_back(
              layerSheets[LayerConverter::eCrossSectionZR]);
        }
      }
    }
  }

  // Run recursively over confined volumes
  if (tVolume.confinedVolumes() != nullptr) {
    for (const auto& volume : tVolume.confinedVolumes()->arrayObjects()) {
      convert(gctx, *volume, cOptions, cState);
    }
  }
}

std::array<actsvg::svg::object, 2> TrackingGeometryProjections::convert(
    const GeometryContext& gctx, const Acts::TrackingGeometry& tGeometry,
    const TrackingGeometryProjections::Options& cOptions) {
  // The projections
  actsvg::svg::object xyView;
  actsvg::svg::object zrView;

  // Get the world volume
  const Acts::TrackingVolume* world = tGeometry.highestTrackingVolume();
  if (world != nullptr) {
    // Initiate the cache
    Acts::Svg::TrackingGeometryConverter::State cState;

    // Run the conversion recursively
    Acts::Svg::TrackingGeometryConverter::convert(
        gctx, *world, cOptions.trackingGeometryOptions, cState);

    xyView = Acts::Svg::group(cState.xyCrossSection,
                              cOptions.prefix + "projection_xy");
    zrView = Acts::Svg::group(cState.zrCrossSection,
                              cOptions.prefix + "projection_zr");
  }
  return {xyView, zrView};
}

// Gen3

namespace {
void convertPortalLink(const GeometryContext& gctx,
                       const PortalLinkBase& portalLink, Direction direction,
                       std::vector<Svg::ProtoPortal::link>& links) {
  auto getCenters = [&](const GridPortalLink& grid, AxisDirection axis0) {
    assert((grid.dim() == 2 || grid.dim() == 1) &&
           "Grid has unexpected number of dimension");

    std::vector<Vector2> centers;

    const auto localBins = grid.grid().numLocalBinsAny();
    if (grid.dim() == 2) {
      for (std::size_t i = 1; i <= localBins[0]; ++i) {
        for (std::size_t j = 1; j <= localBins[1]; ++j) {
          auto center = grid.grid().binCenterAny({i, j});
          assert(center.size() == 2);
          centers.push_back(Vector2(center[0], center[1]));
        }
      }
    } else {
      for (std::size_t i = 1; i <= localBins[0]; ++i) {
        auto center = grid.grid().binCenterAny({i});
        assert(center.size() == 1);
        if (axis0 == grid.direction()) {
          centers.push_back(Vector2(center[0], 0.));
        } else {
          centers.push_back(Vector2(0., center[0]));
        }
      }
    }

    return centers;
  };

  if (const auto* discBounds =
          dynamic_cast<const DiscBounds*>(&portalLink.surface().bounds());
      discBounds != nullptr) {
    if (const auto* grid = dynamic_cast<const GridPortalLink*>(&portalLink);
        grid != nullptr) {
      std::vector<Vector2> centers = getCenters(*grid, AxisDirection::AxisR);
      for (const auto& center : centers) {
        Svg::ProtoPortal::link link;
        link._start = portalLink.surface().localToGlobal(gctx, center);
        Vector3 normal = portalLink.surface().normal(gctx, link._start);
        link._end = link._start + normal * 10. * direction.sign();
        links.push_back(link);
      }
    } else {
      double rMin = discBounds->rMin();
      double rMax = discBounds->rMax();
      double rMid = 0.5 * (rMin + rMax);

      Svg::ProtoPortal::link link;
      link._start = portalLink.surface().center(gctx) + Vector3(rMid, 0., 0.);
      Vector3 normal = portalLink.surface().normal(gctx, link._start);
      link._end = link._start + normal * 10. * direction.sign();
      links.push_back(link);
    }

  } else if (const auto* cylBounds = dynamic_cast<const CylinderBounds*>(
                 &portalLink.surface().bounds());
             cylBounds != nullptr) {
    if (const auto* grid = dynamic_cast<const GridPortalLink*>(&portalLink);
        grid != nullptr) {
      std::vector<Vector2> centers = getCenters(*grid, AxisDirection::AxisRPhi);
      for (const auto& center : centers) {
        Svg::ProtoPortal::link link;
        link._start = portalLink.surface().localToGlobal(gctx, center);
        Vector3 normal = portalLink.surface().normal(gctx, link._start);
        link._end = link._start + normal * 10. * direction.sign();
        links.push_back(link);
      }
    } else {
      double r = cylBounds->get(CylinderBounds::eR);
      Svg::ProtoPortal::link link;
      link._start = portalLink.surface().center(gctx) + Vector3(r, 0., 0.);
      Vector3 normal = portalLink.surface().normal(gctx, link._start);
      link._end = link._start + normal * 10. * direction.sign();
      links.push_back(link);
    }
  } else {
    throw std::invalid_argument("Unknown bounds type");
  }
}

struct Visitor : TrackingGeometryVisitor {
  explicit Visitor(const GeometryContext& gctxIn) : gctx(gctxIn) {}

  void visitSurface(const Surface& surface) override {
    auto proto = Svg::SurfaceConverter::convert(gctx, surface, {});
    surfaces.push_back(std::pair{&surface, proto});
  }

  void visitPortal(const Portal& portal) override {
    auto pIt = std::ranges::find_if(
        portals, [&](const auto& p) { return p.first == &portal; });

    Svg::ProtoPortal pPortal;
    if (pIt == portals.end()) {
      pPortal._name = "portal_" + std::to_string(portals.size());
      auto surface = Svg::SurfaceConverter::convert(gctx, portal.surface(), {});
      pPortal._surface = surface;

      if (auto link = portal.getLink(Direction::AlongNormal());
          link != nullptr) {
        convertPortalLink(gctx, *link, Direction::AlongNormal(),
                          pPortal._volume_links);
      }

      if (auto link = portal.getLink(Direction::OppositeNormal());
          link != nullptr) {
        convertPortalLink(gctx, *link, Direction::OppositeNormal(),
                          pPortal._volume_links);
      }
      portals.push_back(std::pair{&portal, pPortal});
    } else {
      pPortal = pIt->second;
    }

    auto& [volume, pVolume] = volumes.back();
    pVolume._portals.push_back(pPortal);
  }

  void visitVolume(const TrackingVolume& volume) override {
    Svg::ProtoVolume pVolume;
    pVolume._name = volume.volumeName();

    auto volumeTransform = volume.transform();

    // https://github.com/acts-project/actsvg/blob/2f1aaa58365a1dd1af62dc27aea5039459a65a38/meta/include/actsvg/display/geometry.hpp#L687-L692
    enum svgBv : unsigned int {
      rInner = 0u,
      rOuter = 1u,
      zPos = 2u,
      zHalf = 3u,
      phiSec = 4u,
      avgPhi = 5u,
    };

    using enum CylinderVolumeBounds::BoundValues;
    const auto& boundValues = volume.volumeBounds().values();
    if (volume.volumeBounds().type() == Acts::VolumeBounds::eCylinder) {
      pVolume._bound_values.resize(6);
      pVolume._bound_values.at(svgBv::rInner) =
          static_cast<actsvg::scalar>(boundValues[eMinR]);
      pVolume._bound_values.at(svgBv::rOuter) =
          static_cast<actsvg::scalar>(boundValues[eMaxR]);
      pVolume._bound_values.at(svgBv::zPos) =
          static_cast<actsvg::scalar>(volumeTransform.translation().z());
      pVolume._bound_values.at(svgBv::zHalf) =
          static_cast<actsvg::scalar>(boundValues[eHalfLengthZ]);
      pVolume._bound_values.at(svgBv::phiSec) =
          static_cast<actsvg::scalar>(boundValues[eHalfPhiSector]);
      pVolume._bound_values.at(svgBv::avgPhi) =
          static_cast<actsvg::scalar>(boundValues[eAveragePhi]);
    } else {
      throw std::invalid_argument("Unknown bounds type");
    }

    volumes.push_back(std::pair{&volume, pVolume});
  }

  const GeometryContext& gctx;

  std::vector<std::pair<const Portal*, Svg::ProtoPortal>> portals;

  std::vector<std::pair<const TrackingVolume*, Svg::ProtoVolume>> volumes;

  std::vector<std::pair<const Surface*, Svg::ProtoSurface>> surfaces;
};
}  // namespace

std::vector<actsvg::svg::object> drawTrackingGeometry(
    const GeometryContext& gctx, const TrackingGeometry& tGeometry,
    std::variant<actsvg::views::x_y, actsvg::views::z_r> view,
    bool drawSurfaces, bool highlightMaterial) {
  if (tGeometry.geometryVersion() != TrackingGeometry::GeometryVersion::Gen3) {
    throw std::invalid_argument{
        "Input tracking geometry needs to have been built in Gen3 mode"};
  }

  Visitor visitor(gctx);
  tGeometry.apply(visitor);

  std::vector<actsvg::svg::object> objects;

  std::visit(
      [&](auto& _view) {
        for (const auto& [tv, volume] : visitor.volumes) {
          std::string id = tv->volumeName();
          auto object = actsvg::display::volume(id, volume, _view);
          objects.push_back(object);

          std::vector<std::string> lines;

          std::stringstream ss;
          ss << "ID: " << tv->geometryId();
          lines.push_back(ss.str());

          ss.str("");
          ss << tv->volumeBounds();
          std::string bounds = ss.str();
          // split at first space after 40 characters
          ss.str("");
          for (std::size_t i = 0; i < bounds.size(); ++i) {
            ss << bounds[i];
            if (ss.str().size() > 40 && bounds[i] == ' ') {
              lines.push_back(ss.str());
              ss.str("");
            }
          }

          auto text = actsvg::draw::connected_info_box(
              "info_volume_" + volume._name, {0, 0}, volume._name,
              {{._rgb{200, 200, 200}}}, {._fc{._rgb{0, 0, 0}}, ._size = 24},
              lines, {{._rgb{220, 220, 220}}},
              {._fc{._rgb{0, 0, 0}}, ._size = 24}, {}, object);

          objects.push_back(text);
        }

        if (drawSurfaces) {
          for (const auto& [surface, proto] : visitor.surfaces) {
            std::stringstream ss;
            ss << surface->geometryId();
            auto object = actsvg::display::surface(ss.str(), proto, _view);

            if (highlightMaterial && surface->surfaceMaterial() != nullptr) {
              object._stroke._sc._rgb = {255, 0, 0};
              object._stroke._width = 1.5;
            }

            objects.push_back(object);
          }
        }
      },
      view);

  return objects;
}

}  // namespace Acts::Svg
