// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Geometry/CylinderVolumeBounds.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Geometry/PortalLinkBase.hpp"
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Geometry/TrackingVolume.hpp"
#include "Acts/Navigation/SurfaceArrayNavigationPolicy.hpp"
#include "Acts/Surfaces/CylinderBounds.hpp"
#include "Acts/Surfaces/DiscBounds.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "ActsPlugins/ActSVG/LayerSvgConverter.hpp"
#include "ActsPlugins/ActSVG/SurfaceArraySvgConverter.hpp"
#include "ActsPlugins/ActSVG/SurfaceSvgConverter.hpp"
#include "ActsPlugins/ActSVG/SvgUtils.hpp"
#include "ActsPlugins/ActSVG/TrackingGeometrySvgConverter.hpp"
#include "ActsPython/Utilities/Helpers.hpp"
#include "ActsPython/Utilities/Macros.hpp"
#include <actsvg/core/draw.hpp>

#include <algorithm>
#include <memory>
#include <ranges>
#include <sstream>
#include <string>
#include <tuple>
#include <vector>

#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

namespace py = pybind11;
using namespace pybind11::literals;

using namespace Acts;
using namespace ActsExamples;
using namespace ActsPlugins;

PYBIND11_MODULE(ActsPluginsPythonBindingsSvg, svg) {
  using namespace Acts;
  using namespace ActsPlugins;

  // Primitives, should be dropped in favour of actsvg pybind11 bindings
  {
    py::class_<actsvg::svg::object>(svg, "object")
        .def_readwrite("id", &actsvg::svg::object::_id);

    py::class_<actsvg::svg::file>(svg, "file")
        .def(py::init<>())
        .def("add_object", &actsvg::svg::file::add_object)
        .def("add_objects", &actsvg::svg::file::add_objects)
        .def("clip",
             [](actsvg::svg::file& self, std::array<actsvg::scalar, 4> box) {
               self.set_view_box(box);
             })
        .def("write",
             [](const actsvg::svg::file& self, const std::string& filename) {
               std::ofstream file(filename);
               file << self;
               file.close();
             });

    svg.def("toFile", &Svg::toFile, py::arg("objects"), py::arg("filename"));
  }

  // Core components, added as an acts.svg submodule
  {
    auto c = py::class_<Svg::Style>(svg, "Style").def(py::init<>());
    ACTS_PYTHON_STRUCT(c, fillColor, fillOpacity, highlightColor, highlights,
                       strokeWidth, strokeColor, highlightStrokeWidth,
                       highlightStrokeColor, fontSize, fontColor,
                       quarterSegments);
  }

  // How surfaces should be drawn: Svg Surface options & drawning
  {
    auto c = py::class_<Svg::SurfaceConverter::Options>(svg, "SurfaceOptions")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, style, templateSurface);

    // Define the proto surface
    py::class_<Svg::ProtoSurface>(svg, "ProtoSurface");
    // Convert an Surface object into an svg::proto::surface

    svg.def("convertSurface", &Svg::SurfaceConverter::convert);

    // Define the view functions
    svg.def("viewSurface", [](const Svg::ProtoSurface& pSurface,
                              const std::string& identification,
                              const std::string& view = "xy") {
      if (view == "xy") {
        return Svg::View::xy(pSurface, identification);
      } else if (view == "zr") {
        return Svg::View::zr(pSurface, identification);
      } else if (view == "zphi") {
        return Svg::View::zphi(pSurface, identification);
      } else if (view == "zrphi") {
        return Svg::View::zrphi(pSurface, identification);
      } else {
        throw std::invalid_argument("Unknown view type");
      }
    });
  }

  // Draw primitives
  {
    svg.def("drawArrow", &actsvg::draw::arrow);

    svg.def("drawText", &actsvg::draw::text);

    svg.def("drawInfoBox", &Svg::infoBox);
  }

  // Draw Eta Lines
  {
    svg.def(
        "drawEtaLines",
        [](const std::string& id, actsvg ::scalar z, actsvg::scalar r,
           const std::vector<actsvg::scalar>& etaMain,
           actsvg::scalar strokeWidthMain, unsigned int sizeMain,
           bool labelMain, const std::vector<actsvg::scalar>& etaSub,
           actsvg::scalar strokeWidthSub, const std::vector<int>& strokeDashSub,
           unsigned int sizeSub, bool labelSub) {
          // The main eta lines
          actsvg::style::stroke strokeMain;
          strokeMain._width = strokeWidthMain;
          actsvg::style::font fontMain;
          fontMain._size = sizeMain;

          actsvg::style::stroke strokeSub;
          strokeSub._width = strokeWidthSub;
          strokeSub._dasharray = strokeDashSub;
          actsvg::style::font fontSub;
          fontSub._size = sizeSub;

          return actsvg::display::eta_lines(
              id, z, r,
              {std::tie(etaMain, strokeMain, labelMain, fontMain),
               std::tie(etaSub, strokeSub, labelSub, fontSub)});
        });
  }

  // How detector volumes are drawn: Svg DetectorVolume options & drawning
  {
    py::class_<Svg::ProtoVolume>(svg, "ProtoVolume");

    py::class_<Svg::ProtoGrid>(svg, "ProtoGrid");

    py::class_<Svg::ProtoIndexedSurfaceGrid>(svg, "ProtoIndexedSurfaceGrid");
  }

  { svg.def("drawSurfaceArrays", &Svg::drawSurfaceArrays); }

  // Legacy geometry drawing
  {
    using DefinedStyle = std::pair<GeometryIdentifier, Svg::Style>;
    using DefinedStyleSet = std::vector<DefinedStyle>;

    auto sm = py::class_<GeometryHierarchyMap<Svg::Style>>(svg, "StyleMap")
                  .def(py::init<DefinedStyleSet>(), py::arg("elements"));

    auto c = py::class_<Svg::LayerConverter::Options>(svg, "LayerOptions")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, name, surfaceStyles, zRange, phiRange, gridInfo,
                       moduleInfo, projectionInfo, labelProjection, labelGauge);
  }

  {
    using DefinedLayerOptions =
        std::pair<GeometryIdentifier, Svg::LayerConverter::Options>;
    using DefinedLayerOptionsSet = std::vector<DefinedLayerOptions>;

    auto lom =
        py::class_<GeometryHierarchyMap<Svg::LayerConverter::Options>>(
            svg, "LayerOptionMap")
            .def(py::init<DefinedLayerOptionsSet>(), py::arg("elements"));

    auto c = py::class_<Svg::TrackingGeometryConverter::Options>(
                 svg, "TrackingGeometryOptions")
                 .def(py::init<>());
    ACTS_PYTHON_STRUCT(c, prefix, layerOptions);
  }

  // Tracking geometry drawing
  {
    svg.def(
        "drawTrackingGeometry",
        [](const GeometryContext& gctx, const TrackingGeometry& tGeometry,
           const std::string& view, bool drawSurfaces, bool highlightMaterial) {
          std::variant<actsvg::views::x_y, actsvg::views::z_r> v;
          if (view == "xy") {
            v = actsvg::views::x_y();
          } else if (view == "zr") {
            v = actsvg::views::z_r();
          } else {
            throw std::invalid_argument("Unknown view type");
          }

          return Svg::drawTrackingGeometry(gctx, tGeometry, v, drawSurfaces,
                                           highlightMaterial);
        },
        py::arg("gctx"), py::arg("tGeometry"), py::arg("view"),
        py::arg("drawSurfaces") = true, py::arg("highlightMaterial") = false);
  }
}
