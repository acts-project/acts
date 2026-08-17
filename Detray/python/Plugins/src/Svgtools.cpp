// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Detray include(s)
#include "detray/algebra/concepts.hpp"
#include "detray/core/detector.hpp"
#include "detray/plugins/svgtools/illustrator.hpp"
#include "detray/plugins/svgtools/styling/styling.hpp"
#include "detray/plugins/svgtools/writer.hpp"

// Detray algebra plugin + detector metadata
#include "algebra/array.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/detectors/default_metadata.hpp"

// Actsvg include(s)
#include <actsvg/core/defs.hpp>
#include <actsvg/core/draw.hpp>
#include <actsvg/core/style.hpp>
#include <actsvg/core/svg.hpp>
#include <actsvg/core/views.hpp>

// Pybind11 include(s)
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>

// System include(s)
#include <array>
#include <string>
#include <vector>

namespace py = pybind11;

namespace {

using scalar_t = DETRAY_CUSTOM_SCALARTYPE;
using algebra_t = detray::array<scalar_t>;
using detector_t = detray::detector<detray::default_metadata<algebra_t>>;
using geometry_context_t = detector_t::geometry_context;
using style_t = detray::svgtools::styling::style;
using illustrator_t = detray::svgtools::illustrator<detector_t>;

/// Bind the draw methods of the illustrator for @tparam view_t .
template <typename view_t>
void bind_draw_methods(py::class_<illustrator_t> &c) {
  c.def(
      "drawDetector",
      [](const illustrator_t &il, const view_t &view,
         const geometry_context_t &gctx, float r_length, float z_length) {
        return il.draw_detector(view, gctx, r_length, z_length);
      },
      py::arg("view"), py::arg("gctx"), py::arg("rLength") = 1100.f,
      py::arg("zLength") = 3100.f, "Draw the whole detector");

  c.def(
      "drawVolumes",
      [](const illustrator_t &il, const std::vector<detray::dindex> &indices,
         const view_t &view, const geometry_context_t &gctx) {
        return il.draw_volumes(indices, view, gctx);
      },
      py::arg("indices"), py::arg("view"), py::arg("gctx"),
      "Draw the volumes with the given indices. Returns their svg together "
      "with the svgs of their surface grids");
  c.def(
      "drawVolumes",
      [](const illustrator_t &il, const std::vector<std::string> &names,
         const view_t &view, const geometry_context_t &gctx) {
        return il.draw_volumes(names, view, gctx);
      },
      py::arg("names"), py::arg("view"), py::arg("gctx"),
      "Draw the volumes with the given names. Returns their svg together with "
      "the svgs of their surface grids");

  c.def(
      "drawSurfaces",
      [](const illustrator_t &il, const std::vector<detray::dindex> &indices,
         const view_t &view, const geometry_context_t &gctx) {
        return il.draw_surfaces(indices, view, gctx);
      },
      py::arg("indices"), py::arg("view"), py::arg("gctx"),
      "Draw the surfaces with the given indices. Returns their svg together "
      "with the svg of their material");
}

}  // namespace

PYBIND11_MODULE(DetraySvgtoolsPythonBindings, m) {
  m.doc() = "Detray svgtools bindings";

  // Svg primitives
  {
    // Duplicates actsvg's binding. Registered module-locally, so that importing
    // this module alongside actsvg does not clash over the same C++ type
    py::class_<actsvg::svg::object>(m, "SvgObject", py::module_local())
        .def_readwrite("id", &actsvg::svg::object::_id,
                       "Identification string");

    m.def(
        "writeSvg",
        [](const std::string &path,
           const std::vector<actsvg::svg::object> &svgs,
           bool replace) { detray::svgtools::write_svg(path, svgs, replace); },
        py::arg("path"), py::arg("svgs"), py::arg("replace") = true,
        "Write svg objects to '<path>.svg'. Their ids must be unique");

    m.def(
        "drawAxes",
        [](const std::string &name,
           const std::array<actsvg::scalar, 2> &x_range,
           const std::array<actsvg::scalar, 2> &y_range,
           const std::string &x_label, const std::string &y_label,
           unsigned int font_size) {
          actsvg::style::font font{};
          font._size = font_size;

          // The default stroke is black
          return actsvg::draw::x_y_axes(name, x_range, y_range,
                                        actsvg::style::stroke{}, x_label,
                                        y_label, font);
        },
        py::arg("name"), py::arg("xRange"), py::arg("yRange"),
        py::arg("xLabel") = "", py::arg("yLabel") = "",
        py::arg("fontSize") = 12u, "Draw a pair of x-y axes");
  }

  // Svg views
  {
    py::class_<actsvg::views::x_y>(m, "ViewXY", py::module_local())
        .def(py::init<>(), "Transverse view");
    py::class_<actsvg::views::z_r>(m, "ViewZR", py::module_local())
        .def(py::init<>(), "Longitudinal view");
    py::class_<actsvg::views::z_phi>(m, "ViewZPhi", py::module_local())
        .def(py::init<>(), "Unrolled cylinder view");
    py::class_<actsvg::views::z_rphi>(m, "ViewZRPhi", py::module_local())
        .def(py::init<>(), "Unrolled cylinder view at a fixed radius");
  }

  // Svg styling
  {
    py::class_<style_t>(m, "Style")
        .def_property(
            "fontSize",
            [](const style_t &s) {
              return s._detector_style._volume_style._sensitive_surface_style
                  ._font_size;
            },
            [](style_t &s, unsigned int font_size) {
              auto &vol_style = s._detector_style._volume_style;
              vol_style._sensitive_surface_style._font_size = font_size;
              vol_style._passive_surface_style._font_size = font_size;
              vol_style._portal_style._surface_style._font_size = font_size;
              s._eta_lines_style._font_size = font_size;
            },
            "Font size of all surface and eta line labels")
        .def_property(
            "materialGradientPos",
            [](const style_t &s) {
              return s._detector_style._volume_style._sensitive_surface_style
                  ._material_style._gradient_pos;
            },
            [](style_t &s, const actsvg::point2 &pos) {
              auto &vol_style = s._detector_style._volume_style;
              vol_style._sensitive_surface_style._material_style._gradient_pos =
                  pos;
              vol_style._passive_surface_style._material_style._gradient_pos =
                  pos;
              vol_style._portal_style._surface_style._material_style
                  ._gradient_pos = pos;
            },
            "Position of the material gradient box in all surface styles");

    m.def(
        "tableauColorblindStyle",
        []() { return detray::svgtools::styling::tableau_colorblind::style; },
        "Style that matches the colors used in the data plotting");
  }

  // Svg illustrator
  {
    auto c =
        py::class_<illustrator_t>(m, "Illustrator")
            .def(py::init<const detector_t &, const detray::name_map &,
                          const style_t &>(),
                 py::arg("detector"), py::arg("names"), py::arg("style"),
                 py::keep_alive<1, 2>(), py::keep_alive<1, 3>(),
                 "Build the svg generator of a detector. The detector and its "
                 "names are kept alive for as long as the illustrator is used, "
                 "as it only references them. The style is copied")
            .def("showInfo", &illustrator_t::show_info,
                 py::arg("toggle") = true, "Toggle the information boxes")
            .def("hideEtaLines", &illustrator_t::hide_eta_lines,
                 py::arg("toggle") = true,
                 "Toggle the eta lines of the z-r view")
            .def("hideGrids", &illustrator_t::hide_grids,
                 py::arg("toggle") = true, "Toggle the surface grids")
            .def("hideMaterial", &illustrator_t::hide_material,
                 py::arg("toggle") = true, "Toggle the surface material")
            .def("hidePortals", &illustrator_t::hide_portals,
                 py::arg("toggle") = true, "Toggle the portal surfaces")
            .def("hidePassives", &illustrator_t::hide_passives,
                 py::arg("toggle") = true, "Toggle the passive surfaces")
            .def("searchWindow", &illustrator_t::search_window,
                 py::arg("window"),
                 "Neighborhood search window for the grid display")
            .def_property_readonly("detName", &illustrator_t::det_name,
                                   "Name of the detector");

    bind_draw_methods<actsvg::views::x_y>(c);
    bind_draw_methods<actsvg::views::z_r>(c);
    bind_draw_methods<actsvg::views::z_phi>(c);
    bind_draw_methods<actsvg::views::z_rphi>(c);
  }
}
