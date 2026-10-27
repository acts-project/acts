// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/plugins/svgtools/meta/proto/eta_lines.hpp"
#include "detray/plugins/svgtools/meta/proto/intersection.hpp"
#include "detray/plugins/svgtools/meta/proto/landmark.hpp"
#include "detray/plugins/svgtools/meta/proto/trajectory.hpp"
#include "detray/plugins/svgtools/styling/colors.hpp"

// Actsvg include(s)
#include "actsvg/core.hpp"
#include "actsvg/core/style.hpp"
#include "actsvg/proto/detector.hpp"
#include "actsvg/proto/grid.hpp"
#include "actsvg/proto/material.hpp"
#include "actsvg/proto/portal.hpp"
#include "actsvg/proto/surface.hpp"
#include "actsvg/proto/volume.hpp"
#include "actsvg/styles/defaults.hpp"

// System include(s)
#include <array>
#include <cstdlib>
#include <vector>

namespace detray::svgtools::styling {

/// Style applied to an actsvg proto grid
struct grid_style {
  /// Color used to draw grid lines.
  actsvg::style::color _stroke_color{};
  /// Width of grid lines.
  actsvg::scalar _stroke_width{};
};

/// Style applied to an actsvg proto surface_material
struct surface_material_style {
  /// Position of the material information annotation.
  actsvg::point2 _info_pos{};

  /// Position of the material gradient.
  actsvg::point2 _gradient_pos{};
  /// Dimensions of the material gradient box.
  std::array<actsvg::scalar, 2u> _gradient_box{};
  /// Font used for the material gradient.
  actsvg::style::font _gradient_font{};
  /// Label shown with the material gradient.
  actsvg::style::label _gradient_label{};
  /// Color of the material gradient outline.
  actsvg::style::color _gradient_stroke_color{};
  /// Width of the material gradient outline.
  actsvg::scalar _gradient_stroke_width{0.f};
  /// Colors sampled to construct the material gradient.
  std::vector<actsvg::style::color> _gradient_color_scale{};
  /// Number of stops generated for the material gradient.
  unsigned int _gradient_n_stops{0u};

  /// Style of the grid overlaid on material surfaces.
  grid_style _grid_style{};
};

/// Style applied to an actsvg proto surface
struct surface_style {
  /// Construct a surface style with fill, stroke width, and material settings.
  /// @param fill_color Fill color for the surface and its outline.
  /// @param stroke_width Width of the surface outline.
  /// @param mat_style Material annotation and gradient settings.
  explicit surface_style(
      const actsvg::style::color& fill_color =
          colors::tableau_colorblind10::red_tones(0.4f).front(),
      actsvg::scalar stroke_width = 0.4f, surface_material_style mat_style = {})
      : _fill_color{fill_color},
        _stroke_width{stroke_width},
        _material_style{std::move(mat_style)} {
    _stroke_color = fill_color;
    _highlight_stroke_rgb = _highlight_rgb;
    _highlight_stroke_width = stroke_width;
    _material_style._gradient_font._size = 26u;
    _material_style._gradient_label._font._size = 26u;
  }
  /// Surface fill color.
  actsvg::style::color _fill_color{};

  /// Width of the surface outline.
  actsvg::scalar _stroke_width{10.f};
  /// Color of the surface outline.
  actsvg::style::color _stroke_color{};

  /// RGB color used to highlight the surface fill.
  std::array<int, 3> _highlight_rgb{colors::macaroni_and_cheese};
  /// Names of surface features that should be highlighted.
  std::vector<std::string> _highlights{};
  /// RGB color used to highlight the surface outline.
  std::array<int, 3> _highlight_stroke_rgb{};
  /// Width of the highlighted surface outline.
  actsvg::scalar _highlight_stroke_width{};

  /// Font size used for surface labels.
  unsigned int _font_size{22u};
  /// Number of segments used to approximate curved surface edges.
  unsigned int _n_segments{72u};

  /// Material annotation and gradient style for the surface.
  surface_material_style _material_style{};
};

/// Style applied to an actsvg proto portal link
struct link_style {
  /// Size of the marker at a portal link endpoint.
  actsvg::scalar _marker_size{1.f};
};

/// Style applied to an actsvg proto portal
struct portal_style {
  /// Style applied to the portal surface.
  surface_style _surface_style{};
  /// Style applied to links attached to the portal.
  link_style _link_style{};
};

/// Style applied to an actsvg proto volume
struct volume_style {
  /// Style for sensitive surfaces in the volume.
  surface_style _sensitive_surface_style{};
  /// Style for passive surfaces in the volume.
  surface_style _passive_surface_style{};
  /// Style for portals in the volume.
  portal_style _portal_style{};
  /// Style for the volume surface grid.
  grid_style _grid_style{};
};

/// Style applied to an actsvg proto detector
struct detector_style {
  /// Style applied to volumes in the detector.
  volume_style _volume_style{};
  /// Style applied to the detector grid.
  grid_style _grid_style{};
};

/// Style applied to proto eta lines
struct eta_lines_style {
  /// Color of the main eta lines.
  actsvg::style::color _fill_color_main{};
  /// Color of the half eta lines.
  actsvg::style::color _fill_color_half{};
  /// Width of the main eta lines.
  actsvg::scalar _stroke_width_main{10.f};
  /// Width of the half eta lines.
  actsvg::scalar _stroke_width_half{10.f};
  /// Font size of eta line labels.
  unsigned int _font_size{11u};
};

/// Style applied to a proto landmark
struct landmark_style {
  /// Fill color for landmark markers.
  actsvg::style::color _fill_color{};
  /// Outline width for landmark markers.
  actsvg::scalar _stroke_width{10.f};
  /// Size of landmark markers.
  actsvg::scalar _marker_size{10.f};
  /// Marker glyph used to draw landmarks.
  std::string _marker_type{};
};

/// Style applied to a proto trajectory
struct trajectory_style {
  /// Color used to distinguish tracks.
  actsvg::style::color _fill_color{};
  /// Width of trajectory lines.
  actsvg::scalar _stroke_width{10.f};
};

/// Global styling options
struct style {
  /// Styling for detector geometry.
  detector_style _detector_style{};
  /// Styling for eta guide lines.
  eta_lines_style _eta_lines_style{};
  /// Styling for trajectories.
  trajectory_style _trajectory_style{};
  /// Styling for track landmarks.
  landmark_style _landmark_style{};
  /// Styling for intersections.
  landmark_style _intersection_style{};
};

/// The default styling
namespace svg_default {

// Material style
/// Default material annotation and gradient style.
const styling::surface_material_style material_style{
    {-1000.f, 1000.f},
    {1000.f, -200.f},
    {20.f, 400.f},
    actsvg::style::font{},
    actsvg::style::label{"t/X0", actsvg::style::label::horizontal::right,
                         actsvg::style::label::vertical::top},
    {colors::black, 1.f},
    1.f,
    colors::gradient::rainbow_scale,
    4u,
    {{colors::black, 1.f}, 0.4f}};

// Surface styles
/// Default style for sensitive detector surfaces.
const styling::surface_style surface_style_sensitive{
    colors::blue_theme(0.5f).front(), 3.f, material_style};
/// Default style for passive detector surfaces.
const styling::surface_style surface_style_passive{
    {colors::dark_grey, 0.8f}, 3.f, material_style};
/// Default style for portal surfaces.
const styling::surface_style surface_style_portal{
    colors::red_theme(0.5f).front(), 3.f, material_style};

// Portal style
/// Default style for links attached to portals.
const styling::link_style link_style{1.2f};
/// Default portal style combining portal surfaces and links.
const styling::portal_style portal_style{surface_style_portal, link_style};

// Volume style
/// Default style for detector and volume grids.
const styling::grid_style grid_style{colors::red_theme(1.f).front(), 1.f};
/// Default volume style.
const styling::volume_style volume_style{
    surface_style_sensitive, surface_style_passive, portal_style, grid_style};

// Detector style
/// Default detector style.
const styling::detector_style detector_style{volume_style, grid_style};

// Eta lines style
/// Default eta guide line style.
const styling::eta_lines_style eta_lines_style{
    {colors::black, 1.f}, {colors::black, 1.f}, 1.f, 1.f, 30u};

// Intersection points and track state style
/// Default style for trajectory intersections.
const styling::landmark_style intersection_style{
    colors::black_theme(1.f).front(), 0.8f, 5.f, "x"};
/// Default style for trajectory landmarks.
const styling::landmark_style landmark_style{colors::green_theme(1.f).front(),
                                             0.8f, 3.f, "o"};

// Trajectory style
/// Default trajectory style.
const styling::trajectory_style trajectory_style{
    colors::green_theme(1.f).front(), 1.f};

// Full style
/// Complete default SVG style configuration.
const styling::style style{detector_style, eta_lines_style, trajectory_style,
                           landmark_style, intersection_style};
}  // namespace svg_default

/// Styling that matches data plotting
namespace tableau_colorblind {

// Surface material style
/// Material annotation and gradient style using the viridis scale.
const styling::surface_material_style material_style{
    {-1000.f, 1000.f},
    {1000.f, -200.f},
    {20.f, 400.f},
    actsvg::style::font{},
    actsvg::style::label{"t/X0", actsvg::style::label::horizontal::right,
                         actsvg::style::label::vertical::top},
    {colors::black, 1.f},
    1.f,
    colors::gradient::viridis_scale,
    8u,
    {{colors::black, 1.f}, 0.4f}};

// Surface styles
/// Passive surface style in the Tableau colorblind palette.
const styling::surface_style surface_style_passive{
    colors::tableau_colorblind10::grey_tones(0.8f).front(), 1.f,
    material_style};
/// Portal surface style in the Tableau colorblind palette.
const styling::surface_style surface_style_portal{
    colors::tableau_colorblind10::blue_tones(0.2f).front(), 1.5f,
    material_style};
/// Sensitive surface style in the Tableau colorblind palette.
const styling::surface_style surface_style_sensitive{
    colors::tableau_colorblind10::red_tones(0.3f).front(), 1.f, material_style};

// Portal style
/// Portal style using the Tableau portal surface and default link style.
const styling::portal_style portal_style{surface_style_portal,
                                         svg_default::link_style};

// Volume style
/// Grid style used by the Tableau detector and volumes.
const styling::grid_style grid_style{{colors::red, 1.f}, 0.6f};
/// Volume style using the Tableau colorblind palette.
const styling::volume_style volume_style{
    surface_style_sensitive, surface_style_passive, portal_style, grid_style};

// Detector style
/// Detector style using the Tableau colorblind palette.
const styling::detector_style detector_style{volume_style, grid_style};

// Full style
/// Complete style configuration using the Tableau colorblind palette.
const styling::style style{
    detector_style, svg_default::eta_lines_style, svg_default::trajectory_style,
    svg_default::landmark_style, svg_default::intersection_style};
}  // namespace tableau_colorblind

/// @brief Sets the style of the proto surface.
template <typename point3_container_t>
inline void apply_style(actsvg::proto::surface<point3_container_t>& p_surface,
                        const surface_style& styling) {
  // Fill color
  p_surface._fill = actsvg::style::fill(styling._fill_color);
  p_surface._fill._fc._hl_rgb = styling._highlight_rgb;
  p_surface._fill._fc._highlight = styling._highlights;

  // Stroke
  p_surface._stroke =
      actsvg::style::stroke(styling._stroke_color, styling._stroke_width);
  p_surface._stroke._sc._hl_rgb = styling._highlight_stroke_rgb;
  p_surface._stroke._hl_width = styling._highlight_stroke_width;
}

/// @brief Sets the style of the proto link.
template <typename point3_container_t>
inline void apply_style(
    typename actsvg::proto::portal<point3_container_t>::link& p_link,
    const link_style& styling) {
  p_link._end_marker._size = styling._marker_size;
}

/// @brief Sets the style of the proto portal.
template <typename point3_container_t>
inline void apply_style(actsvg::proto::portal<point3_container_t>& p_portal,
                        const portal_style& styling) {
  apply_style(p_portal._surface, styling._surface_style);

  for (auto& volume_link : p_portal._volume_links) {
    apply_style<point3_container_t>(volume_link, styling._link_style);
  }
}

/// @brief Sets the style of the proto grid.
inline void apply_style(actsvg::proto::grid& p_grid,
                        const grid_style& styling) {
  p_grid._stroke._sc = styling._stroke_color;
  p_grid._stroke._width = styling._stroke_width;
  // Use dashed lines for the grid
  // p_grid._stroke._dasharray = {4};
}

/// @brief Sets the style of the surface proto material.
inline void apply_style(actsvg::proto::surface_material& p_material,
                        const surface_material_style& styling) {
  p_material._info_pos = styling._info_pos;

  p_material._gradient_pos = styling._gradient_pos;
  p_material._gradient_box = styling._gradient_box;
  p_material._gradient_font = styling._gradient_font;
  p_material._gradient_label = styling._gradient_label;
  p_material._gradient_stroke._sc = styling._gradient_stroke_color;
  p_material._gradient_stroke._width = styling._gradient_stroke_width;
  p_material._gradient._stops = colors::gradient::generate_stops(
      styling._gradient_color_scale, styling._gradient_n_stops);

  apply_style(p_material._grid, styling._grid_style);
}

/// @brief Sets the style of the proto volume.
template <typename point3_container_t>
inline void apply_style(actsvg::proto::volume<point3_container_t>& p_volume,
                        const volume_style& styling) {
  using proto_sf_t = actsvg::proto::surface<point3_container_t>;

  for (auto& p_surface : p_volume._v_surfaces) {
    const auto& sf_style =
        (p_surface._sf_type == proto_sf_t::sf_type::e_sensitive)
            ? styling._sensitive_surface_style
            : styling._passive_surface_style;
    apply_style(p_surface, sf_style);
  }
  for (auto& p_portals : p_volume._portals) {
    apply_style(p_portals, styling._portal_style);
  }
  apply_style(p_volume._surface_grid, styling._grid_style);
}

/// @brief Sets the style of the proto volume.
template <typename point3_container_t>
inline void apply_style(actsvg::proto::detector<point3_container_t>& p_detector,
                        const detector_style& styling) {
  for (auto& p_volume : p_detector._volumes) {
    apply_style(p_volume, styling._volume_style);
  }
}

/// @brief Sets the style of the proto eta_lines.
inline void apply_style(meta::proto::eta_lines& p_eta_lines,
                        const eta_lines_style& styling) {
  p_eta_lines._stroke_main = actsvg::style::stroke(styling._fill_color_main,
                                                   styling._stroke_width_main);

  p_eta_lines._stroke_half = actsvg::style::stroke(styling._fill_color_half,
                                                   styling._stroke_width_half);
  p_eta_lines._stroke_half._dasharray = {2, 2};
  p_eta_lines._label_font._size = styling._font_size;
}

/// @brief Sets the style of the proto landmark.
template <concepts::point3D point3_t>
inline void apply_style(meta::proto::landmark<point3_t>& p_landmark,
                        const landmark_style& styling) {
  const auto fill = actsvg::style::fill(styling._fill_color);
  const auto stroke =
      actsvg::style::stroke(styling._fill_color, styling._stroke_width);

  const auto marker = actsvg::style::marker{styling._marker_type,
                                            styling._marker_size, fill, stroke};
  p_landmark._marker = marker;
}

/// @brief Sets the style of the proto intersection record.
template <concepts::point3D point3_t>
inline void apply_style(meta::proto::intersection<point3_t>& p_intersection,
                        const landmark_style& styling) {
  const auto stroke =
      actsvg::style::stroke(styling._fill_color, styling._stroke_width);

  auto fill_color = styling._fill_color;
  fill_color._opacity *= 0.5f;
  const auto fill = actsvg::style::fill(fill_color);

  const auto marker = actsvg::style::marker{styling._marker_type,
                                            styling._marker_size, fill, stroke};

  for (auto& p_landmark : p_intersection._landmarks) {
    p_landmark._marker = marker;
  }
}

/// @brief Sets the style of the proto trajectory.
template <concepts::point3D point3_t>
inline void apply_style(meta::proto::trajectory<point3_t>& p_trajectory,
                        const trajectory_style& styling) {
  p_trajectory._stroke =
      actsvg::style::stroke(styling._fill_color, styling._stroke_width);
}

}  // namespace detray::svgtools::styling
