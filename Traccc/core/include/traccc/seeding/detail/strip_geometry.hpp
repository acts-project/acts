/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

#include <cmath>
#include <cstdint>
#include <limits>

#include "traccc/definitions/qualifiers.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/seeding/strip_spacepoint_formation_data.hpp"

namespace traccc::details {

/// Direction spans the full strip; trajectory is twice center - beam_spot.
struct strip_material {
  point3 center{};
  vector3 direction{};
  vector3 trajectory{};
  vector3 normal{};
  scalar half_length{0.f};
  unsigned int valid{0u};
};

template <typename surface_info_collection_t>
TRACCC_HOST_DEVICE inline bool find_strip_surface_info(
    const surface_info_collection_t& surface_infos,
    const std::uint64_t surface_link, strip_measurement_surface_info& result) {
  unsigned int first = 0u;
  unsigned int last = static_cast<unsigned int>(surface_infos.size());
  while (first < last) {
    const unsigned int middle = first + (last - first) / 2u;
    const auto candidate = surface_infos.at(middle);
    if (candidate.surface_link < surface_link) {
      first = middle + 1u;
    } else {
      last = middle;
    }
  }
  if ((first >= surface_infos.size()) ||
      (surface_infos.at(first).surface_link != surface_link)) {
    return false;
  }
  result = surface_infos.at(first);
  return true;
}

struct exact_point2 {
  double u{0.};
  double v{0.};
};

struct exact_point3 {
  double x{0.};
  double y{0.};
  double z{0.};
};

TRACCC_HOST_DEVICE inline exact_point3 surface_global_position(
    const strip_measurement_surface_info& info, const exact_point2& p) {
  return {info.origin[0] + p.u * info.local_u[0] + p.v * info.local_v[0],
          info.origin[1] + p.u * info.local_u[1] + p.v * info.local_v[1],
          info.origin[2] + p.u * info.local_u[2] + p.v * info.local_v[2]};
}

template <typename measurement_backend_t>
TRACCC_HOST_DEVICE inline int measured_strip_index(
    const edm::measurement<measurement_backend_t>& measurement,
    const strip_index_mapping& mapping) {
  if ((mapping.component > 1u) || (mapping.count == 0u) ||
      !(mapping.pitch > 0.)) {
    return -1;
  }
  const double local =
      static_cast<double>(measurement.local_position()[mapping.component]);
  const double value =
      mapping.scale * std::floor(local / mapping.pitch) + mapping.offset;
  // Check before the integer conversion, including NaN and infinities.
  if (!(value >= static_cast<double>(std::numeric_limits<int>::lowest())) ||
      !(value <= static_cast<double>(std::numeric_limits<int>::max()))) {
    return -1;
  }
  if (mapping.check_before_truncation &&
      !((value >= 0.) && (value < static_cast<double>(mapping.count)))) {
    return -1;
  }
  const int index = static_cast<int>(value);
  return ((index >= 0) && (static_cast<unsigned int>(index) < mapping.count))
             ? index
             : -1;
}

struct local_strip_endpoints {
  exact_point2 first{};
  exact_point2 second{};
  unsigned int valid{0u};
};

TRACCC_HOST_DEVICE inline local_strip_endpoints linear_strip_endpoints(
    const int strip, const linear_strip_geometry& geometry) {
  local_strip_endpoints result{};
  if ((strip < 0) || (geometry.n_rows == 0u) ||
      (static_cast<unsigned int>(strip) >= geometry.n_strips) ||
      !(geometry.pitch > 0.) || !(geometry.row_length > 0.)) {
    return result;
  }
  int row = 0;
  if (geometry.n_rows > 1u) {
    const double row_value =
        std::floor(geometry.row_coordinate / geometry.row_length) +
        static_cast<double>(geometry.n_rows) * 0.5;
    if (!(row_value >=
          static_cast<double>(std::numeric_limits<int>::lowest())) ||
        !(row_value <= static_cast<double>(std::numeric_limits<int>::max()))) {
      return result;
    }
    row = static_cast<int>(row_value);
  }
  if ((row < 0) || (static_cast<unsigned int>(row) >= geometry.n_rows)) {
    return result;
  }
  const double start =
      (static_cast<double>(row) - static_cast<double>(geometry.n_rows) * 0.5) *
      geometry.row_length;
  const double end = start + geometry.row_length;
  const double transverse =
      (static_cast<double>(strip) -
       static_cast<double>(geometry.n_strips) * 0.5 + 0.5) *
      geometry.pitch;
  result.first = {start, transverse};
  result.second = {end, transverse};
  result.valid = 1u;
  return result;
}

TRACCC_HOST_DEVICE inline exact_point2 radial_strip_position_at_radius(
    const int strip, const double radius,
    const radial_strip_geometry& geometry) {
  const double phi = (static_cast<double>(strip) -
                      static_cast<double>(geometry.n_strips) * 0.5 + 0.5) *
                     geometry.angular_pitch;
  const double b =
      -2. * geometry.frame_chord * std::sin(0.5 * geometry.stereo_angle + phi);
  const double c =
      geometry.frame_chord * geometry.frame_chord - radius * radius;
  const double strip_r = 0.5 * (-b + std::sqrt(b * b - 4. * c));
  const double strip_x = strip_r * std::cos(phi);
  const double strip_y = strip_r * std::sin(phi);
  const double x = geometry.cos_stereo * (strip_x - geometry.focal_radius) -
                   geometry.sin_stereo * strip_y + geometry.focal_radius;
  const double y = geometry.sin_stereo * (strip_x - geometry.focal_radius) +
                   geometry.cos_stereo * strip_y;
  if (geometry.local_frame == strip_local_frame::polar) {
    return {std::sqrt(x * x + y * y), std::atan2(y, x)};
  }
  return {x, y};
}

template <typename measurement_backend_t>
TRACCC_HOST_DEVICE inline strip_material make_strip_material(
    const edm::measurement<measurement_backend_t>& measurement,
    const strip_measurement_surface_info& info, const point3& beam_spot) {
  strip_material result{};
  if (measurement.dimensions() != 1u) {
    return result;
  }
  const int strip = measured_strip_index(measurement, info.index_mapping);
  if (strip < 0) {
    return result;
  }
  exact_point3 first{};
  exact_point3 second{};
  switch (info.geometry_model) {
    case strip_geometry_model::linear: {
      const auto endpoints = linear_strip_endpoints(strip, info.linear);
      if (endpoints.valid == 0u) {
        return result;
      }
      first = surface_global_position(info, endpoints.first);
      second = surface_global_position(info, endpoints.second);
      break;
    }
    case strip_geometry_model::radial: {
      if ((static_cast<unsigned int>(strip) >= info.radial.n_strips) ||
          !(info.radial.angular_pitch > 0.) || !(info.radial.min_radius > 0.) ||
          !(info.radial.max_radius > info.radial.min_radius)) {
        return result;
      }
      first = surface_global_position(
          info, radial_strip_position_at_radius(strip, info.radial.min_radius,
                                                info.radial));
      second = surface_global_position(
          info, radial_strip_position_at_radius(strip, info.radial.max_radius,
                                                info.radial));
      break;
    }
    default:
      return result;
  }

  // Keep endpoint arithmetic in double precision until the completed
  // material is converted to the tracking scalar type.
  const double center_x = 0.5 * (first.x + second.x);
  const double center_y = 0.5 * (first.y + second.y);
  const double center_z = 0.5 * (first.z + second.z);
  const double direction_x = first.x - second.x;
  const double direction_y = first.y - second.y;
  const double direction_z = first.z - second.z;
  const double trajectory_x =
      2. * (center_x - static_cast<double>(beam_spot[0]));
  const double trajectory_y =
      2. * (center_y - static_cast<double>(beam_spot[1]));
  const double trajectory_z =
      2. * (center_z - static_cast<double>(beam_spot[2]));
  const double normal_x =
      direction_y * trajectory_z - direction_z * trajectory_y;
  const double normal_y =
      direction_z * trajectory_x - direction_x * trajectory_z;
  const double normal_z =
      direction_x * trajectory_y - direction_y * trajectory_x;
  const double length =
      std::sqrt(direction_x * direction_x + direction_y * direction_y +
                direction_z * direction_z);

  result.center = {static_cast<scalar>(center_x), static_cast<scalar>(center_y),
                   static_cast<scalar>(center_z)};
  result.direction = {static_cast<scalar>(direction_x),
                      static_cast<scalar>(direction_y),
                      static_cast<scalar>(direction_z)};
  result.trajectory = {static_cast<scalar>(trajectory_x),
                       static_cast<scalar>(trajectory_y),
                       static_cast<scalar>(trajectory_z)};
  result.normal = {static_cast<scalar>(normal_x), static_cast<scalar>(normal_y),
                   static_cast<scalar>(normal_z)};
  result.half_length = static_cast<scalar>(0.5 * length);
  result.valid = length > 0. ? 1u : 0u;
  return result;
}

}  // namespace traccc::details
