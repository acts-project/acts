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

template <typename collection_t, typename value_t, typename comparison_t>
TRACCC_HOST_DEVICE inline unsigned int collection_bound(
    const collection_t& collection, const value_t& value,
    comparison_t comparison) {
  unsigned int first = 0u;
  unsigned int last = static_cast<unsigned int>(collection.size());
  while (first < last) {
    const unsigned int middle = first + (last - first) / 2u;
    if (comparison(collection.at(middle), value)) {
      first = middle + 1u;
    } else {
      last = middle;
    }
  }
  return first;
}

struct strip_surface_before {
  TRACCC_HOST_DEVICE bool operator()(
      const strip_measurement_surface_info& candidate,
      const std::uint64_t surface_link) const {
    return candidate.surface_link < surface_link;
  }
};

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
  const unsigned int first =
      collection_bound(surface_infos, surface_link, strip_surface_before{});
  if ((first >= surface_infos.size()) ||
      (surface_infos.at(first).surface_link != surface_link)) {
    return false;
  }
  result = surface_infos.at(first);
  return true;
}

TRACCC_HOST_DEVICE inline point3 surface_global_position(
    const strip_measurement_surface_info& info, const point2& p) {
  return {info.origin[0] + p[0] * info.local_u[0] + p[1] * info.local_v[0],
          info.origin[1] + p[0] * info.local_u[1] + p[1] * info.local_v[1],
          info.origin[2] + p[0] * info.local_u[2] + p[1] * info.local_v[2]};
}

template <typename measurement_backend_t>
TRACCC_HOST_DEVICE inline int measured_strip_index(
    const edm::measurement<measurement_backend_t>& measurement,
    const strip_index_mapping& mapping) {
  if ((mapping.component > 1u) || (mapping.count == 0u) ||
      !(mapping.pitch > 0.f)) {
    return -1;
  }
  const scalar local = measurement.local_position()[mapping.component];
  const scalar value =
      mapping.scale * std::floor(local / mapping.pitch) + mapping.offset;
  // Check before the integer conversion, including NaN and infinities.
  if (!(value >= static_cast<scalar>(std::numeric_limits<int>::lowest())) ||
      !(value <= static_cast<scalar>(std::numeric_limits<int>::max()))) {
    return -1;
  }
  if (mapping.check_before_truncation &&
      !((value >= 0.f) && (value < static_cast<scalar>(mapping.count)))) {
    return -1;
  }
  const int index = static_cast<int>(value);
  return ((index >= 0) && (static_cast<unsigned int>(index) < mapping.count))
             ? index
             : -1;
}

struct local_strip_endpoints {
  point2 first{};
  point2 second{};
  unsigned int valid{0u};
};

TRACCC_HOST_DEVICE inline local_strip_endpoints linear_strip_endpoints(
    const int strip, const linear_strip_geometry& geometry) {
  local_strip_endpoints result{};
  if ((strip < 0) || (geometry.n_rows == 0u) ||
      (static_cast<unsigned int>(strip) >= geometry.n_strips) ||
      !(geometry.pitch > 0.f) || !(geometry.row_length > 0.f)) {
    return result;
  }
  int row = 0;
  if (geometry.n_rows > 1u) {
    const scalar row_value =
        std::floor(geometry.row_coordinate / geometry.row_length) +
        static_cast<scalar>(geometry.n_rows) * 0.5f;
    if (!(row_value >=
          static_cast<scalar>(std::numeric_limits<int>::lowest())) ||
        !(row_value <= static_cast<scalar>(std::numeric_limits<int>::max()))) {
      return result;
    }
    row = static_cast<int>(row_value);
  }
  if ((row < 0) || (static_cast<unsigned int>(row) >= geometry.n_rows)) {
    return result;
  }
  const scalar start =
      (static_cast<scalar>(row) - static_cast<scalar>(geometry.n_rows) * 0.5f) *
      geometry.row_length;
  const scalar end = start + geometry.row_length;
  const scalar transverse =
      (static_cast<scalar>(strip) -
       static_cast<scalar>(geometry.n_strips) * 0.5f + 0.5f) *
      geometry.pitch;
  result.first = {start, transverse};
  result.second = {end, transverse};
  result.valid = 1u;
  return result;
}

TRACCC_HOST_DEVICE inline point2 radial_strip_position_at_radius(
    const int strip, const scalar radius,
    const radial_strip_geometry& geometry) {
  const scalar phi = (static_cast<scalar>(strip) -
                      static_cast<scalar>(geometry.n_strips) * 0.5f + 0.5f) *
                     geometry.angular_pitch;
  const scalar b = -2.f * geometry.frame_chord *
                   std::sin(0.5f * geometry.stereo_angle + phi);
  const scalar c =
      geometry.frame_chord * geometry.frame_chord - radius * radius;
  const scalar strip_r = 0.5f * (-b + std::sqrt(b * b - 4.f * c));
  const scalar strip_x = strip_r * std::cos(phi);
  const scalar strip_y = strip_r * std::sin(phi);
  const scalar x = geometry.cos_stereo * (strip_x - geometry.focal_radius) -
                   geometry.sin_stereo * strip_y + geometry.focal_radius;
  const scalar y = geometry.sin_stereo * (strip_x - geometry.focal_radius) +
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
  point3 first{};
  point3 second{};
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
          !(info.radial.angular_pitch > 0.f) ||
          !(info.radial.min_radius > 0.f) ||
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

  const scalar center_x = 0.5f * (first[0] + second[0]);
  const scalar center_y = 0.5f * (first[1] + second[1]);
  const scalar center_z = 0.5f * (first[2] + second[2]);
  const scalar direction_x = first[0] - second[0];
  const scalar direction_y = first[1] - second[1];
  const scalar direction_z = first[2] - second[2];
  const scalar trajectory_x = 2.f * (center_x - beam_spot[0]);
  const scalar trajectory_y = 2.f * (center_y - beam_spot[1]);
  const scalar trajectory_z = 2.f * (center_z - beam_spot[2]);
  const scalar normal_x =
      direction_y * trajectory_z - direction_z * trajectory_y;
  const scalar normal_y =
      direction_z * trajectory_x - direction_x * trajectory_z;
  const scalar normal_z =
      direction_x * trajectory_y - direction_y * trajectory_x;
  const scalar length =
      std::sqrt(direction_x * direction_x + direction_y * direction_y +
                direction_z * direction_z);

  result.center = {center_x, center_y, center_z};
  result.direction = {direction_x, direction_y, direction_z};
  result.trajectory = {trajectory_x, trajectory_y, trajectory_z};
  result.normal = {normal_x, normal_y, normal_z};
  result.half_length = 0.5f * length;
  result.valid = length > 0.f ? 1u : 0u;
  return result;
}

}  // namespace traccc::details
