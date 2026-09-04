// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/grid_axis.hpp"
#include "detray/io/frontend/definitions.hpp"
#include "detray/utils/invalid_values.hpp"

// System include(s)
#include <array>
#include <cstdint>
#include <iomanip>
#include <limits>
#include <map>
#include <optional>
#include <ostream>
#include <string>
#include <vector>

//
// Intermediate data representation as an abstraction to different IO formats
//

/// Raw indices (std::size_t) denote links between data components in different
/// files, while links used in detray detector objects are modelled as e.g.
/// @c single_link_payload
///
/// @note All indices are with respect to volume that the data belongs to,
/// starting at zero
namespace detray::io {

/// @brief a payload for common information
struct common_header_payload {
  /// Detray version number
  std::string version{};
  /// Detector name
  std::string detector{};
  /// Type of data the file contains (e.g. geometry vs material_maps)
  std::string tag{};
  /// Date of file creation
  std::string date{};
};

/// @brief a payload for common and extra information
template <typename sub_header_payload_t = bool>
struct header_payload {
  common_header_payload common{};
  /// Information that is specific to a tag, i.e. type of detector payload
  std::optional<sub_header_payload_t> sub_header;
};

/// @brief A payload for a single object link (volume local)
struct single_link_payload {
  std::size_t link{std::numeric_limits<std::size_t>::max()};
};

/// @brief A payload for a typed object link (@c dtyped_index)
template <typename type_id_t>
struct typed_link_payload {
  using type_id = type_id_t;

  /// Type id of the object, e.g. @c shape_id
  type_id type{type_id::unknown};
  /// Index of the object within its volume
  std::size_t index{std::numeric_limits<std::size_t>::max()};
};

/// Geometry payloads
/// @{

/// @brief a payload for the geometry specific part of the file header
struct geo_sub_header_payload {
  /// Number of volumes in the detector
  std::size_t n_volumes{0ul};
  /// Total number of surfaces in the detector
  std::size_t n_surfaces{0ul};
};

/// @brief a payload for the geometry file header
using geo_header_payload = header_payload<geo_sub_header_payload>;

/// @brief A payload object to link a surface to its material
using material_link_payload = typed_link_payload<io::material_id>;

/// @brief A payload object to link a volume to its acceleration data structures
using acc_links_payload = typed_link_payload<io::accel_id>;

/// @brief A payload for an affine transformation in homogeneous coordinates
struct transform_payload {
  /// Translation
  std::array<io::scalar, 3u> tr{};
  // Column major rotation matrix
  std::array<io::scalar, 9u> rot{};

 private:
  DETRAY_HOST friend inline std::ostream& operator<<(
      std::ostream& os, const transform_payload& transform) {
    const auto& r = transform.rot;
    os << "rot: ";
    os << std::fixed << std::setw(4);
    auto print_line = [&](std::size_t i) {
      os << r[i] << " " << r[i + 3] << " " << r[i + 6] << "\n";
    };
    print_line(0);
    os << "     ";
    print_line(1);
    os << "     ";
    print_line(2);

    const auto& t = transform.tr;

    os << "tr:  ";
    os << t[0] << " " << t[1] << " " << t[2];
    return os;
  }
};

/// @brief A payload object for surface masks
struct mask_payload {
  using mask_shape = io::shape_id;

  mask_shape shape{mask_shape::unknown};
  /// Sensitive/passive surface: Mother volume index, Portal: Volume the
  /// portal leads into
  single_link_payload volume_link{};
  /// Boundary values according to @c boundaries enum of the shape in @c shape
  std::vector<io::scalar> boundaries{};
};

/// @brief  A payload for surfaces
struct surface_payload {
  /// Position of the surface in the volume surface container.
  /// Optional, if no specific ordering
  std::optional<std::size_t> index_in_coll;
  transform_payload transform{};
  std::vector<mask_payload> masks{};
  /// Write the material link of the surface as additional information
  std::optional<material_link_payload> material;
  /// Identifier of the source object (e.g. value of ACTS::GeometryIdentifier)
  std::uint64_t source{};
  /// Write the surface identifier as an additional information
  std::optional<std::uint64_t> identifier{
      detray::detail::invalid_value<std::uint64_t>()};
  /// Either sensitive, passive or portal
  detray::surface_id type{detray::surface_id::e_sensitive};
};

/// @brief A payload for volumes
struct volume_payload {
  std::string name{};
  detray::volume_id type{detray::volume_id::e_cylinder};
  transform_payload transform{};
  std::vector<surface_payload> surfaces{};
  /// Index of the volume in the detector volume container
  single_link_payload index{};
  /// Optional accelerator data structures
  std::optional<std::vector<acc_links_payload>> acc_links{};
};

/// @}

/// Material payloads
/// @{

/// @brief a payload for the material specific part of the file header
struct homogeneous_material_sub_header_payload {
  /// Total number of material slabs in the detector
  std::size_t n_slabs{0ul};
  /// Total number of surface with material slabs
  std::size_t n_slab_surfaces{0ul};
  /// Total number of material rods (wire material) in the detector
  std::size_t n_rods{0ul};
  /// Total number of surface with material rods
  std::size_t n_rod_surfaces{0ul};
};

/// @brief a payload for the homogeneous material file header
using homogeneous_material_header_payload =
    header_payload<homogeneous_material_sub_header_payload>;

/// @brief A payload object for a material parametrization
struct material_param_payload {
  /// Material parametrization
  /// 0: X0
  /// 1: L0
  /// 2: Ar
  /// 3: Z
  /// 4: Mass density
  /// 5: Molar densitty (unused)
  /// 6: @c material_state enum (solid, liquid, gaseous)
  std::array<io::scalar, 7u> params{};
};

/// @brief A payload object for a material slab/rod
struct surface_material_payload {
  using mat_type = io::material_id;

  /// Either 'slab' or 'rod'
  mat_type type{mat_type::unknown};
  /// Index of the material in the material container
  /// Optional, if no specific sorting required
  std::optional<std::size_t> index_in_coll;
  /// Index of the surface in the volume the material belongs to
  single_link_payload surface{};
  /// Thickness of the material slab
  io::scalar thickness{std::numeric_limits<io::scalar>::max()};
  /// Material parameters
  material_param_payload mat{};
};

/// @brief A payload for the material of a volume
struct volume_material_payload {
  // Not implemented
};

/// @brief A payload object for all of the [surface] material contained in a
/// volume
struct material_volume_payload {
  /// Volume index the payload belongs to
  single_link_payload volume_link{};
  /// Material of the volume
  std::optional<volume_material_payload> volume_mat;
  /// Material slabs/rods of contained surfaces
  std::vector<surface_material_payload> surface_mat{};
};

/// @brief A payload for the homogeneous material description of a detector
struct detector_homogeneous_material_payload {
  std::vector<material_volume_payload> volumes{};
};

/// @}

/// Payloads for a uniform grid
/// @{

/// @brief a payload for the grid specific part of the file header
struct grid_sub_header_payload {
  std::size_t n_grids{0u};
};

/// @brief a payload for the grid file header
using grid_header_payload = header_payload<grid_sub_header_payload>;

/// @brief axis definition and bin edges
struct axis_payload {
  /// Axis binning type
  axis::binning binning{axis::binning::e_regular};
  /// Axis behaviour at final bin
  axis::bounds bounds{axis::bounds::e_closed};
  /// Label to retrieve the axis
  axis::label label{axis::label::e_r};

  /// Number of bins
  std::size_t bins{0u};
  /// Axis span for reg. binning or distinct bin edge values for irr. binning
  std::vector<io::scalar> edges{};

 private:
  DETRAY_HOST friend inline std::ostream& operator<<(std::ostream& os,
                                                     const axis_payload& axis) {
    os << "axis_payload{binning: " << axis.binning
       << ", bounds: " << axis.bounds << ", label: " << axis.label
       << ", bins: " << axis.bins << ", edges: [";
    for (std::size_t i = 0; i < axis.edges.size(); ++i) {
      if (i > 0) {
        os << ", ";
      }
      os << axis.edges[i];
    }
    os << "]}";
    return os;
  }
};

/// @brief A payload for a single grid bin
/// @note The local indices and bin content container need to have the same
/// ordering for the data to be matched
template <typename content_t = std::size_t>
struct grid_bin_payload {
  /// Bin index per axis (global index is determined by serializer)
  std::vector<unsigned int> loc_index{};
  /// Bin content for the bin (e.g. coll. of surface indices or material slab)
  std::vector<content_t> content{};
};

/// @brief A payload for a grid definition
/// @tparam bin_content_t the value type of the grid, e.g. surface indices or
/// surface_material_payloads
/// @tparam grid_id_t grid type enum, e.g. @c accel_id or @c material_id
template <typename bin_content_t = std::size_t,
          typename grid_id_t = io::accel_id>
struct grid_payload {
  using grid_type = grid_id_t;
  // volume-local surface or global volume index the grids belongs to
  single_link_payload owner_link{};
  // Type and global index of the grid (index is optional only for debugging)
  typed_link_payload<grid_id_t> grid_link{};

  /// Must match the number of axis in the grid type indicated by @c grid_link
  std::vector<axis_payload> axes{};
  std::vector<grid_bin_payload<bin_content_t>> bins{};
};

/// @brief A payload for the grid collections of a detector
template <typename bin_content_t = std::size_t,
          typename grid_id_t = io::accel_id>
struct detector_grids_payload {
  // A collection of grids (can be the surface search grid, volume or surface
  // material grids) per volume index (not all volumes may contain grids)
  std::map<std::size_t, std::vector<grid_payload<bin_content_t, grid_id_t>>>
      grids = {};
};

/// @}

/// @brief A payload for a detector geometry
struct detector_payload {
  std::vector<volume_payload> volumes = {};
  /// Volume acceleration structure
  std::optional<grid_payload<std::size_t, io::accel_id>> volume_grid;
};

}  // namespace detray::io
