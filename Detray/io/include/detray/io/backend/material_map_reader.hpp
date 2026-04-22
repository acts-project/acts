// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2023-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Project include(s)
#include "detray/builders/detector_builder.hpp"
#include "detray/builders/material_map_builder.hpp"
#include "detray/geometry/concepts.hpp"
#include "detray/io/backend/detail/basic_converter.hpp"
#include "detray/io/backend/detail/type_info.hpp"
#include "detray/io/backend/homogeneous_material_reader.hpp"
#include "detray/io/frontend/payloads.hpp"
#include "detray/material/material_slab.hpp"
#include "detray/utils/type_list.hpp"

// System include(s)
#include <memory>
#include <stdexcept>
#include <string_view>
#include <vector>

namespace detray::io {

/// @brief Material map reader backend
///
/// Fills a @c detector_builder from a material @c detector_grids_payload
template <typename DIM = std::integral_constant<std::size_t, 2u>>
class material_map_reader {
  using material_reader_t = homogeneous_material_reader;

 public:
  static constexpr std::size_t dim{DIM()};

  /// Payload type that the reader processes
  using payload_type =
      detector_grids_payload<surface_material_payload, io::material_id>;

  using bin_index_type = axis::multi_bin<dim>;

  /// Tag the reader as "material_maps"
  static constexpr std::string_view tag = "material_maps";

  /// Convert the material grids @param grids_data from their IO
  /// payload
  template <typename detector_t>
  static void from_payload(detector_builder<typename detector_t::metadata,
                                            volume_builder> &det_builder,
                           payload_type &&grids_data) {
    DETRAY_VERBOSE_HOST("Reading payload object...");

    using scalar_t = dscalar<typename detector_t::algebra_type>;
    using mat_factory_t = material_map_factory<detector_t, bin_index_type>;
    using mat_data_t = typename mat_factory_t::data_type;
    using mat_id = typename detector_t::material::id;

    // Map the io::shape_id in the payload to the shape types known by
    // the detector under construction
    using mat_registry_t =
        io::detail::filtered_material_map_registry<detector_t>;
    static_assert((!types::contains<mat_registry_t, io::detail::unknown_type> &&
                   types::size<mat_registry_t> > 0) ||
                  (types::contains<mat_registry_t, io::detail::unknown_type> &&
                   types::size<mat_registry_t> > 1));

    DETRAY_VERBOSE_HOST("Converting material grids for "
                        << grids_data.grids.size() << " volumes");
    // Convert the material volume by volume
    for (const auto &[vol_idx, mat_grids] : grids_data.grids) {
      if (!det_builder.has_volume(vol_idx)) {
        std::stringstream err_stream;
        err_stream << "Volume " << vol_idx << ": "
                   << "Cannot build material map for volume "
                   << "(volume not registered in detector builder)";
        DETRAY_FATAL_HOST(err_stream.str());
        throw std::invalid_argument(err_stream.str());
      }

      DETRAY_VERBOSE_HOST("Reading material maps in volume: "
                          << det_builder[static_cast<dindex>(vol_idx)]->name());

      // Decorate the current volume builder with material maps
      auto vm_builder =
          det_builder.template decorate<material_map_builder<detector_t, dim>>(
              static_cast<dindex>(vol_idx));

      DETRAY_VERBOSE_HOST("-> Found " << mat_grids.size()
                                      << " material grids:");

      // Add the material data to the factory
      auto mat_factory =
          std::make_shared<material_map_factory<detector_t, bin_index_type>>();

      // Convert the material grid of each surface
      for (const auto &[idx, grid_data] : detray::views::enumerate(mat_grids)) {
        DETRAY_VERBOSE_HOST("Reading material map payload #" << idx << "...");
        // Get the number of bins per axis
        std::vector<std::size_t> n_bins{};
        for (const auto &axis_data : grid_data.axes) {
          n_bins.push_back(axis_data.bins);
        }

        DETRAY_VERBOSE_HOST("-> Belongs to "
                            << (dim == 2 ? "surface" : "volume") << " = "
                            << grid_data.owner_link.link);

        mat_id map_id =
            from_payload<mat_registry_t, detector_t>(grid_data.grid_link.type);
        DETRAY_VERBOSE_HOST("-> Type id: " << map_id);

        DETRAY_VERBOSE_HOST("-> Reading axis bins: Dims = " << dim);
        for ([[maybe_unused]] const auto &[i, count] :
             detray::views::enumerate(n_bins)) {
          DETRAY_VERBOSE_HOST("--> Axis " << i << ": " << count << " bins");
        }
        // Get the axis spans
        DETRAY_VERBOSE_HOST("-> Reading axis spans:");
        std::vector<std::vector<scalar_t>> axis_spans = {};
        for (const auto &[i, axis_data] :
             detray::views::enumerate(grid_data.axes)) {
          axis_spans.push_back({static_cast<scalar_t>(axis_data.edges.front()),
                                static_cast<scalar_t>(axis_data.edges.back())});
          DETRAY_VERBOSE_HOST("--> Axis " << i << " ["
                                          << axis_spans.back().at(0) << ", "
                                          << axis_spans.back().at(1) << "]");
        }

        // Get the local bin indices and the material parametrization
        DETRAY_VERBOSE_HOST("-> Reading bin content...");

        std::vector<bin_index_type> loc_bins{};
        mat_data_t mat_data{
            detail::basic_converter::from_payload(grid_data.owner_link)};
        DETRAY_VERBOSE_HOST("--> Found " << grid_data.bins.size() << " bins");
        for (const auto &bin_data : grid_data.bins) {
          assert(dim == bin_data.loc_index.size() &&
                 "Dimension of local bin indices in input file does "
                 "not match material grid dimension");

          // The local bin indices for the bin to be filled
          bin_index_type mbin;
          for (const auto &[i, bin_idx] :
               detray::views::enumerate(bin_data.loc_index)) {
            mbin[i] = bin_idx;
          }
          loc_bins.push_back(std::move(mbin));

          // Add the material slab per bin
          for (const auto &slab_data : bin_data.content) {
            mat_data.append(
                material_reader_t::template from_payload<scalar_t>(slab_data));
          }
        }

        DETRAY_VERBOSE_HOST("Finished reading data for material map #" << idx);
        DETRAY_DEBUG_HOST(mat_data);

        DETRAY_VERBOSE_HOST("Add material data to factory...");
        mat_factory->add_material(map_id, std::move(mat_data),
                                  std::move(n_bins), std::move(axis_spans),
                                  std::move(loc_bins));
      }

      // Add the material maps to the volume
      DETRAY_VERBOSE_HOST("Add the material maps to the volume:");
      vm_builder->add_surfaces(mat_factory);
    }
  }

 private:
  /// Get the detector material id from the payload material type id
  template <typename mat_registry_t, typename detector_t, std::size_t I = 0u>
  static typename detector_t::material::id from_payload(
      io::material_id type_id) {
    // Get the next mask shape type
    using frame_t = types::at<mat_registry_t, I>;
    using algebra_t = typename detector_t::algebra_type;

    // Material id of map data found
    if constexpr (!std::is_same_v<frame_t, io::detail::unknown_type> &&
                  ::detray::concepts::coordinate_frame<frame_t>) {
      // Get the corresponding material id for this detector
      constexpr auto mat_id{
          types::id<io::material_registry<algebra_t>, frame_t>};
      if (type_id == mat_id) {
        using mat_frame_registry_t = io::detail::mat_frame_registry<detector_t>;

        constexpr std::size_t mapped_idx{
            types::position<mat_frame_registry_t, frame_t>};

        return types::id_cast<typename detector_t::material,
                              mat_frame_registry_t::original_index(mapped_idx)>;
      }
    }
    // Test next material type id
    if constexpr (I < types::size<mat_registry_t> - 1u) {
      return from_payload<mat_registry_t, detector_t, I + 1u>(type_id);
    } else {
      return detector_t::material::id::e_none;
    }
  }
};

}  // namespace detray::io
