// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/builders/detector_builder.hpp"
#include "detray/builders/homogeneous_material_builder.hpp"
#include "detray/builders/homogeneous_material_factory.hpp"
#include "detray/core/concepts.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/io/backend/detail/basic_converter.hpp"
#include "detray/io/backend/detail/type_info.hpp"
#include "detray/io/frontend/payloads.hpp"
#include "detray/material/material.hpp"

// System include(s)
#include <memory>
#include <string_view>
#include <vector>

namespace detray::io {

/// @brief Homogeneous material reader backend
///
/// Fills a @c detector_builder from a @c detector_homogeneous_material_payload
class homogeneous_material_reader {
  /// IO material ids do not need to coincide with the detector ids,
  /// they are shared with ACTS
  using material_type = io::material_id;

 public:
  /// Tag the reader as "homogeneous material"
  static constexpr std::string_view tag = "homogeneous_material";

  /// Payload type that the reader processes
  using payload_type = detector_homogeneous_material_payload;

  /// Convert the detector material @param det_mat_data from its IO
  /// payload
  template <class detector_t>
  static void from_payload(detector_builder<typename detector_t::metadata,
                                            volume_builder>& det_builder,
                           const payload_type& det_mat_data) {
    DETRAY_VERBOSE_HOST("Reading payload object...");

    using scalar_t = dscalar<typename detector_t::algebra_type>;
    using mat_id = typename detector_t::material::id;

    DETRAY_DEBUG_HOST("Converting material for " << det_mat_data.volumes.size()
                                                 << " volumes");

    // Convert the material volume by volume
    for (const auto& mv_data : det_mat_data.volumes) {
      const auto vol_idx{
          detail::basic_converter::from_payload(mv_data.volume_link)};

      DETRAY_DEBUG_HOST(" - volume index from payload is " << vol_idx);

      if (!det_builder.has_volume(vol_idx)) {
        std::stringstream err_stream;
        err_stream << "Volume " << vol_idx << ": "
                   << "Cannot build homogeneous material for volume "
                   << "(volume not registered in detector builder)";
        DETRAY_FATAL_HOST(err_stream.str());
        throw std::invalid_argument(err_stream.str());
      }

      // Decorate the current volume builder with material
      auto vm_builder =
          det_builder
              .template decorate<homogeneous_material_builder<detector_t>>(
                  vol_idx);
      DETRAY_DEBUG_HOST(
          "Retrieved volume builder from detector build using vol_idx="
          << vol_idx);

      // Add the material data to the factory
      auto mat_factory =
          std::make_shared<homogeneous_material_factory<detector_t>>();

      DETRAY_DEBUG_HOST("Adding " << mv_data.surface_mat.size()
                                  << " material slabs to material factory");
      for (const auto& mat_data : mv_data.surface_mat) {
        // Determine the material type id
        mat_id mat_type{mat_id::e_none};
        if constexpr (detray::concepts::has_material_rods<detector_t>) {
          if (mat_data.type == io::material_id::rod) {
            mat_type = mat_id::e_material_rod;
          }
        }
        if constexpr (detray::concepts::has_material_slabs<detector_t>) {
          if (mat_data.type == io::material_id::slab) {
            mat_type = mat_id::e_material_slab;
          }
        }
        assert(mat_type != mat_id::e_none);

        const auto mat_idx{mat_data.index_in_coll.has_value()
                               ? mat_data.index_in_coll.value()
                               : detray::detail::invalid_value<std::size_t>()};

        DETRAY_DEBUG_HOST("-> Surface link is: " << mat_idx);

        mat_factory->add_material(mat_type, from_payload<scalar_t>(mat_data),
                                  mat_idx);
      }

      // Add the material to the volume
      DETRAY_DEBUG_HOST("Adding material factory to the builder");
      vm_builder->add_surfaces(mat_factory);
    }
  }

  /// @returns material data for a material factory from a slab io payload
  /// @param slab_data
  template <detray::concepts::scalar scalar_t>
  static material_data<scalar_t> from_payload(
      const surface_material_payload& slab_data) {
    return {static_cast<scalar_t>(slab_data.thickness),
            from_payload<scalar_t>(slab_data.mat),
            detail::basic_converter::from_payload(slab_data.surface)};
  }

  /// @returns the material from its IO payload @param mat_data
  template <detray::concepts::scalar scalar_t>
  static auto from_payload(const material_param_payload& mat_data) {
    return material<scalar_t>{
        static_cast<scalar_t>(mat_data.params[0]),
        static_cast<scalar_t>(mat_data.params[1]),
        static_cast<scalar_t>(mat_data.params[2]),
        static_cast<scalar_t>(mat_data.params[3]),
        static_cast<scalar_t>(mat_data.params[4]),
        // The molar density (params[5]) is calculated on the fly
        static_cast<material_state>(mat_data.params[6])};
  }
};

}  // namespace detray::io
