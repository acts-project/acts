// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/core/concepts.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/geometry/surface.hpp"
#include "detray/geometry/tracking_volume.hpp"
#include "detray/io/backend/detail/basic_converter.hpp"
#include "detray/io/backend/detail/type_info.hpp"
#include "detray/io/frontend/payloads.hpp"
#include "detray/io/frontend/writer_interface.hpp"
#include "detray/material/concepts.hpp"
#include "detray/material/material_rod.hpp"
#include "detray/material/material_slab.hpp"
#include "detray/utils/type_list.hpp"

// System include(s)
#include <string_view>

namespace detray::io {

/// @brief Homogeneous material writer backend
///
/// Fills a @c detector_homogeneous_material_payload from a @c detector instance
template <detray::concepts::detector detector_t>
class homogeneous_material_writer final : public writer_interface<detector_t> {
 public:
  /// Tag the writer as "homogeneous_material"
  static constexpr std::string_view s_tag = "homogeneous_material";

  /// Payload type that the writer processes
  using payload_type = detector_homogeneous_material_payload;

  /// @returns the tag of the writer: "homogeneous_material"
  std::string_view tag() const override { return s_tag; }

  /// Convert the homogeneous material of the detector @param det into the
  /// homogeneous material payload of @param det_data
  void to_payload(const detector_t& det,
                  const typename detector_t::name_map& names,
                  detector_payload& det_data) const override {
    det_data.homogeneous_material = to_payload_impl(det, names);
  }

  /// Convert the material description of a detector @param det into its io
  /// payload
  static payload_type to_payload_impl(
      const detector_t& det, const typename detector_t::name_map& /*unused*/) {
    payload_type dm_data;
    dm_data.volumes.reserve(det.volumes().size());

    for (const auto& vol : det.volumes()) {
      dm_data.volumes.push_back(to_payload_impl(vol, det));
    }

    return dm_data;
  }

  /// Convert the material description of a volume @param vol_desc into its
  /// io payload
  static material_volume_payload to_payload_impl(
      const typename detector_t::volume_type& vol_desc, const detector_t& det) {
    using material_type = surface_material_payload::mat_type;

    material_volume_payload mv_data;
    mv_data.volume_link = detail::basic_converter::to_payload(vol_desc.index());

    // Return early if the stores for homogeneous materials are empty
    using mat_id = typename detector_t::material::id;
    using algebra_t = typename detector_t::algebra_type;

    bool contains_rods{true};
    bool contains_slabs{true};
    if constexpr (detray::concepts::has_material_rods<detector_t>) {
      if (det.material_store().template empty<mat_id::e_material_rod>()) {
        contains_rods = false;
      }
    }
    if constexpr (detray::concepts::has_material_slabs<detector_t>) {
      if (det.material_store().template empty<mat_id::e_material_slab>()) {
        contains_slabs = false;
      }
    }
    if (!contains_rods && !contains_slabs) {
      return mv_data;
    }

    // Find all surfaces that belong to the volume and count them
    std::size_t sf_idx{0u};
    std::size_t slab_idx{0u};
    std::size_t rod_idx{0u};
    auto vol = tracking_volume{det, vol_desc};
    for (const auto& sf_desc : vol.surfaces()) {
      // Convert material slabs and rods
      const auto sf = geometry::surface{det, sf_desc};

      if (surface_material_payload mslp =
              sf.has_material()
                  ? sf.template visit_material<get_material_payload<algebra_t>>(
                        sf_idx)
                  : surface_material_payload{};
          mslp.type == material_type::slab) {
        mslp.index_in_coll = slab_idx++;
        mv_data.surface_mat.push_back(mslp);
      } else if (mslp.type == material_type::rod) {
        mslp.index_in_coll = rod_idx++;
        mv_data.surface_mat.push_back(mslp);
      } else { /* material maps are handled by another writer */
      }
      ++sf_idx;
    }

    return mv_data;
  }

  /// Convert surface material @param mat into its io payload
  template <detray::concepts::material_params material_t>
  static material_param_payload to_payload_impl(const material_t& mat) {
    material_param_payload mat_data;

    mat_data.params = {mat.X0(),
                       mat.L0(),
                       mat.Ar(),
                       mat.Z(),
                       mat.mass_density(),
                       mat.molar_density(),
                       static_cast<io::scalar>(mat.state())};
    return mat_data;
  }

  /// Convert a surface material slab @param mat_slab into its io payload
  template <detray::concepts::algebra algebra_t,
            detray::concepts::homogeneous_material material_t>
  static surface_material_payload to_payload_impl(const material_t& mat,
                                                  std::size_t sf_idx) {
    surface_material_payload mat_data;

    mat_data.type = io::detail::get_id<algebra_t, material_t>();
    mat_data.surface = detail::basic_converter::to_payload(sf_idx);
    mat_data.thickness = mat.thickness();
    mat_data.mat = to_payload_impl(mat.get_material());

    return mat_data;
  }

 private:
  /// Retrieve @c surface_material_payload from a material store element
  template <detray::concepts::algebra algebra_t>
  struct get_material_payload {
    template <typename material_group_t, typename index_t>
    inline auto operator()(const material_group_t& material_group,
                           const index_t& index,
                           [[maybe_unused]] std::size_t sf_index) const {
      using material_t = typename material_group_t::value_type;

      if constexpr (detray::concepts::material_slab<material_t> ||
                    detray::concepts::material_rod<material_t>) {
        return homogeneous_material_writer::template to_payload_impl<algebra_t>(
            material_group[index], sf_index);
      } else {
        return surface_material_payload{};
      }
    }
  };
};

}  // namespace detray::io
