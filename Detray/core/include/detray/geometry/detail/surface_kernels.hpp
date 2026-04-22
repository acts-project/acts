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
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/material/concepts.hpp"
#include "detray/material/detail/material_accessor.hpp"
#include "detray/material/material.hpp"
#include "detray/utils/concepts.hpp"
#include "detray/utils/ranges.hpp"

// System include(s)
#include <limits>
#include <ostream>
#include <ranges>

namespace detray::detail {

/// Functors to be used in the @c surface class
template <concepts::algebra algebra_t>
struct surface_kernels {
  using algebra_type = algebra_t;
  using scalar_type = dscalar<algebra_t>;
  using point2_type = dpoint2D<algebra_t>;
  using point3_type = dpoint3D<algebra_t>;
  using vector3_type = dvector3D<algebra_t>;
  using transform3_type = dtransform3D<algebra_t>;

  /// A functor to retrieve the masks shape name
  struct get_shape_name {
    template <typename mask_group_t, typename index_t>
    DETRAY_HOST inline std::string operator()(const mask_group_t&,
                                              const index_t) const {
      return std::string(mask_group_t::value_type::shape::name);
    }
  };

  /// A functor that checks if a local point @param loc_p is within the
  /// surface mask with tolerance @param tol
  struct is_inside {
    template <typename mask_group_t, typename idx_range_t>
    DETRAY_HOST_DEVICE constexpr bool operator()(const mask_group_t& mask_group,
                                                 const idx_range_t idx_range,
                                                 const point3_type& loc_p,
                                                 const scalar_type tol) const {
      for (const auto& mask : detray::ranges::subrange(mask_group, idx_range)) {
        if (mask.is_inside(loc_p, tol)) {
          return true;
        }
      }

      return false;
    }
  };

  /// A functor to run the mask self check. Puts error messages into @param os
  struct mask_self_check {
    template <typename mask_group_t, typename idx_range_t>
    DETRAY_HOST_DEVICE constexpr auto operator()(const mask_group_t& mask_group,
                                                 const idx_range_t idx_range,
                                                 std::ostream& os) const {
      bool is_good{true};
      for (const auto& mask : detray::ranges::subrange(mask_group, idx_range)) {
        is_good &= mask.self_check(os);
      }

      return is_good;
    }
  };

  /// A functor to retrieve the masks volume link
  struct get_volume_links {
    template <typename mask_group_t, typename idx_range_t>
    DETRAY_HOST inline auto operator()(const mask_group_t& mask_group,
                                       const idx_range_t idx_range) const {
      using volume_link_t = typename mask_group_t::value_type::links_type;

      dvector<volume_link_t> volume_links{};

      for (const auto& mask : detray::ranges::subrange(mask_group, idx_range)) {
        volume_links.push_back(mask.volume_link());
      }

      return volume_links;
    }
  };

  /// A functor to retrieve a mask boundary, determined by @param i
  struct get_mask_value {
    template <typename mask_group_t, typename idx_range_t>
    DETRAY_HOST inline auto operator()(const mask_group_t& mask_group,
                                       const idx_range_t idx_range,
                                       const std::size_t i) const {
      dvector<scalar_type> values{};

      for (const auto& mask : detray::ranges::subrange(mask_group, idx_range)) {
        values.push_back(mask[i]);
      }

      return values;
    }
  };

  /// A functor to retrieve the mask boundaries (host only)
  struct get_mask_values {
    template <typename mask_group_t, typename idx_range_t>
    DETRAY_HOST inline auto operator()(const mask_group_t& mask_group,
                                       const idx_range_t idx_range) const {
      dvector<dvector<scalar_type>> values{};

      for (const auto& mask : detray::ranges::subrange(mask_group, idx_range)) {
        dvector<scalar_type>& sub_mask_values = values.emplace_back();
        std::ranges::copy(mask.values(), std::back_inserter(sub_mask_values));
      }

      return values;
    }
  };

  /// A functor to retrieve the material parameters
  struct get_material_params {
    template <typename mat_group_t, concepts::index index_t>
    DETRAY_HOST_DEVICE constexpr auto operator()(
        const mat_group_t& mat_group, const index_t idx,
        const point2_type& loc_p) const {
      using material_t = typename mat_group_t::value_type;

      if constexpr (concepts::surface_material<material_t>) {
        return &(detail::material_accessor::get(mat_group, idx, loc_p)
                     .get_material());
      } else {
        using scalar_t = typename material_t::scalar_type;
        // Volume material (cannot be reached from a surface)
        return static_cast<const material<scalar_t>*>(nullptr);
      }
    }
  };

  /// A functor get the surface normal at a given local/bound position
  struct normal {
    template <typename mask_group_t, typename index_t>
    DETRAY_HOST_DEVICE constexpr point3_type operator()(
        const mask_group_t& /*mask_group*/, const index_t /*index*/,
        const transform3_type& trf3, const point3_type& local) const {
      using mask_t = typename mask_group_t::value_type;

      return mask_t::get_local_frame().normal(trf3, local);
    }

    template <typename mask_group_t, concepts::index index_t>
    DETRAY_HOST_DEVICE constexpr point3_type operator()(
        const mask_group_t& mask_group, const index_t index,
        const transform3_type& trf3, const point2_type& bound) const {
      using mask_t = typename mask_group_t::value_type;

      return mask_t::get_local_frame().normal(trf3, bound, mask_group[index]);
    }

    template <typename mask_group_t, concepts::interval idx_range_t>
    DETRAY_HOST_DEVICE constexpr point3_type operator()(
        const mask_group_t& mask_group, const idx_range_t idx_range,
        const transform3_type& trf3, const point2_type& bound) const {
      using mask_t = typename mask_group_t::value_type;

      // The normal is defined by the coordinate system of the surface and
      // the local position: take any mask
      const mask_t& mask = mask_group[idx_range.lower()];
      return mask_t::get_local_frame().normal(trf3, bound, mask);
    }
  };

  /// A functor get the mask centroid in local cartesian coordinates
  struct centroid {
    template <typename mask_group_t, concepts::index index_t>
    DETRAY_HOST_DEVICE constexpr point3_type operator()(
        const mask_group_t& mask_group, const index_t index) const {
      return mask_group[index].centroid();
    }

    template <typename mask_group_t, concepts::interval idx_range_t>
    DETRAY_HOST_DEVICE constexpr point3_type operator()(
        const mask_group_t& mask_group, const idx_range_t idx_range) const {
      using mask_t = typename mask_group_t::value_type;
      using index_t = typename idx_range_t::index_type;

      if (idx_range.size() > 1u) {
        // Find the true surface extent over all masks
        mask_t sf_mask = mask_group[idx_range.lower()];

        const idx_range_t other_masks{
            idx_range.lower() + 1u,
            static_cast<index_t>(idx_range.size() - 1u)};

        // Merge sub-masks
        for (const auto& sub_mask :
             detray::ranges::subrange(mask_group, other_masks)) {
          sf_mask = sf_mask + sub_mask;
        }

        return sf_mask.centroid();
      } else {
        return mask_group[idx_range.lower()].centroid();
      }
    }
  };

  /// A functor to perform global to local transformation
  struct global_to_local {
    template <typename mask_group_t, typename index_t>
    DETRAY_HOST_DEVICE constexpr point3_type operator()(
        const mask_group_t& /*mask_group*/, const index_t& /*index*/,
        const transform3_type& trf3, const point3_type& global,
        const vector3_type& dir) const {
      using mask_t = typename mask_group_t::value_type;

      return mask_t::to_local_frame3D(trf3, global, dir);
    }
  };

  /// A functor to perform local to global transformation
  struct local_to_global {
    template <typename mask_group_t, concepts::index index_t>
    DETRAY_HOST_DEVICE constexpr point3_type operator()(
        const mask_group_t& mask_group, const index_t index,
        const transform3_type& trf3, const point2_type& bound,
        const vector3_type& dir) const {
      using mask_t = typename mask_group_t::value_type;

      return mask_t::get_local_frame().local_to_global(trf3, mask_group[index],
                                                       bound, dir);
    }

    template <typename mask_group_t, concepts::interval idx_range_t>
    DETRAY_HOST_DEVICE constexpr point3_type operator()(
        const mask_group_t& mask_group, const idx_range_t idx_range,
        const transform3_type& trf3, const point2_type& bound,
        const vector3_type& dir) const {
      using mask_t = typename mask_group_t::value_type;

      // The global position is defined by the coordinate system of the
      // surface and the local position: take any mask
      const mask_t& mask = mask_group[idx_range.lower()];
      return mask_t::get_local_frame().local_to_global(trf3, mask, bound, dir);
    }

    template <typename mask_group_t, typename index_t>
    DETRAY_HOST_DEVICE constexpr point3_type operator()(
        const mask_group_t& /*mask_group*/, const index_t& /*index*/,
        const transform3_type& trf3, const point3_type& local,
        const vector3_type&) const {
      using mask_t = typename mask_group_t::value_type;

      return mask_t::to_global_frame(trf3, local);
    }
  };
};

}  // namespace detray::detail
