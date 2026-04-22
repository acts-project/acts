// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/builders/surface_factory_interface.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/definitions/units.hpp"

// Example include(s)
#include "detray/tutorial/detector_metadata.hpp"
#include "detray/tutorial/my_square2D.hpp"

// System include(s)
#include <vector>

namespace detray::tutorial {

/// @brief Generates a sequence of square surfaces for the tutorial detector
class square_surface_generator final
    : public surface_factory_interface<detector<tutorial::my_metadata>> {
 public:
  using detector_t = detector<tutorial::my_metadata>;
  using scalar_t = dscalar<typename detector_t::algebra_type>;

  /// Generate @param n square surfaces with half length @param hl .
  DETRAY_HOST
  square_surface_generator(std::size_t n, scalar_t hl)
      : m_n_squares{static_cast<dindex>(n)}, m_half_length{hl} {}

  /// @returns the number of surfaces this factory will produce
  DETRAY_HOST
  auto size() const -> dindex override { return m_n_squares; }

  /// Generator, does not aggregate any data
  /// @{
  DETRAY_HOST
  void clear() override { /*Do nothing*/ };

  DETRAY_HOST
  void push_back(surface_data<detector_t> &&) override { /*Do nothing*/ }
  DETRAY_HOST
  auto push_back(std::vector<surface_data<detector_t>> &&)
      -> void override { /*Do nothing*/ }
  /// @}

  /// Generate the surfaces and add them to given data collections.
  ///
  /// @param volume the volume that will be added to the detector.
  /// @param surfaces the resulting surface descriptors.
  /// @param transforms the transforms of the surfaces.
  /// @param masks the masks of the surfaces (all of the same shape).
  /// @param ctx the geometry context.
  DETRAY_HOST
  auto operator()(typename detector_t::volume_type &volume,
                  typename detector_t::surface_lookup_container &surfaces,
                  typename detector_t::transform_container &transforms,
                  typename detector_t::mask_container &masks,
                  typename detector_t::geometry_context ctx = {})
      -> dindex_range override {
    using surface_t = typename detector_t::surface_type;
    using mask_link_t = typename surface_t::mask_link;
    using material_link_t = typename surface_t::material_link;

    // Position in the detector mask tuple for square surface masks
    constexpr auto mask_id{detector_t::masks::id::e_square2D};
    // The material will be added in a later step
    constexpr auto no_material = surface_t::material_id::e_none;
    // In case the surfaces container is prefilled with other surfaces
    auto surfaces_offset = static_cast<dindex>(surfaces.size());
    // No ACTS source surface
    constexpr auto invalid_src_link{detail::invalid_value<std::uint64_t>()};

    // Produce a series of square surfaces,
    scalar_t z_translation{0.f};
    for (unsigned int i = 0u; i < m_n_squares; ++i) {
      // Surface placement: no rotation, just translation
      typename detector_t::point3_type translation{0.f, 0.f, z_translation};
      z_translation += 10.f * unit<scalar_t>::mm;
      typename detector_t::transform3_type trf{translation};
      transforms.push_back(trf, ctx);

      // Construct the mask
      masks.template emplace_back<mask_id>(empty_context{}, volume.index(),
                                           m_half_length);

      // Add surface descriptor with all links set (relative to the given
      // containers)
      mask_link_t mask_link{mask_id, masks.template size<mask_id>() - 1u};
      material_link_t material_link{no_material, dindex_invalid};

      surfaces.push_back({transforms.size(ctx) - 1u, mask_link, material_link,
                          volume.index(), surface_id::e_sensitive},
                         invalid_src_link);
    }

    // The range of surfaces that was added to the volume builder
    return {surfaces_offset, static_cast<dindex>(surfaces.size())};
  }

 private:
  /// How many surfaces should be produced
  dindex m_n_squares;
  /// Half length of square
  scalar_t m_half_length;
};

}  // namespace detray::tutorial
