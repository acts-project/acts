// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/core/detail/container_buffers.hpp"
#include "detray/core/detail/container_views.hpp"
#include "detray/definitions/grid_axis.hpp"
#include "detray/definitions/indexing.hpp"
#include "detray/geometry/surface_descriptor.hpp"
#include "detray/utils/concepts.hpp"
#include "detray/utils/ranges/ranges.hpp"

// System include(s)
#include <concepts>
#include <utility>

namespace detray::concepts {

template <typename A>
concept axis = requires(const A ax) {
  typename A::bounds_type;
  typename A::binning_type;
  typename A::scalar_type;

  { ax.label() } -> std::same_as<axis::label>;

  { ax.bounds() } -> std::same_as<axis::bounds>;

  { ax.binning() } -> std::same_as<axis::binning>;

  { ax.nbins() } -> std::same_as<dindex>;

  { ax.bin(typename A::scalar_type()) } -> std::same_as<dindex>;

  {
    ax.range(typename A::scalar_type(), darray<dindex, 2>())
  } -> std::same_as<darray<int, 2>>;

  {
    ax.range(typename A::scalar_type(), darray<typename A::scalar_type, 2>())
  } -> std::same_as<darray<int, 2>>;

  {
    ax.bin_edges(dindex())
  } -> std::same_as<darray<typename A::scalar_type, 2>>;

  { ax.bin_edges() } -> detray::ranges::range_of<typename A::scalar_type>;

  { ax.span() } -> std::same_as<darray<typename A::scalar_type, 2>>;

  { ax.min() } -> std::same_as<typename A::scalar_type>;

  { ax.max() } -> std::same_as<typename A::scalar_type>;
};

template <typename G>
concept grid = viewable<G> && bufferable<G> && requires(const G g) {
  typename G::bin_type;
  typename G::value_type;
  typename G::glob_bin_index;
  typename G::loc_bin_index;
  typename G::local_frame_type;
  typename G::algebra_type;
  typename G::point_type;

  G::dim;
  G::is_owning;

  // TODO: Implement coordinate frame concept
  { g.get_local_frame() } -> std::same_as<typename G::local_frame_type>;

  { g.template get_axis<0>() } -> concepts::axis;

  { g.nbins() } -> std::same_as<dindex>;

  { g.size() } -> std::same_as<dindex>;

  {
    g.deserialize(typename G::glob_bin_index())
  } -> std::same_as<typename G::loc_bin_index>;

  {
    g.serialize(typename G::loc_bin_index())
  } -> std::same_as<typename G::glob_bin_index>;

  { g.bins() } -> detray::ranges::range_of<typename G::bin_type>;

  {
    g.bin(typename G::glob_bin_index())
  } -> concepts::same_as_cvref<typename G::bin_type>;

  {
    g.bin(typename G::loc_bin_index())
  } -> concepts::same_as_cvref<typename G::bin_type>;

  {
    g.at(typename G::loc_bin_index(), dindex())
  } -> concepts::same_as_cvref<typename G::value_type>;

  {
    g.at(typename G::glob_bin_index(), dindex())
  } -> concepts::same_as_cvref<typename G::value_type>;
};

/// Grid that contains surfaces (used e.g. for grid acceleration structures)
/// TODO: Add surface descriptor concept to geometry package
template <class G>
concept surface_grid =
    concepts::grid<G> &&
    std::same_as<typename G::value_type,
                 surface_descriptor<typename G::value_type::mask_link,
                                    typename G::value_type::material_link,
                                    typename G::value_type::transform_link,
                                    typename G::value_type::navigation_link>>;

/// Acceleration structure that contains volumes (volume indices)
template <class G>
concept volume_grid =
    concepts::grid<G> && std::same_as<typename G::value_type, dindex>;

}  // namespace detray::concepts
