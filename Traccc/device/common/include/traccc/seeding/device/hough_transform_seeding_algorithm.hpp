// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Local include(s).
#include "traccc/device/abstract_awaitable.hpp"
#include "traccc/device/algorithm_base.hpp"

// Project include(s).
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/edm/houghplane_collection.hpp"
#include "traccc/geometry/detector_buffer.hpp"
#include "traccc/utils/algorithm.hpp"
#include "traccc/utils/memory_resource.hpp"
#include "traccc/utils/messaging.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

static constexpr unsigned int N_COLS = 256u;
static constexpr unsigned int N_ROWS = 7000u;

// we need to template the class i guess
class hough_transform_seeding_algorithm
    : public algorithm<edm::houghplane_collection<N_COLS, N_ROWS>::buffer(
          const edm::spacepoint_collection::const_view&)>,
      public messaging,
      public algorithm_base,
      public virtual abstract_awaitable {
 public:
  hough_transform_seeding_algorithm(
      const traccc::memory_resource& mr, const vecmem::copy& copy,
      std::unique_ptr<const Logger> logger = getDummyLogger().clone());

  output_type operator()(const edm::spacepoint_collection::const_view&
                             spacepoints) const override;

 protected:
  struct fill_hough_planes_kernel_payload {
    unsigned int n_cols = N_COLS;
    unsigned int n_rows = N_ROWS;
    /// The number of spacepoints in the event
    edm::spacepoint_collection::const_view::size_type n_spacepoints;
    /// The input spacepoints
    const edm::spacepoint_collection::const_view& spacepoints;
    /// The output hough planes
    edm::houghplane_collection<N_COLS, N_ROWS>::view& houghplanes;
  };

  struct define_plane_edges_kernel_payload {
    unsigned int n_cols = N_COLS;
    unsigned int n_rows = N_ROWS;
    edm::houghplane_collection<N_COLS, N_ROWS>::view& houghplanes;
  };

  virtual void define_plane_edges_kernel(const define_plane_edges_kernel_payload& payload) const = 0;
  virtual void fill_hough_planes_kernel(const fill_hough_planes_kernel_payload& payload) const = 0;

  virtual void synchronize() const = 0;
  };

}  // namespace traccc::device

#include "traccc/seeding/device/impl/define_plane_edges.ipp"

