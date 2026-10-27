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
#include "traccc/seeding/detail/hough_seeding_config.hpp"

// Project include(s).
#include "traccc/bfield/magnetic_field.hpp"
#include "traccc/edm/measurement_collection.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/edm/houghplane_collection.hpp"
#include "traccc/geometry/detector_buffer.hpp"
#include "traccc/utils/algorithm.hpp"
#include "traccc/utils/memory_resource.hpp"
#include "traccc/utils/messaging.hpp"
#include "traccc/device/global_index.hpp"

// VecMem include(s).
#include <vecmem/containers/data/vector_view.hpp>

namespace traccc::device {

static constexpr unsigned int N_COLS = 256u;
static constexpr unsigned int N_ROWS = 7000u;

// we need to template the class i guess
class hough_transform_seeding_algorithm
    : public algorithm<edm::houghplane_collection<unsigned short, N_COLS, N_ROWS>::buffer(
	  const magnetic_field&,
	  const detector_buffer&,
	  const edm::measurement_collection::const_view&,
	  const edm::spacepoint_collection::const_view&)>,
      public messaging,
      public algorithm_base,
      public virtual abstract_awaitable {
 public:
   using config_type = traccc::hough_seeding_config;
	 
   hough_transform_seeding_algorithm(
      const traccc::memory_resource& mr, const vecmem::copy& copy,
      const traccc::hough_seeding_config& config,
      std::unique_ptr<const Logger> logger = getDummyLogger().clone());

  output_type operator()(
		  const magnetic_field&,
		  const detector_buffer& det,
		  const edm::measurement_collection::const_view& measurements,
		  const edm::spacepoint_collection::const_view& spacepoints) const override;

 protected:
  struct define_plane_edges_kernel_payload {
    const hough_seeding_config& config;
    edm::houghplane_collection<unsigned short, N_COLS, N_ROWS>::view& houghplanes;
  };

  struct precompute_spacepoint_variables_payload {
    const hough_seeding_config& config;
    /// The input measurements
    const edm::measurement_collection::const_view& measurements;
    /// The number of spacepoints in the event
    const edm::spacepoint_collection::const_view::size_type n_spacepoints;
    /// The number of spacepoints in the event
    const edm::spacepoint_collection::const_view& spacepoints;
    // association to wedge
    vecmem::data::vector_view<unsigned long long>& spacepoint_wedge_association;
    // the detector id to use in the HT filling
    vecmem::data::vector_view<unsigned long long>& spacepoint_detector_id;
    // the space point r
    vecmem::data::vector_view<scalar>& spacepoint_r;
    // the space point phi
    vecmem::data::vector_view<scalar>& spacepoint_phi;
  };

  struct fill_hough_planes_kernel_payload {
    const hough_seeding_config& config;
    /// The magnetic field object
    const magnetic_field& bfield;
    /// The number of spacepoints in the event
    const edm::spacepoint_collection::const_view::size_type n_spacepoints;
    /// The input spacepoints
    const edm::spacepoint_collection::const_view& spacepoints;
    // association to wedge
    vecmem::data::vector_view<unsigned long long>& spacepoint_wedge_association;
    // the detector id to use in the HT filling
    vecmem::data::vector_view<unsigned long long>& spacepoint_detector_id;
    // the space point r
    vecmem::data::vector_view<scalar>& spacepoint_r;
    // the space point phi
    vecmem::data::vector_view<scalar>& spacepoint_phi;
    /// The output hough planes
    edm::houghplane_collection<unsigned short, N_COLS, N_ROWS>::view& houghplanes;
  };

  struct count_maxima_payload {
    // The hough planes
    const edm::houghplane_collection<unsigned short, N_COLS, N_ROWS>::view& houghplanes;
    // The count vector
    vecmem::data::vector_view<unsigned int>& count_vector;
  };
 
  struct find_maxima_payload {
    // The hough planes
    const edm::houghplane_collection<unsigned short, N_COLS, N_ROWS>::view& houghplanes;
    // The global index vector
    vecmem::data::jagged_vector_view<global_index_t>& global_index;
  };
 
  struct take_picture_payload {
    // Number of maxima
    unsigned int n_maxima;
    // The hough planes
    const edm::houghplane_collection<unsigned short, N_COLS, N_ROWS>::view& houghplanes;
    // The maxima location
    const vecmem::data::jagged_vector_view<global_index_t>& global_index;
    // The images
    edm::houghplane_collection<scalar, 32, 32>::view& images;
  };

  virtual void define_plane_edges_kernel(const define_plane_edges_kernel_payload& payload) const = 0;
  virtual void precompute_spacepoint_variables_kernel(const precompute_spacepoint_variables_payload& payload) const = 0;
  virtual void fill_hough_plane_kernel(const fill_hough_planes_kernel_payload& payload) const = 0;
  virtual void count_maxima_kernel(const count_maxima_payload& payload) const = 0;
  virtual void find_maxima_kernel(const find_maxima_payload& payload) const = 0;
  virtual void take_picture_kernel(const take_picture_payload& payload) const = 0;

  virtual void synchronize() const = 0;
  
  config_type m_config;    
  };

}  // namespace traccc::device

#include "traccc/seeding/device/impl/hough_seeding_utilities.ipp"
#include "traccc/seeding/device/impl/define_plane_edges.ipp"
#include "traccc/seeding/device/impl/fill_hough_plane.ipp"
#include "traccc/seeding/device/impl/precompute_spacepoint_variables.ipp"
#include "traccc/seeding/device/impl/count_maxima.ipp"
#include "traccc/seeding/device/impl/find_maxima.ipp"
#include "traccc/seeding/device/impl/take_pictures.ipp"

