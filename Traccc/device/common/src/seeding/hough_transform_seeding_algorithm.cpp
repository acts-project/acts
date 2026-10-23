// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/seeding/device/hough_transform_seeding_algorithm.hpp"

namespace traccc::device {


hough_transform_seeding_algorithm::hough_transform_seeding_algorithm(
    const traccc::memory_resource& mr, const vecmem::copy& copy,
    const traccc::hough_seeding_config& config,
    std::unique_ptr<const Logger> logger)
    : messaging(std::move(logger)),
      algorithm_base{mr, copy},
      m_config(config) {}

auto hough_transform_seeding_algorithm::operator()(
		const magnetic_field& bfield,
		const detector_buffer& /*det*/,
		const edm::measurement_collection::const_view& measurements,
		const edm::spacepoint_collection::const_view& spacepoints) const -> output_type {
      	// Get the number of space points. In an asynchronous way if possible.
  edm::spacepoint_collection::const_view::size_type n_spacepoints = 0u;
  if (mr().host) {
    vecmem::async_size size = copy().get_size(spacepoints, *(mr().host));
    // Block or suspend execution until the size is available.
    await(size);
    n_spacepoints = size.unsafe_get();
  } else {
    n_spacepoints = copy().get_size(spacepoints);
  }

  // If there are no spacepoints, return right away.
  if (n_spacepoints == 0) {
    return {};
  }
 
  edm::houghplane_collection<unsigned short, N_COLS, N_ROWS>::const_view::size_type n_slices = m_config.wedges.size(); 
  vecmem::vector<unsigned int> plane_sizes_host =
        ((mr().host != nullptr) ? vecmem::vector<unsigned int>(n_slices, mr().host)
                                : vecmem::vector<unsigned int>(n_slices));
  for (unsigned int i = 0; i < n_slices; ++i) plane_sizes_host[i] = N_COLS * N_ROWS;

  // Set up the output buffer.
  edm::houghplane_collection<unsigned short, N_COLS, N_ROWS>::buffer result(plane_sizes_host,
		  mr().main, mr().host, 
		  vecmem::data::buffer_type::resizable);
  copy().setup(result)->ignore();

  // Fill plane edges
  // we will move to a config at a certain point and template the algorithm
  define_plane_edges_kernel({m_config, result});
  
  // Pre-compute some quantities
  // These quantities can be computed once per spacepoint
  // they do not depend on the wedge definition
  vecmem::data::vector_buffer<unsigned long long> spacepoint_wedge_association(n_spacepoints, mr().main);
  vecmem::data::vector_buffer<unsigned long long> spacepoint_detector_id(n_spacepoints, mr().main);
  vecmem::data::vector_buffer<scalar> spacepoint_r(n_spacepoints, mr().main);
  vecmem::data::vector_buffer<scalar> spacepoint_phi(n_spacepoints, mr().main);
  //  no idea what these do ...
  copy().setup(vecmem::ignore_event, spacepoint_wedge_association);
  copy().memset(vecmem::ignore_event, spacepoint_wedge_association, 0);
  copy().setup(vecmem::ignore_event, spacepoint_detector_id);
  copy().memset(vecmem::ignore_event, spacepoint_detector_id, 0);
  copy().setup(vecmem::ignore_event, spacepoint_r);
  copy().memset(vecmem::ignore_event, spacepoint_r, 0);
  copy().setup(vecmem::ignore_event, spacepoint_phi);
  copy().memset(vecmem::ignore_event, spacepoint_phi, 0);

  precompute_spacepoint_variables_kernel({m_config,
		  measurements,
		  n_spacepoints, spacepoints,
		  spacepoint_wedge_association,
		  spacepoint_detector_id,
		  spacepoint_r,
		  spacepoint_phi});
  
  // Fill the planes with space points
  fill_hough_plane_kernel({m_config, bfield,
                  n_spacepoints, spacepoints,
                  spacepoint_wedge_association, spacepoint_detector_id,
		  spacepoint_r, spacepoint_phi,
		  result});

/*  
  synchronize();


  edm::houghplane_collection<unsigned short, N_COLS, N_ROWS>::const_view host_view;
  edm::houghplane_collection<unsigned short, N_COLS, N_ROWS>::buffer host_buffer;
  if (mr().host) {
    host_buffer = copy().to(result, *(mr().host));
    host_view = host_buffer;
  } else {
    host_view = result;
  }

  edm::houghplane_collection<unsigned short, N_COLS, N_ROWS>::const_device host_device{host_view};
//  std::cout << std::endl;
  std::cout << "----------------- SLICE 13 -----------------" << std::endl;
  for (unsigned int i(0); i<N_COLS; ++i) {
      for (unsigned int j(0); j<N_ROWS; ++j) {
	// if (host_device[13].plane()[i * N_ROWS + j] != 0) std::cout << "SOMETHING !!!" << std::endl;
	std::cout << host_device[13].plane()[i * N_ROWS + j] << " ";
      }
      std::cout << std::endl;
  }
  std::cout << std::endl;

  /*
  edm::houghplane_collection<unsigned short, N_COLS, N_ROWS>::const_device host_device{host_view};
  for (unsigned int i = 0; i < host_device.size(); ++i) {
    const edm::houghplane plane = host_device.at(i);
    const auto& edges_x = plane.edges_x();
    const auto& edges_y = plane.edges_y();
    std::cout << "----------------" << std::endl; 
    std::cout << "n edges x: " << edges_x.size() << ": " << edges_x[0] << " -> " << edges_x[N_COLS] << std::endl; 
    std::cout << "n edges y: " << edges_y.size() << ": " << edges_y[0] << " -> " << edges_y[N_ROWS] << std::endl;
  }
*/

  synchronize();

  // Return the result.
  return result;
}

}  // namespace traccc::device
