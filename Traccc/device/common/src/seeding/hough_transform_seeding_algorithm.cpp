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
      m_config{config} {}

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
 
  // Set up the output buffer.
  edm::houghplane_collection<unsigned short, N_COLS, N_ROWS>::const_view::size_type n_slices = m_config.wedges.size();
  edm::houghplane_collection<unsigned short, N_COLS, N_ROWS>::buffer result(
          std::vector<unsigned int>(n_slices, N_COLS * N_ROWS),
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

  // Count and then find the maxima
  vecmem::data::vector_buffer<unsigned int> maxima_sizes(n_slices, mr().main);
  copy().setup(vecmem::ignore_event, maxima_sizes);
  copy().memset(vecmem::ignore_event, maxima_sizes, 0);

  count_maxima_kernel({result, maxima_sizes});
  vecmem::vector<unsigned int> maxima_sizes_host = ((mr().host != nullptr)
                ? vecmem::vector<unsigned int>(mr().host)
                : vecmem::vector<unsigned int>());
  copy()(maxima_sizes, maxima_sizes_host)->wait();
  
  /*
  std::cout << std::endl;
  std::cout << "maxima_sizes_host size: " << maxima_sizes_host.size() << std::endl;
  for (unsigned int i(0); i<maxima_sizes_host.size(); ++i) {
    std::cout << "* n [" << i << "]: " << maxima_sizes_host[i] << std::endl; 
  }
  */

  unsigned int total_maxima = 0u;
  for (unsigned int i(0); i<maxima_sizes_host.size(); ++i) total_maxima += maxima_sizes_host[i];
 
  vecmem::data::jagged_vector_buffer<global_index_t> maxima_global_index(maxima_sizes_host,
        mr().main, mr().host,
        vecmem::data::buffer_type::resizable);
  copy().setup(maxima_global_index)->ignore();

  find_maxima_kernel({result,
                      maxima_global_index});

  /*
  vecmem::data::jagged_vector_view<unsigned int> host_view;
  vecmem::data::jagged_vector_buffer<unsigned int> host_buffer;
  if (mr().host) {
    host_buffer = copy().to(maxima_global_index, *(mr().host));
    host_view = host_buffer;
  } else {
    host_view = maxima_global_index;
  }

  unsigned int totale_trovato = 0u;
  vecmem::jagged_device_vector<unsigned int> host_device{host_view};
  std::cout << "There is a total of " << total_maxima << " maxima" << std::endl;
  std::cout << "Jagged vector (index) has size: " << host_device.size() << std::endl;
  for (const auto& val : host_device) {
    std::cout << "* " << val.size() << std::endl;
    totale_trovato += val.size();
    for (const auto value : val) 
        std::cout << "   * " << value << std::endl;
  }
  std::cout<< ">> total: " << totale_trovato << std::endl;
  */

  // Run clustering on the images

  // Take the pictures from the merged maxima
  edm::houghplane_collection<scalar, 32, 32>::buffer images(
          std::vector<unsigned int>(total_maxima, 32 * 32),
          mr().main, mr().host,
          vecmem::data::buffer_type::resizable);
  vecmem::data::vector_buffer<scalar> maxima_x_value(total_maxima, mr().main);

  copy().setup(vecmem::ignore_event, maxima_x_value);
  copy().memset(vecmem::ignore_event, maxima_x_value, 0);
  copy().setup(images)->ignore();

  take_picture_kernel({total_maxima,
                       result,
                       maxima_global_index,
                       images,
                       maxima_x_value});

 /* 
  edm::houghplane_collection<scalar, 32, 32>::const_view host_view;
  edm::houghplane_collection<scalar, 32, 32>::buffer host_buffer;
  if (mr().host) {
    host_buffer = copy().to(images, *(mr().host));
    host_view = host_buffer;
  } else {
    host_view = images;
  }

  edm::houghplane_collection<scalar, 32, 32>::const_device host_device{host_view};
  std::cout << std::endl;
  std::cout << "----------------- IMAGE 3 -----------------" << std::endl;
  for (unsigned int i(0); i<32; ++i) {
      for (unsigned int j(0); j<32; ++j) {
	    std::cout << host_device[3].plane()[i * 32 + j] << " ";
      }
      std::cout << std::endl;
  }
  std::cout << std::endl;
*/

  // Run NN
  vecmem::data::vector_buffer<unsigned int> nn_outcome(total_maxima, mr().main);
  copy().setup(vecmem::ignore_event, nn_outcome);
  copy().memset(vecmem::ignore_event, nn_outcome, 0);

  // no idea
  synchronize();

  // Return the result.
  return result;
}

}  // namespace traccc::device
