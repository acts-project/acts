// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Local include(s).
#include "traccc/seeding/device/hough_transform_seeding_algorithm.hpp"

#include <bitset>

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
    const edm::measurement_collection::const_view& measurements,
    const edm::spacepoint_collection::const_view& spacepoints) const
    -> output_type {
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
  const unsigned int n_slices = m_config.n_wedges;
  const unsigned int n_cols = m_config.n_cols;
  const unsigned int n_rows = m_config.n_rows;

  vecmem::data::vector_buffer<scalar> wedges(n_slices, mr().main);
  copy()(vecmem::get_data(m_config.wedges), wedges,
         vecmem::copy::type::host_to_device)
      ->wait();

  vecmem::data::vector_buffer<unsigned short> thresholds(n_slices, mr().main);
  copy()(vecmem::get_data(m_config.min_contributions), thresholds,
         vecmem::copy::type::host_to_device)
      ->wait();

  vecmem::data::vector_buffer<unsigned int> shift_map(m_config.shift_map.size(),
                                                      mr().main);
  copy()(vecmem::get_data(m_config.shift_map), shift_map,
         vecmem::copy::type::host_to_device)
      ->wait();

  // ---------------------------------------------------
  // Fill plane edges
  vecmem::data::vector_buffer<scalar> edges_x(n_cols + 1, mr().main);
  vecmem::data::vector_buffer<scalar> edges_y(n_rows + 1, mr().main);

  copy().setup(vecmem::ignore_event, edges_x);
  copy().memset(vecmem::ignore_event, edges_x, 0);
  copy().setup(vecmem::ignore_event, edges_y);
  copy().memset(vecmem::ignore_event, edges_y, 0);

  define_plane_edges_kernel({m_config, edges_x, edges_y});

  // ---------------------------------------------------
  // Pre-compute some quantities
  // ---------------------------------------------------
  // These quantities can be computed once per spacepoint
  // they do not depend on the wedge definition
  vecmem::data::vector_buffer<unsigned long long> spacepoint_wedge_association(
      n_spacepoints, mr().main);
  vecmem::data::vector_buffer<unsigned long long> spacepoint_detector_id(
      n_spacepoints, mr().main);
  vecmem::data::vector_buffer<scalar> spacepoint_r(n_spacepoints, mr().main);
  vecmem::data::vector_buffer<scalar> spacepoint_phi(n_spacepoints, mr().main);
  vecmem::data::vector_buffer<scalar> spacepoint_bfield_z(n_spacepoints,
                                                          mr().main);
  //  no idea what these do ...
  copy().setup(vecmem::ignore_event, spacepoint_wedge_association);
  copy().memset(vecmem::ignore_event, spacepoint_wedge_association, 0);
  copy().setup(vecmem::ignore_event, spacepoint_detector_id);
  copy().memset(vecmem::ignore_event, spacepoint_detector_id, 0);
  copy().setup(vecmem::ignore_event, spacepoint_r);
  copy().memset(vecmem::ignore_event, spacepoint_r, 0);
  copy().setup(vecmem::ignore_event, spacepoint_phi);
  copy().memset(vecmem::ignore_event, spacepoint_phi, 0);
  copy().setup(vecmem::ignore_event, spacepoint_bfield_z);
  copy().memset(vecmem::ignore_event, spacepoint_bfield_z, 0);

  precompute_spacepoint_variables_kernel(
      {m_config, bfield, wedges, m_config.invalid, shift_map, measurements,
       n_spacepoints, spacepoints, spacepoint_wedge_association,
       spacepoint_detector_id, spacepoint_r, spacepoint_phi,
       spacepoint_bfield_z});

  // ---------------------------------------------------
  // Fill the planes with space points
  // ---------------------------------------------------
  vecmem::data::vector_buffer<unsigned long long> houghplane_bins(
      n_slices * n_cols * n_rows, mr().main);
  copy().setup(vecmem::ignore_event, houghplane_bins);
  copy().memset(vecmem::ignore_event, houghplane_bins, 0);

  fill_hough_plane_kernel({m_config, n_spacepoints,
                           spacepoint_wedge_association, spacepoint_detector_id,
                           spacepoint_r, spacepoint_phi, spacepoint_bfield_z,
                           edges_x, edges_y, houghplane_bins});

  // ---------------------------------------------------
  // COUNT BITS
  // ---------------------------------------------------

  vecmem::data::vector_buffer<unsigned short> houghplane_entries(
      n_slices * n_cols * n_rows, mr().main);
  copy().setup(vecmem::ignore_event, houghplane_entries);
  copy().memset(vecmem::ignore_event, houghplane_entries, 0);

  count_bits_kernel(
      {n_slices * n_cols * n_rows, houghplane_bins, houghplane_entries});

  // ---------------------------------------------------
  // Count and then find the maxima
  // ---------------------------------------------------
  vecmem::data::vector_buffer<unsigned int> maxima_sizes(n_slices, mr().main);
  copy().setup(vecmem::ignore_event, maxima_sizes);
  copy().memset(vecmem::ignore_event, maxima_sizes, 0);

  count_maxima_kernel({m_config, houghplane_entries, thresholds, maxima_sizes});

  vecmem::vector<unsigned int> maxima_sizes_host =
      ((mr().host != nullptr) ? vecmem::vector<unsigned int>(mr().host)
                              : vecmem::vector<unsigned int>());
  copy()(maxima_sizes, maxima_sizes_host)->wait();

  unsigned int total_maxima = 0u;
  for (unsigned int i(0); i < maxima_sizes_host.size(); ++i)
    total_maxima += maxima_sizes_host[i];

  vecmem::data::jagged_vector_buffer<global_index_t> maxima_global_index(
      maxima_sizes_host, mr().main, mr().host,
      vecmem::data::buffer_type::resizable);
  copy().setup(maxima_global_index)->ignore();

  find_maxima_kernel(
      {m_config, houghplane_entries, thresholds, maxima_global_index});

  // ---------------------------------------------------
  // Take the pictures from the merged maxima
  // ---------------------------------------------------
  vecmem::data::jagged_vector_buffer<scalar> images(
      std::vector<unsigned int>(total_maxima, 32 * 32), mr().main, mr().host,
      vecmem::data::buffer_type::resizable);
  vecmem::data::vector_buffer<scalar> maxima_x_value(total_maxima, mr().main);

  copy().setup(vecmem::ignore_event, maxima_x_value);
  copy().memset(vecmem::ignore_event, maxima_x_value, 0);
  copy().setup(images)->ignore();

  take_picture_kernel({m_config, total_maxima, houghplane_entries, edges_x,
                       maxima_global_index, images, maxima_x_value});

  // ---------------------------------------------------
  // Run NN
  // ---------------------------------------------------
  vecmem::data::vector_buffer<unsigned int> nn_outcome(total_maxima, mr().main);
  copy().setup(vecmem::ignore_event, nn_outcome);
  copy().memset(vecmem::ignore_event, nn_outcome, 0);

  evaluate_nn_kernel({total_maxima, images, maxima_x_value, nn_outcome});

  // ---------------------------------------------------
  // Count and find the contributing hits
  // ---------------------------------------------------
  vecmem::data::vector_buffer<unsigned int> hits_sizes(total_maxima, mr().main);
  copy().setup(vecmem::ignore_event, hits_sizes);
  copy().memset(vecmem::ignore_event, hits_sizes, 0);

  count_contributing_hits_kernel(
      {m_config, total_maxima, bfield, edges_x, edges_y, spacepoint_r,
       spacepoint_phi, spacepoint_bfield_z, spacepoint_wedge_association,
       maxima_global_index, nn_outcome, hits_sizes});

  vecmem::vector<unsigned int> hits_sizes_host =
      ((mr().host != nullptr) ? vecmem::vector<unsigned int>(mr().host)
                              : vecmem::vector<unsigned int>());
  copy()(hits_sizes, hits_sizes_host)->wait();

  vecmem::data::jagged_vector_buffer<global_index_t> hits(
      hits_sizes_host, mr().main, mr().host,
      vecmem::data::buffer_type::resizable);
  copy().setup(hits)->ignore();

  find_contributing_hits_kernel(
      {m_config, total_maxima, bfield, edges_x, edges_y, spacepoint_r,
       spacepoint_phi, spacepoint_bfield_z, spacepoint_wedge_association,
       maxima_global_index, nn_outcome, hits});

  // Compute cot-theta
  std::vector<unsigned int> matrices_sizes(hits_sizes_host.size());
  std::transform(hits_sizes_host.begin(), hits_sizes_host.end(),
                 matrices_sizes.begin(), [](unsigned int x) { return x * x; });

  vecmem::data::jagged_vector_buffer<scalar> matrices(
      matrices_sizes, mr().main, mr().host,
      vecmem::data::buffer_type::resizable);
  copy().setup(matrices)->ignore();

  unsigned int max_hits = 0u;
  for (unsigned int value : hits_sizes_host)
    max_hits = std::max(max_hits, value);
  compute_cottheta_matrix_kernel({total_maxima, max_hits, m_config.z_window,
                                  spacepoints, spacepoint_r, nn_outcome, hits,
                                  matrices});

  // Filter hits
  vecmem::data::jagged_vector_buffer<unsigned int> iotas(matrices_sizes,
                                                         mr().main, mr().host);
  copy().setup(iotas)->ignore();

  vecmem::data::vector_buffer<unsigned int> n_filtered_hits(total_maxima,
                                                            mr().main);
  copy().setup(vecmem::ignore_event, n_filtered_hits);
  copy().memset(vecmem::ignore_event, n_filtered_hits, 0);

  vecmem::data::vector_buffer<unsigned int> best_hit_indexes(total_maxima,
                                                             mr().main);
  copy().setup(vecmem::ignore_event, best_hit_indexes);
  copy().memset(vecmem::ignore_event, best_hit_indexes, 0);

  vecmem::data::vector_buffer<scalar> mode_values(total_maxima, mr().main);
  copy().setup(vecmem::ignore_event, mode_values);
  copy().memset(vecmem::ignore_event, mode_values, 0);

  select_hits_kernel({total_maxima, nn_outcome, hits, matrices, iotas,
                      mode_values, n_filtered_hits, best_hit_indexes});

  vecmem::vector<unsigned int> n_filtered_hits_host =
      ((mr().host != nullptr) ? vecmem::vector<unsigned int>(mr().host)
                              : vecmem::vector<unsigned int>());
  copy()(n_filtered_hits, n_filtered_hits_host)->wait();

  vecmem::data::jagged_vector_buffer<global_index_t> selected_hits(
      n_filtered_hits_host, mr().main, mr().host,
      vecmem::data::buffer_type::resizable);
  copy().setup(selected_hits)->ignore();

  // Create the result object.
  edm::seed_collection::buffer seed_buffer(total_maxima, mr().main);
  copy().setup(vecmem::ignore_event, seed_buffer);

  make_seeds_kernel({total_maxima, spacepoint_r, nn_outcome, hits, matrices,
                     mode_values, best_hit_indexes, selected_hits,
                     seed_buffer});

  // synchronization
  synchronize();

  // Return the result.
  return seed_buffer;
}

}  // namespace traccc::device
