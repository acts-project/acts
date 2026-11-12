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
  const unsigned int n_slices = m_config.n_wedges;
  const unsigned int n_cols = m_config.n_cols;
  const unsigned int n_rows = m_config.n_rows;

  vecmem::data::vector_buffer<scalar> wedges(n_slices, mr().main);
  copy()(vecmem::get_data(m_config.wedges), wedges, vecmem::copy::type::host_to_device)->wait();

  vecmem::data::vector_buffer<unsigned short> thresholds(n_slices, mr().main);
  copy()(vecmem::get_data(m_config.min_contributions), thresholds, vecmem::copy::type::host_to_device)->wait();
 
  // ---------------------------------------------------
  // Fill plane edges
  vecmem::data::vector_buffer<scalar> edges_x(n_cols + 1, mr().main);
  vecmem::data::vector_buffer<scalar> edges_y(n_rows + 1, mr().main);

  copy().setup(vecmem::ignore_event, edges_x);
  copy().memset(vecmem::ignore_event, edges_x, 0);
  copy().setup(vecmem::ignore_event, edges_y);
  copy().memset(vecmem::ignore_event, edges_y, 0);

  define_plane_edges_kernel({m_config,
                             edges_x, edges_y});

  
  // ---------------------------------------------------
  // Pre-compute some quantities
  // ---------------------------------------------------
  // These quantities can be computed once per spacepoint
  // they do not depend on the wedge definition
  vecmem::data::vector_buffer<unsigned long long> spacepoint_wedge_association(n_spacepoints, mr().main);
  vecmem::data::vector_buffer<unsigned long long> spacepoint_detector_id(n_spacepoints, mr().main);
  vecmem::data::vector_buffer<scalar> spacepoint_r(n_spacepoints, mr().main);
  vecmem::data::vector_buffer<scalar> spacepoint_phi(n_spacepoints, mr().main);
  vecmem::data::vector_buffer<scalar> spacepoint_bfield_z(n_spacepoints, mr().main);
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

  precompute_spacepoint_variables_kernel({m_config, bfield,
          wedges,
		  measurements,
		  n_spacepoints, spacepoints,
		  spacepoint_wedge_association,
		  spacepoint_detector_id,
		  spacepoint_r,
		  spacepoint_phi,
          spacepoint_bfield_z});
  


  // -------------------------------------------------------------------------------------  
  /*
  edm::spacepoint_collection::const_view sp_host_view;
  edm::spacepoint_collection::buffer sp_host_buffer;
  if (mr().host) {
    sp_host_buffer = copy().to(spacepoints, *(mr().host));
    sp_host_view = sp_host_buffer;
  } else {
    sp_host_view = spacepoints;
  }
  const edm::spacepoint_collection::const_device sp_host_device{sp_host_view};
  
  vecmem::vector<unsigned long long> spacepoint_detector_id_host = ((mr().host != nullptr)
                ? vecmem::vector<unsigned long long>(mr().host)
                : vecmem::vector<unsigned long long>());
  copy()(spacepoint_detector_id, spacepoint_detector_id_host)->wait(); 

  vecmem::vector<unsigned long long> spacepoint_wedge_association_host = ((mr().host != nullptr)
                ? vecmem::vector<unsigned long long>(mr().host)
                : vecmem::vector<unsigned long long>());
  copy()(spacepoint_wedge_association, spacepoint_wedge_association_host)->wait();
  
  vecmem::vector<scalar> spacepoint_r_host = ((mr().host != nullptr)
                ? vecmem::vector<scalar>(mr().host)
                : vecmem::vector<scalar>());
  copy()(spacepoint_r, spacepoint_r_host)->wait();

  vecmem::vector<scalar> spacepoint_phi_host = ((mr().host != nullptr)
                ? vecmem::vector<scalar>(mr().host)
                : vecmem::vector<scalar>());
  copy()(spacepoint_phi, spacepoint_phi_host)->wait();

  assert(spacepoint_detector_id_host.size() == spacepoint_wedge_association_host.size());
  assert(n_spacepoints == spacepoint_detector_id_host.size());
  for (unsigned int i(0); i<n_spacepoints; ++i) {
    const auto& sp = sp_host_device.at(i);
    unsigned long long detector_id = spacepoint_detector_id_host.at(i);
    unsigned long long wedge_association = spacepoint_wedge_association_host.at(i);
    unsigned int n_bits = count_bits(wedge_association);
    unsigned int layer_number = 0;
    for (unsigned int j(0); j<64; ++j) {
        if ((detector_id & (1ULL << j)) == 0ULL) continue;
        layer_number = j;
        break;
    }
    std::cout << "sp " << i << " x=" << sp.x() << " y=" << sp.y() << " z=" << sp.z()  << " r=" << spacepoint_r_host.at(i) << " phi=" << spacepoint_phi_host.at(i)
              << " bit count: " << n_bits << " [value was " << std::bitset<64>(wedge_association) << "] (id: " << std::bitset<64>(detector_id) << ") layer_number is " << layer_number << std::endl;
  }
  */
  // -------------------------------------------------------------------------------------


  // ---------------------------------------------------
  // Fill the planes with space points
  // ---------------------------------------------------
  vecmem::data::vector_buffer<unsigned long long> houghplane_bins(n_slices * n_cols * n_rows, mr().main);
  copy().setup(vecmem::ignore_event, houghplane_bins);
  copy().memset(vecmem::ignore_event, houghplane_bins, 0);


  fill_hough_plane_kernel({m_config,
                  n_spacepoints, 
                  spacepoint_wedge_association,
                  spacepoint_detector_id,
		          spacepoint_r, 
                  spacepoint_phi,
                  spacepoint_bfield_z,
                  edges_x,
                  edges_y,
                  houghplane_bins});


  // ---------------------------------------------------
  // COUNT BITS
  // ---------------------------------------------------
  
  vecmem::data::vector_buffer<unsigned short> houghplane_entries(n_slices * n_cols * n_rows, mr().main);
  copy().setup(vecmem::ignore_event, houghplane_entries);
  copy().memset(vecmem::ignore_event, houghplane_entries, 0);

  count_bits_kernel({n_slices * n_cols * n_rows,
              houghplane_bins,
              houghplane_entries}); 

  // ---------------------------------------------------
  // Count and then find the maxima
  // ---------------------------------------------------
  vecmem::data::vector_buffer<unsigned int> maxima_sizes(n_slices, mr().main);
  copy().setup(vecmem::ignore_event, maxima_sizes);
  copy().memset(vecmem::ignore_event, maxima_sizes, 0);

  count_maxima_kernel({m_config, 
                       houghplane_entries, 
                       thresholds,
                       maxima_sizes});

  vecmem::vector<unsigned int> maxima_sizes_host = ((mr().host != nullptr)
                ? vecmem::vector<unsigned int>(mr().host)
                : vecmem::vector<unsigned int>());
  copy()(maxima_sizes, maxima_sizes_host)->wait();
 
 
  unsigned int total_maxima = 0u;
  for (unsigned int i(0); i<maxima_sizes_host.size(); ++i) total_maxima += maxima_sizes_host[i];

  vecmem::data::jagged_vector_buffer<global_index_t> maxima_global_index(maxima_sizes_host,
        mr().main, mr().host,
        vecmem::data::buffer_type::resizable);
  copy().setup(maxima_global_index)->ignore();

  find_maxima_kernel({m_config,
                      houghplane_entries,
                      thresholds,
                      maxima_global_index});
 

  /* 
  // -------------------------------------------------------------------------------------
  // planes
  vecmem::vector<unsigned short> houghplane_entries_host = ((mr().host != nullptr)
                ? vecmem::vector<unsigned short>(mr().host)
                : vecmem::vector<unsigned short>());
  copy()(houghplane_entries, houghplane_entries_host)->wait(); 

  // maxima
  vecmem::data::jagged_vector_view<unsigned int> host_view;
  vecmem::data::jagged_vector_buffer<unsigned int> host_buffer;
  if (mr().host) {
    host_buffer = copy().to(maxima_global_index, *(mr().host));
    host_view = host_buffer;
  } else {
    host_view = maxima_global_index;
  }
  vecmem::jagged_device_vector<unsigned int> host_device{host_view};  

  std::cout << std::endl;
  for (unsigned int i(0); i<n_slices; ++i) {
    const unsigned int plane_index = i * n_cols * n_rows;
    const auto& maxima = host_device[i];

    std::cout << "Slice " << i << ": maxima = ";
    for (const auto maximum : maxima) std::cout << maximum << "  ";
    std::cout << "values = ";
    for (unsigned int idx(0); idx<n_cols * n_rows; ++idx) std::cout << houghplane_entries_host[plane_index + idx] << " ";
    std::cout << std::endl;
  }
  throw;
  */


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
  // -------------------------------------------------------------------------------------
  */


  // ---------------------------------------------------
  // Take the pictures from the merged maxima
  // ---------------------------------------------------
  vecmem::data::jagged_vector_buffer<scalar> images(
          std::vector<unsigned int>(total_maxima, 32 * 32),
          mr().main, mr().host,
          vecmem::data::buffer_type::resizable);
  vecmem::data::vector_buffer<scalar> maxima_x_value(total_maxima, mr().main);

  copy().setup(vecmem::ignore_event, maxima_x_value);
  copy().memset(vecmem::ignore_event, maxima_x_value, 0);
  copy().setup(images)->ignore();

  take_picture_kernel({m_config,
                       total_maxima,
                       houghplane_entries,
                       edges_x,
                       maxima_global_index,
                       images,
                       maxima_x_value});


     
  // -------------------------------------------------------------------------------------
  vecmem::data::jagged_vector_view<scalar> host_view;
  vecmem::data::jagged_vector_buffer<scalar> host_buffer;
  if (mr().host) {
    host_buffer = copy().to(images, *(mr().host));
    host_view = host_buffer;
  } else {
    host_view = images;
  }

  vecmem::jagged_device_vector<scalar> host_device{host_view};
  std::cout << std::endl;
  for (unsigned int i(0); i<32; ++i) {
    const auto& image = host_device[i];
    std::cout << "* IMAGE " << i << ": ";
    for (unsigned int idx(0); idx<32*32; ++idx) std::cout << image[idx] << " ";    
    std::cout << std::endl;
  }
  std::cout << std::endl;
  // -------------------------------------------------------------------------------------
  throw;
  


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
 
  count_contributing_hits_kernel({m_config,
                           total_maxima,
                           bfield,
                           edges_x,
                           edges_y,
                           spacepoint_r,
                           spacepoint_phi,
                           spacepoint_bfield_z,
                           spacepoint_wedge_association,
                           maxima_global_index,
                           nn_outcome,
                           hits_sizes});

  vecmem::vector<unsigned int> hits_sizes_host = ((mr().host != nullptr)
                ? vecmem::vector<unsigned int>(mr().host)
                : vecmem::vector<unsigned int>());
  copy()(hits_sizes, hits_sizes_host)->wait();

  /*
  // -------------------------------------------------------------------------------------
  std::cout << std::endl;
  std::cout << "Number of images: " << total_maxima << std::endl;
  std::cout << "hits_sizes_host size: " << hits_sizes_host.size() << std::endl;
  for (unsigned int i(0); i<hits_sizes_host.size(); ++i) {
    std::cout << "* n [" << i << "]: " << hits_sizes_host[i] << std::endl; 
  }

  // unsigned int total_hits = 0u;
  // for (unsigned int i(0); i<hits_sizes_host.size(); ++i) total_hits += hits_sizes_host[i];
  // -------------------------------------------------------------------------------------
  */

  vecmem::data::jagged_vector_buffer<global_index_t> hits(hits_sizes_host,
        mr().main, mr().host,
        vecmem::data::buffer_type::resizable);
  copy().setup(hits)->ignore();


  find_contributing_hits_kernel({m_config,
                            total_maxima,
                            bfield,
                            edges_x,
                            edges_y,
                            spacepoint_r,
                            spacepoint_phi,
                            spacepoint_bfield_z,
                            spacepoint_wedge_association,
                            maxima_global_index,
                            nn_outcome,
                            hits});


  // -------------------------------------------------------------------------------------
  /*
  synchronize();

  vecmem::data::jagged_vector_view<global_index_t> hits_host_view;
  vecmem::data::jagged_vector_buffer<global_index_t> hits_host_buffer;
  if (mr().host) {
    hits_host_buffer = copy().to(hits, *(mr().host));
    hits_host_view = hits_host_buffer;
  } else {
    hits_host_view = hits;
  }

  vecmem::jagged_device_vector<global_index_t> hits_host_device{hits_host_view};
  std::cout << "Number of collections in jagged vector: " << hits_host_device.size() << std::endl;
  std::cout << "Associated hits (from " << total_maxima << " images) [size is " << hits_host_device.size() << "]:" << std::endl;

  for (unsigned int i(0); i<hits_host_device.size(); ++i) {
    unsigned int n_hits = hits_sizes_host[i];
    const auto& vals = hits_host_device[i];
    std::cout << "* image " << i << " should have " << n_hits << " entrie: ";
    for (unsigned int j(0); j<vals.capacity(); ++j) std::cout << vals[j] << " "; 
    std::cout << std::endl;
  }  
  */

  // Make the seeds
  std::vector<unsigned int> matrices_sizes(hits_sizes_host.size());
  std::transform(hits_sizes_host.begin(), hits_sizes_host.end(), matrices_sizes.begin(),
               [](unsigned int x) { return x * x; });

  vecmem::data::jagged_vector_buffer<scalar> matrices(matrices_sizes,
        mr().main, mr().host,
        vecmem::data::buffer_type::resizable);
  copy().setup(matrices)->ignore();

  make_seeds_kernel({total_maxima,
                     spacepoints,
                     spacepoint_r,
                     nn_outcome,
                     hits,
                     matrices});


  // -------------------------------------------------------------------------------------
  /*
  synchronize();

  
  vecmem::data::jagged_vector_view<scalar> matrices_host_view;
  vecmem::data::jagged_vector_buffer<scalar> matrices_host_buffer;
    
  if (mr().host) {
    matrices_host_buffer = copy().to(matrices, *(mr().host));
    matrices_host_view = matrices_host_buffer;
  } else {
    matrices_host_view = matrices;
  }

  vecmem::jagged_device_vector<scalar> matrices_host_device{matrices_host_view};

  for (unsigned int i(0); i<matrices_host_device.size(); ++i) {
    std::cout << "matrix " << i << ": ";
    for (unsigned int j(0); j<matrices_sizes[i]; ++j) std::cout << matrices_host_device[i][j] << " ";
    std::cout << std::endl;
  };
  throw;
  */
  // -------------------------------------------------------------------------------------

  // no idea
  synchronize();


  // Return the result.
  // return result;
  return {};
}

}  // namespace traccc::device
