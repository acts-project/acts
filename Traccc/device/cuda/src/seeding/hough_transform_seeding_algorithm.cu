// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

// Library include(s).
#include "../utils/cuda_error_handling.hpp"
#include "../utils/global_index.hpp"
#include "../utils/magnetic_field_types.hpp"
#include "../utils/utils.hpp"
#include "traccc/cuda/seeding/hough_transform_seeding_algorithm.hpp"
#include "traccc/device/global_index.hpp"

namespace traccc::cuda {
namespace kernels {
   
   __global__ 
   void define_plane_edges(const traccc::hough_seeding_config config,
                           const unsigned int n_cols,
                           const unsigned int n_rows,
                           vecmem::data::vector_view<scalar> edges_x,
                           vecmem::data::vector_view<scalar> edges_y) {
      scalar steps_x = (config.max_range_x - config.min_range_x) / n_cols;
      scalar steps_y = (config.max_range_y - config.min_range_y) / n_rows;
      device::define_plane_edges(details::global_index1(),
              n_cols, n_rows,
              steps_x, config.min_range_x,
		      steps_y, config.min_range_y,
              edges_x, edges_y); 
   }
   
   template <typename bfield_t>
   __global__
   void precompute_spacepoint_variables(const traccc::hough_seeding_config config,
           const bfield_t bfield,
		   const edm::measurement_collection::const_view measurements,
		   const edm::spacepoint_collection::const_view::size_type n_spacepoints,
		   const edm::spacepoint_collection::const_view spacepoints,
		   vecmem::data::vector_view<unsigned long long> spacepoint_wedge_association,
		   vecmem::data::vector_view<unsigned long long> spacepoint_detector_id,
		   vecmem::data::vector_view<scalar> spacepoint_r,
		   vecmem::data::vector_view<scalar> spacepoint_phi,
           vecmem::data::vector_view<scalar> spacepoint_bfield_z) {
      device::precompute_spacepoint_variables<bfield_t>(details::global_index1(), bfield,
		      measurements, n_spacepoints, spacepoints,
		      config.interaction_center, config.interaction_window, config.eta_width, config.module_veto,
		      spacepoint_wedge_association, spacepoint_detector_id,
		      spacepoint_r, spacepoint_phi,
              spacepoint_bfield_z);
   }

   template <typename T>
   __global__
   void fill_hough_plane(const traccc::hough_seeding_config config,
		   const edm::spacepoint_collection::const_view::size_type n_spacepoints,
		   const vecmem::data::vector_view<const unsigned long long> spacepoint_wedge_association,
           const vecmem::data::vector_view<const unsigned long long> spacepoint_detector_id,
		   const vecmem::data::vector_view<const scalar> spacepoint_r,
           const vecmem::data::vector_view<const scalar> spacepoint_phi,
           const vecmem::data::vector_view<const scalar> spacepoint_bfield_z,
           const vecmem::data::vector_view<const scalar> edges_x,
           const vecmem::data::vector_view<const scalar> edges_y,
           vecmem::data::vector_view<unsigned long long> bins) {
      device::fill_hough_plane<T>(details::global_index1(),
              config.wedges.size(),
              config.n_cols,
              config.n_rows,
		      n_spacepoints,
		      spacepoint_wedge_association,
              spacepoint_detector_id,
		      spacepoint_r,
              spacepoint_phi,
              spacepoint_bfield_z,
              edges_x,
              edges_y,
		      bins);
   }

   __global__
   void count_bits(const unsigned long long n_entries,
            const vecmem::data::vector_view<const unsigned long long> bits,
            vecmem::data::vector_view<unsigned short> entries) {
      device::count_bits(details::global_index1(),
            n_entries, bits, entries);
   }

   __global__
   void count_maxima(const traccc::hough_seeding_config config,
                     const vecmem::data::vector_view<const unsigned short> houghplanes,
                     vecmem::data::vector_view<unsigned int> count_vector) {
      device::count_maxima(details::global_index1(),
                           config.wedges.size(),
                           config.n_cols,
                           config.n_rows, 
                           houghplanes,
                           config.min_contributions,
                           count_vector);
   }

   __global__
   void find_maxima(const traccc::hough_seeding_config config,
                    const vecmem::data::vector_view<const unsigned short> houghplanes,
                    vecmem::data::jagged_vector_view<traccc::device::global_index_t> global_indexes) {
      device::find_maxima(details::global_index1(),
                          config.wedges.size(),
                          config.n_cols,
                          config.n_rows,
                          houghplanes,
                          config.min_contributions,
                          global_indexes);
   }

   __global__
   void take_pictures(const traccc::hough_seeding_config config,
                      const vecmem::data::vector_view<const unsigned short> houghplanes,
                      const vecmem::data::vector_view<const scalar> edges_x,
                      const vecmem::data::jagged_vector_view<const traccc::device::global_index_t> global_indexes,
                      const unsigned int n_maxima,
                      vecmem::data::jagged_vector_view<scalar> images,
                      vecmem::data::vector_view<scalar> maxima_x) {
      device::take_pictures(blockIdx.x,
            threadIdx.y, // cols
            threadIdx.x, // rows
            config.n_cols,
            config.n_rows,
            n_maxima,
            houghplanes, 
            edges_x,
            global_indexes, 
            images, 
            maxima_x);
   }

   __global__
   void evaluate_nn(const unsigned int n_images,
                    const vecmem::data::jagged_vector_view<const scalar> images,
                    const vecmem::data::vector_view<const scalar> maxima_x,
                    vecmem::data::vector_view<unsigned int> nn_outcome) {
        device::evaluate_nn(details::global_index1(), n_images, images, maxima_x, nn_outcome);
    }

   __global__
   void count_contributing_hits(const traccc::hough_seeding_config config,
                                const unsigned int n_images,
                                const vecmem::data::vector_view<const scalar> edges_x,
                                const vecmem::data::vector_view<const scalar> edges_y,
                                const vecmem::data::vector_view<const scalar> spacepoint_r,
                                const vecmem::data::vector_view<const scalar> spacepoint_phi,
                                const vecmem::data::vector_view<const scalar> spacepoint_bfield_z,
                                const vecmem::data::vector_view<const unsigned long long> spacepoint_wedge_association,
                                const vecmem::data::jagged_vector_view<const traccc::device::global_index_t> global_index,
                                const vecmem::data::vector_view<const unsigned int> nn_outcome,
                                vecmem::data::vector_view<unsigned int> count_vector) {
        device::count_contributing_hits(details::global_index1(),
                    config.n_cols,
                    config.n_rows,
                    n_images,
                    edges_x,
                    edges_y,
                    spacepoint_r,
                    spacepoint_phi,
                    spacepoint_bfield_z,
                    spacepoint_wedge_association,
                    global_index,
                    nn_outcome,
                    count_vector);
   }

   __global__
   void find_contributing_hits(const traccc::hough_seeding_config config,
                               const unsigned int n_images,
                               const vecmem::data::vector_view<const scalar> edges_x,
                               const vecmem::data::vector_view<const scalar> edges_y,
                               const vecmem::data::vector_view<const scalar> spacepoint_r,
                               const vecmem::data::vector_view<const scalar> spacepoint_phi,
                               const vecmem::data::vector_view<const scalar> spacepoint_bfield_z,
                               const vecmem::data::vector_view<const unsigned long long> spacepoint_wedge_association,
                               const vecmem::data::jagged_vector_view<const traccc::device::global_index_t> global_index,
                               const vecmem::data::vector_view<const unsigned int> nn_outcome,
                               vecmem::data::jagged_vector_view<traccc::device::global_index_t> hits) {
      device::find_contributing_hits(details::global_index1(),
                    config.n_cols,
                    config.n_rows,
                    n_images,
                    edges_x,
                    edges_y,
                    spacepoint_r,
                    spacepoint_phi,
                    spacepoint_bfield_z,
                    spacepoint_wedge_association,
                    global_index,
                    nn_outcome,
                    hits);
  }

  __global__
  void make_seeds(const unsigned int n_images,
                  const edm::spacepoint_collection::const_view spacepoints,
                  const vecmem::data::vector_view<const scalar> spacepoint_r,
                  const vecmem::data::vector_view<const unsigned int> nn_outcome,
                  const vecmem::data::jagged_vector_view<const traccc::device::global_index_t> hits,
                  vecmem::data::jagged_vector_view<scalar> cot_theta_matrix) {
        device::make_seeds(details::global_index1(), n_images,
                           spacepoints, spacepoint_r,
                           nn_outcome, hits,
                           cot_theta_matrix);
  }
} // namespace kernels

hough_transform_seeding_algorithm::hough_transform_seeding_algorithm(
    const traccc::memory_resource& mr,
    const vecmem::copy& copy, const stream_wrapper& str,
    const hough_seeding_config& seeding_config,
    std::unique_ptr<const Logger> logger, await_function_type await_func)
    : device::hough_transform_seeding_algorithm(mr, copy, seeding_config,
                                        std::move(logger)),
      cuda::algorithm_base{str, std::move(await_func)} {}

void hough_transform_seeding_algorithm::define_plane_edges_kernel(const define_plane_edges_kernel_payload& payload) const {
   unsigned int maxValues = std::max(payload.config.n_cols, payload.config.n_rows);
   const unsigned int n_threads = warp_size() * 4;
   const unsigned int n_blocks = (maxValues  + n_threads - 1) / n_threads;
   kernels::define_plane_edges
	   <<<n_blocks, n_threads, 0, details::get_stream(stream())>>>(payload.config, 
                                                                   payload.config.n_cols,
                                                                   payload.config.n_rows,
                                                                   payload.edges_x,
                                                                   payload.edges_y);
   TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void hough_transform_seeding_algorithm::precompute_spacepoint_variables_kernel(const precompute_spacepoint_variables_payload& payload) const {
   const unsigned int n_threads = warp_size() * 8;
   const unsigned int n_blocks = (payload.n_spacepoints  + n_threads - 1) / n_threads;
   magnetic_field_visitor<bfield_type_list<scalar>>(
     payload.bfield, [&]<typename bfield_view_t>(const bfield_view_t& bfield) {
        kernels::precompute_spacepoint_variables<bfield_view_t>
		    <<<n_blocks, n_threads, 0, details::get_stream(stream())>>>(payload.config, bfield,
				payload.measurements, payload.n_spacepoints, payload.spacepoints,
				payload.spacepoint_wedge_association,
				payload.spacepoint_detector_id,
				payload.spacepoint_r,
				payload.spacepoint_phi,
                payload.spacepoint_bfield_z);
   });
   TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void hough_transform_seeding_algorithm::fill_hough_plane_kernel(const fill_hough_planes_kernel_payload& payload) const {
  const unsigned int n_threads = warp_size() * 8;
  const unsigned int n_blocks = (payload.n_spacepoints + n_threads - 1) / n_threads;
  kernels::fill_hough_plane<unsigned short>
            <<<n_blocks, n_threads, 0, details::get_stream(stream())>>>(payload.config,
			    payload.n_spacepoints,
			    payload.spacepoint_wedge_association,
			    payload.spacepoint_detector_id,
			    payload.spacepoint_r,
                payload.spacepoint_phi,
                payload.spacepoint_bfield_z,
                payload.edges_x,
                payload.edges_y,
			    payload.bins);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void hough_transform_seeding_algorithm::count_bits_kernel(const count_bits_kernel_payload& payload) const {
  const unsigned int n_threads = warp_size() * 8;
  const  unsigned int n_blocks = (payload.n_entries + n_threads - 1) / n_threads;
  kernels::count_bits<<<n_blocks, n_threads, 0, details::get_stream(stream())>>>(payload.n_entries,
        payload.bits,
        payload.values);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void hough_transform_seeding_algorithm::count_maxima_kernel(const count_maxima_payload& payload) const {
  const unsigned int n_threads = warp_size() * 8;
  const unsigned int n_blocks = (payload.config.n_cols * payload.config.n_rows + n_threads - 1) / n_threads;
  kernels::count_maxima
        <<<n_blocks, n_threads, 0, details::get_stream(stream())>>>(payload.config, payload.houghplanes, payload.count_vector);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void hough_transform_seeding_algorithm::find_maxima_kernel(const find_maxima_payload& payload) const {
  const unsigned int n_threads = warp_size() * 8;
  const unsigned int n_blocks = (payload.config.n_cols * payload.config.n_rows + n_threads - 1) / n_threads;
  kernels::find_maxima
        <<<n_blocks, n_threads, 0, details::get_stream(stream())>>>(payload.config, payload.houghplanes, payload.global_index);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void hough_transform_seeding_algorithm::take_picture_kernel(const take_picture_payload& payload) const {
   const dim3 n_threads(32, 32);
   const unsigned int n_blocks = payload.n_maxima;
   kernels::take_pictures
        <<<n_blocks, n_threads, 0, details::get_stream(stream())>>>(payload.config,
                                                                    payload.houghplanes,
                                                                    payload.edges_x,
                                                                    payload.global_index, 
                                                                    payload.n_maxima, 
                                                                    payload.images,
                                                                    payload.maxima_x);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void hough_transform_seeding_algorithm::evaluate_nn_kernel(const evaluate_nn_payload& payload) const {
  const unsigned int n_threads = warp_size() * 8;
  const unsigned int n_blocks = (payload.n_images + n_threads - 1) / n_threads;
  kernels::evaluate_nn
        <<<n_blocks, n_threads, 0, details::get_stream(stream())>>>(payload.n_images,
                                                                    payload.images,
                                                                    payload.maxima_x,
                                                                    payload.nn_outcome);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void hough_transform_seeding_algorithm::count_contributing_hits_kernel(const count_contributing_hits_payload& payload) const {
  const unsigned int n_threads = warp_size() * 8;
  const unsigned int n_blocks = (payload.n_images + n_threads - 1) / n_threads;
  kernels::count_contributing_hits
            <<<n_blocks, n_threads, 0, details::get_stream(stream())>>>(payload.config,
                                                                    payload.n_images,
                                                                    payload.edges_x,
                                                                    payload.edges_y,
                                                                    payload.spacepoint_r,
                                                                    payload.spacepoint_phi,
                                                                    payload.spacepoint_bfield_z,
                                                                    payload.spacepoint_wedge_association,
                                                                    payload.global_index,
                                                                    payload.nn_outcome,
                                                                    payload.count_vector);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void hough_transform_seeding_algorithm::find_contributing_hits_kernel(const find_contributing_hits_payload& payload) const {
  const unsigned int n_threads = warp_size() * 4;
  const unsigned int n_blocks = (payload.n_images + n_threads - 1) / n_threads;
  kernels::find_contributing_hits
            <<<n_blocks, n_threads, 0, details::get_stream(stream())>>>(payload.config,
                                                                    payload.n_images,
                                                                    payload.edges_x,
                                                                    payload.edges_y,
                                                                    payload.spacepoint_r,
                                                                    payload.spacepoint_phi,
                                                                    payload.spacepoint_bfield_z,
                                                                    payload.spacepoint_wedge_association,
                                                                    payload.global_index,
                                                                    payload.nn_outcome,
                                                                    payload.hits);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void hough_transform_seeding_algorithm::make_seeds_kernel(const make_seeds_payload& payload) const {
    const unsigned int n_threads = warp_size() * 4;
    const unsigned int n_blocks = (payload.n_images + n_threads - 1) / n_threads;
    kernels::make_seeds<<<n_blocks, n_threads, 0, details::get_stream(stream())>>>(payload.n_images,
                                                                                   payload.spacepoints,
                                                                                   payload.spacepoint_r,
                                                                                   payload.nn_outcome,                                                                                   
                                                                                   payload.hits,
                                                                                   payload.cot_theta_matrix);
}

void hough_transform_seeding_algorithm::synchronize() const {
  stream().synchronize();
}

}  // namespace traccc::cuda
