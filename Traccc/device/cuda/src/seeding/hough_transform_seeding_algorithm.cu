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
      float steps_x = (config.max_range_x - config.min_range_x) / n_cols;
      float steps_y = (config.max_range_y - config.min_range_y) / n_rows;
      device::define_plane_edges(details::global_index1(),
              edges_x, edges_y,
		      steps_x, config.min_range_x,
		      steps_y, config.min_range_y); 
   }
   
   __global__
   void precompute_spacepoint_variables(const traccc::hough_seeding_config config,
		   const edm::measurement_collection::const_view measurements,
		   const edm::spacepoint_collection::const_view::size_type n_spacepoints,
		   const edm::spacepoint_collection::const_view spacepoints,
		   vecmem::data::vector_view<unsigned long long> spacepoint_wedge_association,
		   vecmem::data::vector_view<unsigned long long> spacepoint_detector_id,
		   vecmem::data::vector_view<scalar> spacepoint_r,
		   vecmem::data::vector_view<scalar> spacepoint_phi) {
      device::precompute_spacepoint_variables(details::global_index1(),
		      measurements, n_spacepoints, spacepoints,
		      config.interaction_center, config.interaction_window, config.eta_width, config.module_veto,
		      spacepoint_wedge_association, spacepoint_detector_id,
		      spacepoint_r, spacepoint_phi);
   }

   template <typename bfield_t,
             typename T,
             unsigned int N_COLS,
             unsigned int N_ROWS>
   __global__
   void fill_hough_plane(const traccc::hough_seeding_config config,
		   const bfield_t bfield,
		   const edm::spacepoint_collection::const_view::size_type n_spacepoints,
		   const edm::spacepoint_collection::const_view spacepoints,
		   const vecmem::data::vector_view<unsigned long long> spacepoint_wedge_association,
           const vecmem::data::vector_view<unsigned long long> spacepoint_detector_id,
		   const vecmem::data::vector_view<scalar> spacepoint_r,
           const vecmem::data::vector_view<scalar> spacepoint_phi,
           const vecmem::data::vector_view<scalar> edges_x,
           const vecmem::data::vector_view<scalar> edges_y,
		   typename edm::houghplane_collection<T>::view houghplanes) {
      device::fill_hough_plane<bfield_t, T, N_COLS, N_ROWS>(details::global_index1(), bfield,
		      n_spacepoints,
              spacepoints,
		      spacepoint_wedge_association,
              spacepoint_detector_id,
		      spacepoint_r,
              spacepoint_phi,
              edges_x,
              edges_y,
		      houghplanes);
   }

   template <typename T,
             unsigned int N_COLS,
             unsigned int N_ROWS>
   __global__
   void count_maxima(const typename edm::houghplane_collection<T>::view houghplanes,
                     vecmem::data::vector_view<unsigned int> count_vector) {
      device::count_maxima<T, N_COLS, N_ROWS>(details::global_index1(), houghplanes, count_vector);
   }

   template <typename T,
             unsigned int N_COLS,
             unsigned int N_ROWS>
   __global__
   void find_maxima(const typename edm::houghplane_collection<T>::view houghplanes,
                    vecmem::data::jagged_vector_view<traccc::device::global_index_t> global_indexes) {
      device::find_maxima<T, N_COLS, N_ROWS>(details::global_index1(), houghplanes, global_indexes);
   }

   template <typename T,
             unsigned int N_COLS,
             unsigned int N_ROWS>
   __global__
   void take_pictures(const typename edm::houghplane_collection<T>::const_view houghplanes,
                      const vecmem::data::vector_view<scalar> edges_x,
                      const vecmem::data::jagged_vector_view<traccc::device::global_index_t> global_indexes,
                      const unsigned int n_maxima,
                      vecmem::data::jagged_vector_view<scalar> images,
                      vecmem::data::vector_view<scalar> maxima_x) {
      device::take_pictures<T, N_COLS, N_ROWS>(blockIdx.x,
            threadIdx.y, // cols
            threadIdx.x, // rows
            houghplanes, edges_x,
            global_indexes, n_maxima, images, maxima_x);
   }

   __global__
   void evaluate_nn(const unsigned int n_images,
                    const vecmem::data::jagged_vector_view<scalar> images,
                    const vecmem::data::vector_view<scalar> maxima_x,
                    vecmem::data::vector_view<unsigned int> nn_outcome) {
        device::evaluate_nn(details::global_index1(), n_images, images, maxima_x, nn_outcome);
    }

   template <typename bfield_t,
             typename T,
             unsigned int N_COLS,
             unsigned int N_ROWS>
   __global__
   void count_contributing_hits(const unsigned int n_images,
                                const bfield_t bfield,
                                const typename edm::houghplane_collection<T>::const_view houghplanes,
                                const vecmem::data::vector_view<scalar> edges_x,
                                const vecmem::data::vector_view<scalar> edges_y,
                                const edm::spacepoint_collection::const_view spacepoints,
                                const vecmem::data::vector_view<scalar> spacepoint_r,
                                const vecmem::data::vector_view<scalar> spacepoint_phi,
                                const vecmem::data::vector_view<unsigned long long> spacepoint_wedge_association,
                                const vecmem::data::jagged_vector_view<traccc::device::global_index_t> global_index,
                                const vecmem::data::vector_view<unsigned int> nn_outcome,
                                vecmem::data::vector_view<unsigned int> count_vector) {
        device::count_contributing_hits<bfield_t, T, N_COLS, N_ROWS>(details::global_index1(),
                    n_images,
                    bfield,
                    houghplanes,
                    edges_x,
                    edges_y,
                    spacepoints,
                    spacepoint_r,
                    spacepoint_phi,
                    spacepoint_wedge_association,
                    global_index,
                    nn_outcome,
                    count_vector);
   }

   template <typename bfield_t,
             typename T,
             unsigned int N_COLS,
             unsigned int N_ROWS>
   __global__
   void find_contributing_hits(const unsigned int n_images,
                                const bfield_t bfield,
                                const typename edm::houghplane_collection<T>::const_view houghplanes,
                                const vecmem::data::vector_view<scalar> edges_x,
                                const vecmem::data::vector_view<scalar> edges_y,
                                const edm::spacepoint_collection::const_view spacepoints,
                                const vecmem::data::vector_view<scalar> spacepoint_r,
                                const vecmem::data::vector_view<scalar> spacepoint_phi,
                                const vecmem::data::vector_view<unsigned long long> spacepoint_wedge_association,
                                const vecmem::data::jagged_vector_view<traccc::device::global_index_t> global_index,
                                const vecmem::data::vector_view<unsigned int> nn_outcome,
                                vecmem::data::jagged_vector_view<traccc::device::global_index_t> hits) {
      device::find_contributing_hits<bfield_t, T, N_COLS, N_ROWS>(details::global_index1(),
                    n_images,
                    bfield,
                    houghplanes,
                    edges_x,
                    edges_y,
                    spacepoints,
                    spacepoint_r,
                    spacepoint_phi,
                    spacepoint_wedge_association,
                    global_index,
                    nn_outcome,
                    hits);
  }

  __global__
  void make_seeds(const unsigned int n_images,
                  const edm::spacepoint_collection::const_view spacepoints,
                  const vecmem::data::vector_view<scalar> spacepoint_r,
                  const vecmem::data::vector_view<unsigned int> nn_outcome,
                  const vecmem::data::jagged_vector_view<traccc::device::global_index_t> hits,
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
   unsigned int maxValues = std::max(device::N_COLS, device::N_ROWS);
   const unsigned int n_threads = warp_size() * 4;
   const unsigned int n_blocks = (maxValues  + n_threads - 1) / n_threads;
   kernels::define_plane_edges
	   <<<n_blocks, n_threads, 0, details::get_stream(stream())>>>(payload.config, 
                                                                   device::N_COLS,
                                                                   device::N_ROWS,
                                                                   payload.edges_x,
                                                                   payload.edges_y);
   TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void hough_transform_seeding_algorithm::precompute_spacepoint_variables_kernel(const precompute_spacepoint_variables_payload& payload) const {
   const unsigned int n_threads = warp_size() * 32;
   const unsigned int n_blocks = (payload.n_spacepoints  + n_threads - 1) / n_threads;
   kernels::precompute_spacepoint_variables
		<<<n_blocks, n_threads, 0, details::get_stream(stream())>>>(payload.config,
				payload.measurements, payload.n_spacepoints, payload.spacepoints,
				payload.spacepoint_wedge_association,
				payload.spacepoint_detector_id,
				payload.spacepoint_r,
				payload.spacepoint_phi);
   TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void hough_transform_seeding_algorithm::fill_hough_plane_kernel(const fill_hough_planes_kernel_payload& payload) const {
  const unsigned int n_threads = warp_size() * 8;
  const unsigned int n_blocks = (device::N_COLS * device::N_ROWS + n_threads - 1) / n_threads;
  magnetic_field_visitor<bfield_type_list<scalar>>(
  	payload.bfield, [&]<typename bfield_view_t>(const bfield_view_t& bfield) {
		kernels::fill_hough_plane<bfield_view_t, unsigned short, device::N_COLS, device::N_ROWS>
            <<<n_blocks, n_threads, 0, details::get_stream(stream())>>>(payload.config, bfield,
			    payload.n_spacepoints, payload.spacepoints,
			    payload.spacepoint_wedge_association,
			    payload.spacepoint_detector_id,
			    payload.spacepoint_r,
                payload.spacepoint_phi,
                payload.edges_x,
                payload.edges_y,
			    payload.houghplanes);
  });
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void hough_transform_seeding_algorithm::count_maxima_kernel(const count_maxima_payload& payload) const {
  const unsigned int n_threads = warp_size() * 8;
  const unsigned int n_blocks = (device::N_COLS * device::N_ROWS + n_threads - 1) / n_threads;
  kernels::count_maxima<unsigned short, device::N_COLS, device::N_ROWS>
        <<<n_blocks, n_threads, 0, details::get_stream(stream())>>>(payload.houghplanes, payload.count_vector);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void hough_transform_seeding_algorithm::find_maxima_kernel(const find_maxima_payload& payload) const {
  const unsigned int n_threads = warp_size() * 8;
  const unsigned int n_blocks = (device::N_COLS * device::N_ROWS + n_threads - 1) / n_threads;
  kernels::find_maxima<unsigned short, device::N_COLS, device::N_ROWS>
        <<<n_blocks, n_threads, 0, details::get_stream(stream())>>>(payload.houghplanes, payload.global_index);
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void hough_transform_seeding_algorithm::take_picture_kernel(const take_picture_payload& payload) const {
   const dim3 n_threads(32, 32);
   const unsigned int n_blocks = payload.n_maxima;
   kernels::take_pictures<unsigned short, device::N_COLS, device::N_ROWS>
        <<<n_blocks, n_threads, 0, details::get_stream(stream())>>>(payload.houghplanes,
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
  magnetic_field_visitor<bfield_type_list<scalar>>(
    payload.bfield, [&]<typename bfield_view_t>(const bfield_view_t& bfield) {
        kernels::count_contributing_hits<bfield_view_t, unsigned short, device::N_COLS, device::N_ROWS>
            <<<n_blocks, n_threads, 0, details::get_stream(stream())>>>(payload.n_images,
                                                                    bfield,
                                                                    payload.houghplanes,
                                                                    payload.edges_x,
                                                                    payload.edges_y,
                                                                    payload.spacepoints,
                                                                    payload.spacepoint_r,
                                                                    payload.spacepoint_phi,
                                                                    payload.spacepoint_wedge_association,
                                                                    payload.global_index,
                                                                    payload.nn_outcome,
                                                                    payload.count_vector);
  });
  TRACCC_CUDA_ERROR_CHECK(cudaGetLastError());
}

void hough_transform_seeding_algorithm::find_contributing_hits_kernel(const find_contributing_hits_payload& payload) const {
  const unsigned int n_threads = warp_size() * 4;
  const unsigned int n_blocks = (payload.n_images + n_threads - 1) / n_threads;
  magnetic_field_visitor<bfield_type_list<scalar>>(
    payload.bfield, [&]<typename bfield_view_t>(const bfield_view_t& bfield) {
        kernels::find_contributing_hits<bfield_view_t, unsigned short, device::N_COLS, device::N_ROWS>
            <<<n_blocks, n_threads, 0, details::get_stream(stream())>>>(payload.n_images,
                                                                    bfield,
                                                                    payload.houghplanes,
                                                                    payload.edges_x,
                                                                    payload.edges_y,
                                                                    payload.spacepoints,
                                                                    payload.spacepoint_r,
                                                                    payload.spacepoint_phi,
                                                                    payload.spacepoint_wedge_association,
                                                                    payload.global_index,
                                                                    payload.nn_outcome,
                                                                    payload.hits);
  });
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
