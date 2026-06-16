// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "detray/definitions/detail/cuda_definitions.hpp"
#include "detray/geometry/tracking_surface.hpp"

// Detray test include(s)
#include "detector_cuda_kernel.hpp"

namespace detray {

// cuda kernel to copy sub-detector objects
__global__ void detector_test_kernel(
    typename detector_host_t::view_type det_data,
    vecmem::data::vector_view<det_volume_t> volumes_data,
    vecmem::data::vector_view<det_surface_t> surfaces_data,
    vecmem::data::vector_view<transform_t> transforms_data,
    vecmem::data::vector_view<rectangle_t> rectangles_data,
    vecmem::data::vector_view<disc_t> discs_data,
    vecmem::data::vector_view<cylinder_t> cylinders_data) {
  // convert toy detector_data into detector w/ device vectors
  detector_device_t det_device(det_data);

  // convert subdetector data objects into objects w/ device vectors
  vecmem::device_vector<det_volume_t> volumes_device(volumes_data);
  vecmem::device_vector<det_surface_t> surfaces_device(surfaces_data);
  vecmem::device_vector<transform_t> transforms_device(transforms_data);
  vecmem::device_vector<rectangle_t> rectangles_device(rectangles_data);
  vecmem::device_vector<disc_t> discs_device(discs_data);
  vecmem::device_vector<cylinder_t> cylinders_device(cylinders_data);

  // copy objects - volume
  for (unsigned int i = 0u; i < det_device.volumes().size(); i++) {
    volumes_device[i] = det_device.volumes()[i];
  }

  // copy objects - surfaces
  for (unsigned int i = 0u; i < det_device.surfaces().size(); i++) {
    surfaces_device[i] = det_device.surfaces()[i];
  }

  // copy objects - transforms
  auto& trfs = det_device.transform_store();
  auto ctx = typename detector_host_t::geometry_context{};
  for (unsigned int i = 0u; i < trfs.size(ctx); i++) {
    transforms_device[i] = trfs.at(i, ctx);
  }

  // copy objects - masks
  auto& masks = det_device.mask_store();
  auto& rectangles =
      masks.template get<detector_host_t::masks::id::e_rectangle2D>();
  for (unsigned int i = 0u; i < rectangles.size(); i++) {
    rectangles_device[i] = rectangles[i];
  }

  auto& discs = masks.template get<detector_host_t::masks::id::e_ring2D>();
  for (unsigned int i = 0u; i < discs.size(); i++) {
    discs_device[i] = discs[i];
  }

  auto& cylinders =
      masks.template get<detector_host_t::masks::id::e_concentric_cylinder2D>();
  for (unsigned int i = 0u; i < cylinders.size(); i++) {
    cylinders_device[i] = cylinders[i];
  }

  // print output test for surface finder
  /*auto& accel_device = det_device.accelerator_store();
  for (unsigned int i_s = 0u; i_s < accel_device.size(); i_s++) {
      auto& grid = accel_device[i_s];
      for (unsigned int i = 0u; i < grid.axis_p0().bins(); i++) {
          for (unsigned int j = 0u; j < grid.axis_p1().bins(); j++) {
              const auto& bin = grid.bin(i, j);
              for (auto& id : bin) {
                  // printf("%d \n", id);
              }
          }
      }
  }*/
}

/// implementation of the test function for detector
void detector_test(typename detector_host_t::view_type det_data,
                   vecmem::data::vector_view<det_volume_t> volumes_data,
                   vecmem::data::vector_view<det_surface_t> surfaces_data,
                   vecmem::data::vector_view<transform_t> transforms_data,
                   vecmem::data::vector_view<rectangle_t> rectangles_data,
                   vecmem::data::vector_view<disc_t> discs_data,
                   vecmem::data::vector_view<cylinder_t> cylinders_data) {
  constexpr int block_dim = 1u;
  constexpr int thread_dim = 1u;

  // run the test kernel
  detector_test_kernel<<<block_dim, thread_dim>>>(
      det_data, volumes_data, surfaces_data, transforms_data, rectangles_data,
      discs_data, cylinders_data);

  // cuda error check
  DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
  DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

// cuda kernel to extract surface transforms from two detector views - static
// and misaligned - and to copy them into vectors
__global__ void detector_alignment_test_kernel(
    typename detector_host_t::view_type det_data_static,
    typename detector_host_t::view_type det_data_aligned,
    vecmem::data::vector_view<transform_t> surfacexf_data_static,
    vecmem::data::vector_view<transform_t> surfacexf_data_aligned) {
  auto ctx = typename detector_host_t::geometry_context{};

  // two instances of device detectors
  detector_device_t det_device_static(det_data_static);
  detector_device_t det_device_aligned(det_data_aligned);

  // device vectors of surface transforms
  vecmem::device_vector<transform_t> surfacexf_device_static(
      surfacexf_data_static);
  vecmem::device_vector<transform_t> surfacexf_device_aligned(
      surfacexf_data_aligned);

  // copy surface transforms into relevant vectors
  for (unsigned int i = 0u; i < det_device_static.surfaces().size(); i++) {
    const auto sf =
        tracking_surface{det_device_static, det_device_static.surfaces()[i]};
    surfacexf_device_static[i] = sf.transform(ctx);
  }

  for (unsigned int i = 0u; i < det_device_aligned.surfaces().size(); i++) {
    const auto sf =
        tracking_surface{det_device_aligned, det_device_aligned.surfaces()[i]};
    surfacexf_device_aligned[i] = sf.transform(ctx);
  }
}

/// implementation of the alignment test function for detector
void detector_alignment_test(
    typename detector_host_t::view_type det_data_static,
    typename detector_host_t::view_type det_data_aligned,
    vecmem::data::vector_view<transform_t> surfacexf_data_static,
    vecmem::data::vector_view<transform_t> surfacexf_data_aligned) {
  constexpr int block_dim = 1u;
  constexpr int thread_dim = 1u;

  // run the test kernel
  detector_alignment_test_kernel<<<block_dim, thread_dim>>>(
      det_data_static, det_data_aligned, surfacexf_data_static,
      surfacexf_data_aligned);

  // cuda error check
  DETRAY_CUDA_ERROR_CHECK(cudaGetLastError());
  DETRAY_CUDA_ERROR_CHECK(cudaDeviceSynchronize());
}

}  // namespace detray
