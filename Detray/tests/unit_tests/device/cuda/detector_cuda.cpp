// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2020-2024 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Detray test include(s)
#include "detector_cuda_kernel.hpp"
#include "detray/core/detail/alignment.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/test/common/build_toy_detector.hpp"
#include "detray/test/framework/assert.hpp"
#include "detray/test/framework/types.hpp"

// Vecmem include(s)
#include <vecmem/memory/cuda/device_memory_resource.hpp>
#include <vecmem/memory/cuda/managed_memory_resource.hpp>
#include <vecmem/memory/host_memory_resource.hpp>
#include <vecmem/utils/cuda/async_copy.hpp>
#include <vecmem/utils/cuda/copy.hpp>

// Google Test include(s)
#include <gtest/gtest.h>

// System include(s)
#include <limits>

using namespace detray;

TEST(detector_cuda, detector) {
  // memory resource
  vecmem::cuda::managed_memory_resource mng_mr;

  // create toy geometry
  auto [toy_det, names] = build_toy_detector<test::algebra>(mng_mr);

  auto ctx0 = typename detector_host_t::geometry_context();

  // host objects
  auto& volumes_host = toy_det.volumes();
  auto& surfaces_host = toy_det.surfaces();
  auto& transforms_host = toy_det.transform_store();
  auto& masks_host = toy_det.mask_store();
  auto& discs_host = masks_host.get<disc_id>();
  auto& cylinders_host = masks_host.get<cylinder_id>();
  auto& rectangles_host = masks_host.get<rectangle_id>();

  // copied output from device side
  vecmem::vector<det_volume_t> volumes_device(volumes_host.size(), &mng_mr);
  vecmem::vector<det_surface_t> surfaces_device(surfaces_host.size(), &mng_mr);
  vecmem::vector<transform_t> transforms_device(transforms_host.size(),
                                                &mng_mr);
  vecmem::vector<rectangle_t> rectangles_device(rectangles_host.size(),
                                                &mng_mr);
  vecmem::vector<disc_t> discs_device(discs_host.size(), &mng_mr);
  vecmem::vector<cylinder_t> cylinders_device(cylinders_host.size(), &mng_mr);

  // get data object for toy detector
  auto toy_det_data = detray::get_data(toy_det);

  // get data object for device outputs
  auto volumes_data = vecmem::get_data(volumes_device);
  auto surfaces_data = vecmem::get_data(surfaces_device);
  auto transforms_data = vecmem::get_data(transforms_device);
  auto rectangles_data = vecmem::get_data(rectangles_device);
  auto discs_data = vecmem::get_data(discs_device);
  auto cylinders_data = vecmem::get_data(cylinders_device);

  // run the test code to copy the objects
  detector_test(toy_det_data, volumes_data, surfaces_data, transforms_data,
                rectangles_data, discs_data, cylinders_data);

  // check if the same volume objects are copied
  for (unsigned int i = 0u; i < volumes_host.size(); i++) {
    EXPECT_EQ(volumes_host[i] == volumes_device[i], true);
  }

  // check if the same surface objects are copied
  for (unsigned int i = 0u; i < surfaces_host.size(); i++) {
    EXPECT_EQ(surfaces_device[i] == surfaces_host[i], true);
  }

  // check if the same transform objects are copied
  for (unsigned int i = 0u; i < transforms_host.size(ctx0); i++) {
    EXPECT_EQ(transforms_host.at(i, ctx0) == transforms_device[i], true);
  }

  // check if the same masks are copied
  for (unsigned int i = 0u; i < rectangles_host.size(); i++) {
    EXPECT_EQ(rectangles_host[i] == rectangles_device[i], true);
  }

  for (unsigned int i = 0u; i < discs_host.size(); i++) {
    EXPECT_EQ(discs_host[i] == discs_device[i], true);
  }

  for (unsigned int i = 0u; i < cylinders_host.size(); i++) {
    EXPECT_EQ(cylinders_host[i] == cylinders_device[i], true);
  }
}

TEST(detector_cuda, detector_alignment) {
  // a few typedefs
  using test_algebra = test::algebra;
  using scalar = dscalar<test_algebra>;
  using point3 = dpoint3D<test_algebra>;

  // memory resources
  vecmem::host_memory_resource host_mr;
  vecmem::cuda::device_memory_resource dev_mr;
  vecmem::cuda::managed_memory_resource mng_mr;

  // helper object for performing memory copies to CUDA devices
  vecmem::cuda::copy cuda_cpy;

  // create toy geometry in host memory
  auto [det_host, names_host] = build_toy_detector<test_algebra>(host_mr);

  // copy static detector data (including the initial set of transforms) to
  // the device
  // use synchronous copy and fixed size buffers
  auto det_buff_static = detray::get_buffer(det_host, dev_mr, cuda_cpy);

  // ---------- construct an "aligned" transform store ---------

  // build a vector of aligned transforms on the host
  // for populating this vector take all transforms of the detector
  // and shift them by the same translation
  typename detector_host_t::transform_container tf_store_aligned_host;

  point3 shift{.1f * unit<scalar>::mm, .2f * unit<scalar>::mm,
               .3f * unit<scalar>::mm};

  tf_store_aligned_host.reserve(
      det_host.transform_store().size(),
      typename decltype(det_host)::transform_container::context_type{});

  for (const auto& tf : det_host.transform_store()) {
    point3 shifted = tf.translation() + shift;
    tf_store_aligned_host.push_back(
        transform_t{shifted, tf.x(), tf.y(), tf.z()});
  }

  // copy the vector of aligned transforms to the device
  // again, use synchronous copy and fixed size buffers
  auto tf_buff_aligned =
      get_buffer(tf_store_aligned_host, dev_mr, cuda_cpy, copy::sync,
                 vecmem::data::buffer_type::fixed_size);

  // Get the view of the aligned detector using the vector of aligned
  // transforms and the static part of the detector copied to the device
  // earlier
  auto detector_view_aligned =
      detail::misaligned_detector_view<detector_host_t>(det_buff_static,
                                                        tf_buff_aligned);
  // Get the view of the static detector
  auto detector_view_static = detray::get_data(det_buff_static);

  // make two vectors for surface transforms copied from device side
  vecmem::vector<transform_t> surfacexf_device_static(
      det_host.surfaces().size(), &mng_mr);
  vecmem::vector<transform_t> surfacexf_device_aligned(
      det_host.surfaces().size(), &mng_mr);
  // views of the above vectors
  auto surfacexf_data_static = vecmem::get_data(surfacexf_device_static);
  auto surfacexf_data_aligned = vecmem::get_data(surfacexf_device_aligned);

  // run the test code to extract the surface transforms for the static
  // and misaligned detector views and to store them into the vectors
  detector_alignment_test(detector_view_static, detector_view_aligned,
                          surfacexf_data_static, surfacexf_data_aligned);

  // check that the relevant transforms have been properly shifted
  for (unsigned int i = 0u; i < surfacexf_device_static.size(); i++) {
    auto translation_static = surfacexf_device_static[i].translation();
    auto translation_aligned = surfacexf_device_aligned[i].translation();
    auto translation_diff = translation_aligned - translation_static;
    EXPECT_POINT3_NEAR(translation_diff, shift, 1e-4);
  }
}
