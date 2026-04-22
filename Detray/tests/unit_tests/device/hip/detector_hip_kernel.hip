/** Detray library, part of the ACTS project (R&D line)
 *
 * (c) 2022-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include "detray/definitions/detail/hip_definitions.hpp"
#include "detray/geometry/tracking_surface.hpp"

// Detray test include(s)
#include "detector_hip_kernel.hpp"

namespace detray {

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

    auto& cylinders = masks.template get<
        detector_host_t::masks::id::e_concentric_cylinder2D>();
    for (unsigned int i = 0u; i < cylinders.size(); i++) {
        cylinders_device[i] = cylinders[i];
    }
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
    hipLaunchKernelGGL(detector_test_kernel, dim3(block_dim), dim3(thread_dim),
                       0, 0, det_data, volumes_data, surfaces_data,
                       transforms_data, rectangles_data, discs_data,
                       cylinders_data);
    // HIP error check
    DETRAY_HIP_ERROR_CHECK(hipGetLastError());
    DETRAY_HIP_ERROR_CHECK(hipDeviceSynchronize());
}
}  // namespace detray
