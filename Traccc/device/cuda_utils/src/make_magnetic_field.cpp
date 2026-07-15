/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/cuda/utils/make_magnetic_field.hpp"

// Project include(s).
#include "traccc/bfield/magnetic_field_types.hpp"

// Covfie include(s).
#include <covfie/core/backend/primitive/constant.hpp>
#include <covfie/core/backend/transformer/affine.hpp>
#include <covfie/core/backend/transformer/clamp.hpp>
#include <covfie/core/backend/transformer/linear.hpp>
#include <covfie/core/backend/transformer/strided.hpp>
#include <covfie/core/concepts.hpp>
#include <covfie/core/field.hpp>
#include <covfie/core/vector.hpp>
#include <covfie/cuda/backend/primitive/cuda_device_array.hpp>
#include <covfie/cuda/backend/primitive/cuda_texture.hpp>

// System include(s).
#include <stdexcept>

namespace traccc::cuda {

/// Inhomogeneous B-field backend type using CUDA global memory
template <typename scalar_t>
using inhom_global_bfield_backend_t =
    covfie::backend::affine<covfie::backend::linear<covfie::backend::clamp<
        covfie::backend::strided<covfie::vector::vector_d<std::size_t, 3>,
                                 covfie::backend::cuda_device_array<
                                     covfie::vector::vector_d<scalar_t, 3>>>>>>;
// Test that the type is a valid backend for a field
static_assert(
    covfie::concepts::field_backend<inhom_global_bfield_backend_t<float>>,
    "cuda::inhom_global_bfield_backend_t is not a valid field backend type");

/// Inhomogeneous B-field backend type using CUDA texture memory
using inhom_texture_bfield_backend_t = covfie::backend::affine<
    covfie::backend::cuda_texture<covfie::vector::vector_d<float, 3>,
                                  covfie::vector::vector_d<float, 3>>>;
// Test that the type is a valid backend for a field
static_assert(
    covfie::concepts::field_backend<inhom_texture_bfield_backend_t>,
    "cuda::inhom_texture_bfield_backend_t is not a valid field backend type");

magnetic_field make_magnetic_field(const magnetic_field& bfield,
                                   const magnetic_field_storage storage) {

    if (bfield.is<const_bfield_backend_t<float>>()) {
        return magnetic_field{covfie::field<const_bfield_backend_t<float>>{
            bfield.as_field<const_bfield_backend_t<float>>()}};
    } else if (bfield.is<const_bfield_backend_t<double>>()) {
        return magnetic_field{covfie::field<const_bfield_backend_t<double>>{
            bfield.as_field<const_bfield_backend_t<double>>()}};
    } else if (bfield.is<host::inhom_bfield_backend_t<float>>()) {
        // Convenience access to the Covfie field object.
        const auto& in_field =
            bfield.as_field<host::inhom_bfield_backend_t<float>>();
        // At single precision we can use either global or texture memory.
        if (storage == magnetic_field_storage::global_memory) {
            return magnetic_field{
                covfie::field<cuda::inhom_global_bfield_backend_t<float>>(
                    in_field)};
        } else if (storage == magnetic_field_storage::texture_memory) {
            return magnetic_field{
                covfie::field<cuda::inhom_texture_bfield_backend_t>(
                    covfie::make_parameter_pack(
                        in_field.backend().get_configuration(),
                        in_field.backend()
                            .get_backend()
                            .get_backend()
                            .get_backend()))};
        } else {
            throw std::invalid_argument(
                "Unsupported storage method chosen for inhomogeneous b-field");
        }
    } else if (bfield.is<host::inhom_bfield_backend_t<double>>()) {
        // At double precision we can only use global memory.
        if (storage == magnetic_field_storage::global_memory) {
            return magnetic_field{
                covfie::field<cuda::inhom_global_bfield_backend_t<double>>(
                    bfield.as_field<host::inhom_bfield_backend_t<double>>())};
        } else {
            throw std::invalid_argument(
                "Unsupported storage method chosen for inhomogeneous b-field");
        }
    } else {
        throw std::invalid_argument("Unsupported b-field type received");
    }
}

}  // namespace traccc::cuda
