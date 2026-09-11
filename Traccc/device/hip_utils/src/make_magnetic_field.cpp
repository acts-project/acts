/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025-2026 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

// Local include(s).
#include "traccc/hip/utils/make_magnetic_field.hpp"

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
#include <covfie/hip/backend/primitive/hip_device_array.hpp>

// System include(s).
#include <stdexcept>

namespace traccc::hip {

/// Inhomogeneous B-field backend type using HIP global memory
template <typename scalar_t>
using inhom_global_bfield_backend_t =
    covfie::backend::affine<covfie::backend::linear<covfie::backend::clamp<
        covfie::backend::strided<covfie::vector::vector_d<std::size_t, 3>,
                                 covfie::backend::hip_device_array<
                                     covfie::vector::vector_d<scalar_t, 3>>>>>>;
// Test that the type is a valid backend for a field
static_assert(
    covfie::concepts::field_backend<inhom_global_bfield_backend_t<float>>,
    "hip::inhom_global_bfield_backend_t is not a valid field backend type");

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
          covfie::field<hip::inhom_global_bfield_backend_t<float>>(in_field)};
    } else {
      throw std::invalid_argument(
          "Unsupported storage method chosen for inhomogeneous b-field");
    }
  } else if (bfield.is<host::inhom_bfield_backend_t<double>>()) {
    // At double precision we can only use global memory.
    if (storage == magnetic_field_storage::global_memory) {
      return magnetic_field{
          covfie::field<hip::inhom_global_bfield_backend_t<double>>(
              bfield.as_field<host::inhom_bfield_backend_t<double>>())};
    } else {
      throw std::invalid_argument(
          "Unsupported storage method chosen for inhomogeneous b-field");
    }
  } else {
    throw std::invalid_argument("Unsupported b-field type received");
  }
}

}  // namespace traccc::hip
