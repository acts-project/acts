/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Covfie include(s).
#include <covfie/core/backend/primitive/constant.hpp>
#include <covfie/core/backend/transformer/affine.hpp>
#include <covfie/core/backend/transformer/clamp.hpp>
#include <covfie/core/backend/transformer/linear.hpp>
#include <covfie/core/backend/transformer/strided.hpp>
#include <covfie/core/concepts.hpp>
#include <covfie/core/field.hpp>
#include <covfie/core/vector.hpp>
#include <covfie/sycl/backend/primitive/sycl_device_array.hpp>

#include "traccc/bfield/magnetic_field.hpp"
#include "traccc/bfield/magnetic_field_types.hpp"
#include "traccc/definitions/primitives.hpp"

namespace traccc::sycl {

/// Inhomogeneous B-field backend type for SYCL
template <typename scalar_t>
using inhom_bfield_backend_t =
    covfie::backend::affine<covfie::backend::linear<covfie::backend::clamp<
        covfie::backend::strided<covfie::vector::vector_d<std::size_t, 3>,
                                 covfie::backend::sycl_device_array<
                                     covfie::vector::vector_d<scalar_t, 3>>>>>>;
// Test that the type is a valid backend for a field
static_assert(covfie::concepts::field_backend<inhom_bfield_backend_t<float>>,
              "sycl::inhom_bfield_backend_t is not a valid field backend type");

/// @brief the standard list of SYCL bfield types to support
template <typename scalar_t>
using bfield_type_list = std::tuple<const_bfield_backend_t<scalar_t>,
                                    inhom_bfield_backend_t<scalar_t>>;

/*
 * SYCL requires a little bit of extra massaging to make the kernel tags
 * work...
 */
struct const_bfield_kernel_tag {};
struct inhom_bfield_kernel_tag {};

template <typename T>
struct bfield_tag_selector {};

template <typename scalar_t>
struct bfield_tag_selector<const_bfield_backend_t<scalar_t>> {
    using type = const_bfield_kernel_tag;
};

template <typename scalar_t>
struct bfield_tag_selector<inhom_bfield_backend_t<scalar_t>> {
    using type = inhom_bfield_kernel_tag;
};

template <typename T>
using bfield_tag_selector_t = typename bfield_tag_selector<T>::type;

template <typename T>
concept bfield_tag_exists_for_backend =
    requires { typename bfield_tag_selector_t<T>; };

template <typename>
struct bfield_tag_existance_validator {};

template <typename... Ts>
struct bfield_tag_existance_validator<std::tuple<Ts...>> {
    static constexpr bool value = (bfield_tag_exists_for_backend<Ts> && ...);
};

static_assert(
    bfield_tag_existance_validator<bfield_type_list<scalar>>::value,
    "Not all B-field types registered for SYCL have an accompanying tag");

}  // namespace traccc::sycl
