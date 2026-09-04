/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#include <covfie/core/algebra/affine.hpp>
#include <covfie/core/backend/primitive/array.hpp>
#include <covfie/core/backend/transformer/affine.hpp>
#include <covfie/core/backend/transformer/nearest_neighbour.hpp>
#include <covfie/core/backend/transformer/strided.hpp>
#include <covfie/core/field.hpp>
#include <covfie/core/parameter_pack.hpp>

using backend_t = covfie::backend::affine<covfie::backend::nearest_neighbour<
    covfie::backend::strided<covfie::vector::size3,
                             covfie::backend::array<covfie::vector::float3>>>>;

using field_t = covfie::field<backend_t>;
