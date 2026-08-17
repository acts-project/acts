// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s).
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"
#include "detray/utils/known_substructure_matrix.hpp"
#include "detray/utils/known_substructure_matrix_types.hpp"

namespace detray::detail {

/// The full Jacobian in its known substructure: 25 of its 36 cells are free.
/// Every specialisation of @c update_full_jacobian produces this type, so the
/// frame the Jacobian was assembled for is not visible to the caller.
///
/// The substructure taken here is the one that holds whether or not the volume
/// carries material; the no-material variant is more structured still, and
/// worth selecting statically for a detector known to have none.
template <concepts::algebra algebra_t>
using full_jacobian_matrix =
    ksm::matrix<ksm::full_jacobian_substructure<true>, dscalar<algebra_t>>;

/// @brief The full Jacobian between two bound track states.
///
/// @c J_full = J_F2B * (D + I) * J_transport * J_B2F, where D is the outer
/// product of the path-to-free and free-to-path derivatives. Every factor is
/// sparse, so the product is carried out in the factors' substructures and
/// skips the terms that are structurally zero.
///
/// @tparam has_dpos_dangle whether the departure surface's frame gives
/// d(pos)/d(angle) a value. Only line frames do.
/// @tparam path_has_direction_terms whether the destination surface's frame
/// sets the direction terms of the path derivative. Again only line frames.
///
/// @note The bracketing is deliberate. Matrix multiplication is associative,
/// but which factors are contracted first decides how early the sparsity is
/// exploited, and the spread across the five bracketings of this chain is
/// large. Contracting the two outer pairs first -- @c (f2b * D_plus_I) and
/// @c (transport * b2f) -- is the cheapest: 612 flops with the field gradient
/// against 692 for the left-to-right order, and 362 against 414 once the
/// frame-dependent zeros are taken into account. All five bracketings derive
/// the same result substructure, so this is purely a matter of cost.
template <bool has_dpos_dangle, bool path_has_direction_terms,
          concepts::algebra algebra_t, typename transport_jacobian_t>
DETRAY_HOST_DEVICE auto update_full_jacobian(
    const transport_jacobian_t& transport_jacobian,
    const dmatrix<algebra_t, 3, 2>& b2f_dpos_dloc,
    const dmatrix<algebra_t, 3, 2>& b2f_ddir_dangle,
    [[maybe_unused]] const dmatrix<algebra_t, 3, 2>& b2f_dpos_dangle,
    const dmatrix<algebra_t, 8, 1>& path_to_free_derivative,
    const dmatrix<algebra_t, 1, 8>& free_to_path_derivative,
    const dmatrix<algebra_t, 2, 3>& f2b_dloc_dpos,
    const dmatrix<algebra_t, 2, 3>& f2b_dangle_ddir) {
  using scalar_t = dscalar<algebra_t>;

  // Bound-to-free. The d(pos)/d(angle) block is structurally zero for every
  // frame but the line frames, and is then not stored at all.
  ksm::matrix<ksm::bound_to_free_substructure<has_dpos_dangle>, scalar_t> b2f{};
  ksm::copy_elements_into<0u, 0u>(b2f, b2f_dpos_dloc);
  ksm::copy_elements_into<4u, 2u>(b2f, b2f_ddir_dangle);
  if constexpr (has_dpos_dangle) {
    ksm::copy_elements_into<0u, 2u>(b2f, b2f_dpos_dangle);
  }

  // Free-to-bound.
  ksm::matrix<ksm::free_to_bound_substructure, scalar_t> f2b{};
  ksm::copy_elements_into<0u, 0u>(f2b, f2b_dloc_dpos);
  ksm::copy_elements_into<2u, 4u>(f2b, f2b_dangle_ddir);

  // The path correction, D = path_to_free * free_to_path, plus the identity.
  ksm::matrix<ksm::path_to_free_substructure, scalar_t> p2f{};
  ksm::copy_elements_into<0u, 0u>(p2f, path_to_free_derivative);

  ksm::matrix<ksm::free_to_path_substructure<path_has_direction_terms>,
              scalar_t>
      f2p{};
  ksm::copy_elements_into<0u, 0u>(f2p, free_to_path_derivative);

  using identity_t =
      ksm::matrix<typename ksm::make_identity_substructure<8>::canonical_type,
                  scalar_t>;

  const auto d_plus_i = p2f * f2p + identity_t{};

  const auto result = (f2b * d_plus_i) * (transport_jacobian * b2f);

  // Every specialisation has to agree on the result type: the frame-dependent
  // zeros make the product cheaper, but they do not survive into the result,
  // so all four instantiations are interchangeable to the caller. That is what
  // lets get_full_jacobian pick one at run time and still have a single return
  // type -- and what lets the accumulated Jacobian keep the same substructure.
  static_assert(
      std::is_same_v<
          std::decay_t<decltype(result)>,
          ksm::matrix<ksm::full_jacobian_substructure<true>, scalar_t>>,
      "the assembled full Jacobian does not have the substructure declared in "
      "known_substructure_matrix_types.hpp");

  return result;
}

}  // namespace detray::detail
