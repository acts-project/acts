// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

// Project include(s)
#include "detray/algebra/common/ksm/concepts.hpp"
#include "detray/algebra/common/ksm/matrix.hpp"

// System include(s)
#include <cstddef>

// This file defines some of detray's traits for the known-substructure
// matrices, so that the standard traits can be used with these types.
namespace detray::traits {

// The index type of any known-substructure matrix is trivially std::size_t,
// at least for now.
template <typename substructure_t, typename scalar_t>
struct index<::detray::ksm::matrix<substructure_t, scalar_t>> {
  using type = std::size_t;
};

// Similarly, a known-substructure matrix always has 2 dimensions, and we can
// fetch the number of rows and columns from the substructure itself.
template <typename substructure_t, typename scalar_t>
struct dimensions<::detray::ksm::matrix<substructure_t, scalar_t>> {
  using size_type = std::size_t;

  static constexpr size_type _dim{2};
  static constexpr size_type _rows{substructure_t::rows};
  static constexpr size_type _columns{substructure_t::columns};
};

// Each known-substructure matrix also embeds a scalar type, which we can
// use for the `scalar` trait.
template <typename substructure_t, typename scalar_t>
struct scalar<::detray::ksm::matrix<substructure_t, scalar_t>> {
  using type = scalar_t;
};

// Symmetry and diagonality of a known-substructure matrix. This allows us to
// make powerful assertions about the shape of a matrix at compile time. We
// define these traits in two parts: whether we can statically assert a
// property at all (e.g. a normal dense matrix has no way to assert symmetry,
// that's a runtime property) and another trait to actually assert it.
//
// First, we define checks for symmetry using KSM's existing concepts.
template <typename substructure_t, typename scalar_t>
struct static_symmetric_assertable<
    ::detray::ksm::matrix<substructure_t, scalar_t>> {
  static constexpr bool value = true;
};

template <typename substructure_t, typename scalar_t>
struct static_symmetric<::detray::ksm::matrix<substructure_t, scalar_t>> {
  static constexpr bool value =
      ::detray::ksm::concepts::is_symmetric<typename ::detray::ksm::matrix<
          substructure_t, scalar_t>::canonical_substructure_type>;
};

// And then we define checks for diagonality. Again, nothing more involved
// than reusing the existing KSM concepts.
template <typename substructure_t, typename scalar_t>
struct static_diagonal_assertable<
    ::detray::ksm::matrix<substructure_t, scalar_t>> {
  static constexpr bool value = true;
};

template <typename substructure_t, typename scalar_t>
struct static_diagonal<::detray::ksm::matrix<substructure_t, scalar_t>> {
  static constexpr bool value =
      ::detray::ksm::concepts::is_diagonal<typename ::detray::ksm::matrix<
          substructure_t, scalar_t>::canonical_substructure_type>;
};

}  // namespace detray::traits
