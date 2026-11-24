// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#if DETRAY_ALGEBRA_ARRAY
#include "algebra/array.hpp"
#elif DETRAY_ALGEBRA_EIGEN
#include "algebra/eigen.hpp"
#elif DETRAY_ALGEBRA_FASTOR
#include "algebra/fastor.hpp"
#elif DETRAY_ALGEBRA_SMATRIX
#include "algebra/smatrix.hpp"
#elif DETRAY_ALGEBRA_VC_AOS || DETRAY_ALGEBRA_VC_SOA

#if DETRAY_ALGEBRA_VC_AOS
#include "algebra/vc_aos.hpp"
#endif

#if DETRAY_ALGEBRA_VC_SOA
#include "algebra/vc_soa.hpp"
#endif

#else
#error "No algebra plugin selected! Please link to one of the algebra plugins."
#endif

// Algebra-plugins include(s)
#include "detray/algebra/common/boolean.hpp"
#include "detray/algebra/common/math.hpp"
#include "detray/algebra/utils/approximately_equal.hpp"
#include "detray/algebra/utils/casts.hpp"
#include "detray/algebra/utils/print.hpp"

namespace detray {

template <typename A, typename T>
/// Select the algebra backend with value type T.
using dalgebra = get_algebra_t<A, T>;

template <typename A>
/// Value scalar type of the algebra backend.
using dvalue = get_value_t<A>;

template <typename A>
/// Boolean type associated with the algebra backend.
using dbool = get_boolean_t<A>;

template <typename A, typename T>
/// SIMD type associated with value type T.
using dsimd = get_simd_t<A, T>;

template <typename A>
/// Index type provided by the algebra backend.
using dindex_type = get_index_t<A>;

template <typename A, typename T = dvalue<A>>
/// Scalar type for algebra A and value type T.
using dscalar = get_scalar_t<A, T>;

template <typename A, typename T = dvalue<A>>
/// Two-dimensional point type for algebra A and value type T.
using dpoint2D = get_point2D_t<A, T>;

template <typename A, typename T = dvalue<A>>
/// Three-dimensional point type for algebra A and value type T.
using dpoint3D = get_point3D_t<A, T>;

template <typename A, typename T = dvalue<A>>
/// Two-dimensional vector type for algebra A and value type T.
using dvector2D = get_vector2D_t<A, T>;

template <typename A, typename T = dvalue<A>>
/// Three-dimensional vector type for algebra A and value type T.
using dvector3D = get_vector3D_t<A, T>;

template <typename A, typename T = dvalue<A>>
/// Three-dimensional transform type for algebra A and value type T.
using dtransform3D = get_transform3D_t<A, T>;

template <typename A, std::size_t R, std::size_t C, typename T = dvalue<A>>
/// R by C matrix type for algebra A and value type T.
using dmatrix = get_matrix_t<A, R, C, T>;

}  // namespace detray
