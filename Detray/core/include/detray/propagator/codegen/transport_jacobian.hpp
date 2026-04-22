// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "detray/algebra/type_traits.hpp"
#include "detray/definitions/algebra.hpp"
#include "detray/definitions/detail/qualifiers.hpp"

namespace detray::detail {
template <concepts::algebra algebra_t>
struct transport_jacobian_matrix_without_gradient {
  using algebra_type = algebra_t;
  using scalar_type = dscalar<algebra_type>;
  using value_type = scalar_type;
  template <std::size_t I, std::size_t J>
  DETRAY_HOST_DEVICE scalar_type element() const
    requires(I < 8 && J < 8)
  {
    if constexpr (I == 0 && J == 0) {
      return static_cast<scalar_type>(1.f);
    } else if constexpr (I == 0 && J == 1) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 0 && J == 2) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 0 && J == 3) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 0 && J == 4) {
      return m_contents[0];
    } else if constexpr (I == 0 && J == 5) {
      return m_contents[1];
    } else if constexpr (I == 0 && J == 6) {
      return m_contents[2];
    } else if constexpr (I == 0 && J == 7) {
      return m_contents[3];
    } else if constexpr (I == 1 && J == 0) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 1 && J == 1) {
      return static_cast<scalar_type>(1.f);
    } else if constexpr (I == 1 && J == 2) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 1 && J == 3) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 1 && J == 4) {
      return m_contents[4];
    } else if constexpr (I == 1 && J == 5) {
      return m_contents[5];
    } else if constexpr (I == 1 && J == 6) {
      return m_contents[6];
    } else if constexpr (I == 1 && J == 7) {
      return m_contents[7];
    } else if constexpr (I == 2 && J == 0) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 2 && J == 1) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 2 && J == 2) {
      return static_cast<scalar_type>(1.f);
    } else if constexpr (I == 2 && J == 3) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 2 && J == 4) {
      return m_contents[8];
    } else if constexpr (I == 2 && J == 5) {
      return m_contents[9];
    } else if constexpr (I == 2 && J == 6) {
      return m_contents[10];
    } else if constexpr (I == 2 && J == 7) {
      return m_contents[11];
    } else if constexpr (I == 3 && J == 0) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 3 && J == 1) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 3 && J == 2) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 3 && J == 3) {
      return static_cast<scalar_type>(1.f);
    } else if constexpr (I == 3 && J == 4) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 3 && J == 5) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 3 && J == 6) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 3 && J == 7) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 4 && J == 0) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 4 && J == 1) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 4 && J == 2) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 4 && J == 3) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 4 && J == 4) {
      return m_contents[12];
    } else if constexpr (I == 4 && J == 5) {
      return m_contents[13];
    } else if constexpr (I == 4 && J == 6) {
      return m_contents[14];
    } else if constexpr (I == 4 && J == 7) {
      return m_contents[15];
    } else if constexpr (I == 5 && J == 0) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 5 && J == 1) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 5 && J == 2) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 5 && J == 3) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 5 && J == 4) {
      return m_contents[16];
    } else if constexpr (I == 5 && J == 5) {
      return m_contents[17];
    } else if constexpr (I == 5 && J == 6) {
      return m_contents[18];
    } else if constexpr (I == 5 && J == 7) {
      return m_contents[19];
    } else if constexpr (I == 6 && J == 0) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 6 && J == 1) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 6 && J == 2) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 6 && J == 3) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 6 && J == 4) {
      return m_contents[20];
    } else if constexpr (I == 6 && J == 5) {
      return m_contents[21];
    } else if constexpr (I == 6 && J == 6) {
      return m_contents[22];
    } else if constexpr (I == 6 && J == 7) {
      return m_contents[23];
    } else if constexpr (I == 7 && J == 0) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 7 && J == 1) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 7 && J == 2) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 7 && J == 3) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 7 && J == 4) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 7 && J == 5) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 7 && J == 6) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 7 && J == 7) {
      return m_contents[24];
    }
  }
  template <std::size_t I, std::size_t J>
  DETRAY_HOST_DEVICE scalar_type& element()
    requires((I == 0 && J == 4) || (I == 0 && J == 5) || (I == 0 && J == 6) ||
             (I == 0 && J == 7) || (I == 1 && J == 4) || (I == 1 && J == 5) ||
             (I == 1 && J == 6) || (I == 1 && J == 7) || (I == 2 && J == 4) ||
             (I == 2 && J == 5) || (I == 2 && J == 6) || (I == 2 && J == 7) ||
             (I == 4 && J == 4) || (I == 4 && J == 5) || (I == 4 && J == 6) ||
             (I == 4 && J == 7) || (I == 5 && J == 4) || (I == 5 && J == 5) ||
             (I == 5 && J == 6) || (I == 5 && J == 7) || (I == 6 && J == 4) ||
             (I == 6 && J == 5) || (I == 6 && J == 6) || (I == 6 && J == 7) ||
             (I == 7 && J == 7))
  {
    if constexpr (I == 0 && J == 4) {
      return m_contents[0];
    } else if constexpr (I == 0 && J == 5) {
      return m_contents[1];
    } else if constexpr (I == 0 && J == 6) {
      return m_contents[2];
    } else if constexpr (I == 0 && J == 7) {
      return m_contents[3];
    } else if constexpr (I == 1 && J == 4) {
      return m_contents[4];
    } else if constexpr (I == 1 && J == 5) {
      return m_contents[5];
    } else if constexpr (I == 1 && J == 6) {
      return m_contents[6];
    } else if constexpr (I == 1 && J == 7) {
      return m_contents[7];
    } else if constexpr (I == 2 && J == 4) {
      return m_contents[8];
    } else if constexpr (I == 2 && J == 5) {
      return m_contents[9];
    } else if constexpr (I == 2 && J == 6) {
      return m_contents[10];
    } else if constexpr (I == 2 && J == 7) {
      return m_contents[11];
    } else if constexpr (I == 4 && J == 4) {
      return m_contents[12];
    } else if constexpr (I == 4 && J == 5) {
      return m_contents[13];
    } else if constexpr (I == 4 && J == 6) {
      return m_contents[14];
    } else if constexpr (I == 4 && J == 7) {
      return m_contents[15];
    } else if constexpr (I == 5 && J == 4) {
      return m_contents[16];
    } else if constexpr (I == 5 && J == 5) {
      return m_contents[17];
    } else if constexpr (I == 5 && J == 6) {
      return m_contents[18];
    } else if constexpr (I == 5 && J == 7) {
      return m_contents[19];
    } else if constexpr (I == 6 && J == 4) {
      return m_contents[20];
    } else if constexpr (I == 6 && J == 5) {
      return m_contents[21];
    } else if constexpr (I == 6 && J == 6) {
      return m_contents[22];
    } else if constexpr (I == 6 && J == 7) {
      return m_contents[23];
    } else if constexpr (I == 7 && J == 7) {
      return m_contents[24];
    }
  }
  DETRAY_HOST_DEVICE explicit operator dmatrix<algebra_t, 8, 8>() const {
    dmatrix<algebra_t, 8, 8> rv;
    getter::element<0, 0>(rv) = element<0, 0>();
    getter::element<0, 1>(rv) = element<0, 1>();
    getter::element<0, 2>(rv) = element<0, 2>();
    getter::element<0, 3>(rv) = element<0, 3>();
    getter::element<0, 4>(rv) = element<0, 4>();
    getter::element<0, 5>(rv) = element<0, 5>();
    getter::element<0, 6>(rv) = element<0, 6>();
    getter::element<0, 7>(rv) = element<0, 7>();
    getter::element<1, 0>(rv) = element<1, 0>();
    getter::element<1, 1>(rv) = element<1, 1>();
    getter::element<1, 2>(rv) = element<1, 2>();
    getter::element<1, 3>(rv) = element<1, 3>();
    getter::element<1, 4>(rv) = element<1, 4>();
    getter::element<1, 5>(rv) = element<1, 5>();
    getter::element<1, 6>(rv) = element<1, 6>();
    getter::element<1, 7>(rv) = element<1, 7>();
    getter::element<2, 0>(rv) = element<2, 0>();
    getter::element<2, 1>(rv) = element<2, 1>();
    getter::element<2, 2>(rv) = element<2, 2>();
    getter::element<2, 3>(rv) = element<2, 3>();
    getter::element<2, 4>(rv) = element<2, 4>();
    getter::element<2, 5>(rv) = element<2, 5>();
    getter::element<2, 6>(rv) = element<2, 6>();
    getter::element<2, 7>(rv) = element<2, 7>();
    getter::element<3, 0>(rv) = element<3, 0>();
    getter::element<3, 1>(rv) = element<3, 1>();
    getter::element<3, 2>(rv) = element<3, 2>();
    getter::element<3, 3>(rv) = element<3, 3>();
    getter::element<3, 4>(rv) = element<3, 4>();
    getter::element<3, 5>(rv) = element<3, 5>();
    getter::element<3, 6>(rv) = element<3, 6>();
    getter::element<3, 7>(rv) = element<3, 7>();
    getter::element<4, 0>(rv) = element<4, 0>();
    getter::element<4, 1>(rv) = element<4, 1>();
    getter::element<4, 2>(rv) = element<4, 2>();
    getter::element<4, 3>(rv) = element<4, 3>();
    getter::element<4, 4>(rv) = element<4, 4>();
    getter::element<4, 5>(rv) = element<4, 5>();
    getter::element<4, 6>(rv) = element<4, 6>();
    getter::element<4, 7>(rv) = element<4, 7>();
    getter::element<5, 0>(rv) = element<5, 0>();
    getter::element<5, 1>(rv) = element<5, 1>();
    getter::element<5, 2>(rv) = element<5, 2>();
    getter::element<5, 3>(rv) = element<5, 3>();
    getter::element<5, 4>(rv) = element<5, 4>();
    getter::element<5, 5>(rv) = element<5, 5>();
    getter::element<5, 6>(rv) = element<5, 6>();
    getter::element<5, 7>(rv) = element<5, 7>();
    getter::element<6, 0>(rv) = element<6, 0>();
    getter::element<6, 1>(rv) = element<6, 1>();
    getter::element<6, 2>(rv) = element<6, 2>();
    getter::element<6, 3>(rv) = element<6, 3>();
    getter::element<6, 4>(rv) = element<6, 4>();
    getter::element<6, 5>(rv) = element<6, 5>();
    getter::element<6, 6>(rv) = element<6, 6>();
    getter::element<6, 7>(rv) = element<6, 7>();
    getter::element<7, 0>(rv) = element<7, 0>();
    getter::element<7, 1>(rv) = element<7, 1>();
    getter::element<7, 2>(rv) = element<7, 2>();
    getter::element<7, 3>(rv) = element<7, 3>();
    getter::element<7, 4>(rv) = element<7, 4>();
    getter::element<7, 5>(rv) = element<7, 5>();
    getter::element<7, 6>(rv) = element<7, 6>();
    getter::element<7, 7>(rv) = element<7, 7>();
    return rv;
  }
  DETRAY_HOST_DEVICE static constexpr transport_jacobian_matrix_without_gradient
  identity() {
    transport_jacobian_matrix_without_gradient rv;
    rv.element<0, 4>() = static_cast<scalar_type>(0.f);
    rv.element<0, 5>() = static_cast<scalar_type>(0.f);
    rv.element<0, 6>() = static_cast<scalar_type>(0.f);
    rv.element<0, 7>() = static_cast<scalar_type>(0.f);
    rv.element<1, 4>() = static_cast<scalar_type>(0.f);
    rv.element<1, 5>() = static_cast<scalar_type>(0.f);
    rv.element<1, 6>() = static_cast<scalar_type>(0.f);
    rv.element<1, 7>() = static_cast<scalar_type>(0.f);
    rv.element<2, 4>() = static_cast<scalar_type>(0.f);
    rv.element<2, 5>() = static_cast<scalar_type>(0.f);
    rv.element<2, 6>() = static_cast<scalar_type>(0.f);
    rv.element<2, 7>() = static_cast<scalar_type>(0.f);
    rv.element<4, 4>() = static_cast<scalar_type>(1.f);
    rv.element<4, 5>() = static_cast<scalar_type>(0.f);
    rv.element<4, 6>() = static_cast<scalar_type>(0.f);
    rv.element<4, 7>() = static_cast<scalar_type>(0.f);
    rv.element<5, 4>() = static_cast<scalar_type>(0.f);
    rv.element<5, 5>() = static_cast<scalar_type>(1.f);
    rv.element<5, 6>() = static_cast<scalar_type>(0.f);
    rv.element<5, 7>() = static_cast<scalar_type>(0.f);
    rv.element<6, 4>() = static_cast<scalar_type>(0.f);
    rv.element<6, 5>() = static_cast<scalar_type>(0.f);
    rv.element<6, 6>() = static_cast<scalar_type>(1.f);
    rv.element<6, 7>() = static_cast<scalar_type>(0.f);
    rv.element<7, 7>() = static_cast<scalar_type>(1.f);
    return rv;
  }

 private:
  std::array<scalar_type, 25> m_contents;
};
template <concepts::algebra algebra_t>
struct transport_jacobian_matrix_with_gradient {
  using algebra_type = algebra_t;
  using scalar_type = dscalar<algebra_type>;
  using value_type = scalar_type;
  template <std::size_t I, std::size_t J>
  DETRAY_HOST_DEVICE scalar_type element() const
    requires(I < 8 && J < 8)
  {
    if constexpr (I == 0 && J == 0) {
      return m_contents[0];
    } else if constexpr (I == 0 && J == 1) {
      return m_contents[1];
    } else if constexpr (I == 0 && J == 2) {
      return m_contents[2];
    } else if constexpr (I == 0 && J == 3) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 0 && J == 4) {
      return m_contents[3];
    } else if constexpr (I == 0 && J == 5) {
      return m_contents[4];
    } else if constexpr (I == 0 && J == 6) {
      return m_contents[5];
    } else if constexpr (I == 0 && J == 7) {
      return m_contents[6];
    } else if constexpr (I == 1 && J == 0) {
      return m_contents[7];
    } else if constexpr (I == 1 && J == 1) {
      return m_contents[8];
    } else if constexpr (I == 1 && J == 2) {
      return m_contents[9];
    } else if constexpr (I == 1 && J == 3) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 1 && J == 4) {
      return m_contents[10];
    } else if constexpr (I == 1 && J == 5) {
      return m_contents[11];
    } else if constexpr (I == 1 && J == 6) {
      return m_contents[12];
    } else if constexpr (I == 1 && J == 7) {
      return m_contents[13];
    } else if constexpr (I == 2 && J == 0) {
      return m_contents[14];
    } else if constexpr (I == 2 && J == 1) {
      return m_contents[15];
    } else if constexpr (I == 2 && J == 2) {
      return m_contents[16];
    } else if constexpr (I == 2 && J == 3) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 2 && J == 4) {
      return m_contents[17];
    } else if constexpr (I == 2 && J == 5) {
      return m_contents[18];
    } else if constexpr (I == 2 && J == 6) {
      return m_contents[19];
    } else if constexpr (I == 2 && J == 7) {
      return m_contents[20];
    } else if constexpr (I == 3 && J == 0) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 3 && J == 1) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 3 && J == 2) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 3 && J == 3) {
      return static_cast<scalar_type>(1.f);
    } else if constexpr (I == 3 && J == 4) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 3 && J == 5) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 3 && J == 6) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 3 && J == 7) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 4 && J == 0) {
      return m_contents[21];
    } else if constexpr (I == 4 && J == 1) {
      return m_contents[22];
    } else if constexpr (I == 4 && J == 2) {
      return m_contents[23];
    } else if constexpr (I == 4 && J == 3) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 4 && J == 4) {
      return m_contents[24];
    } else if constexpr (I == 4 && J == 5) {
      return m_contents[25];
    } else if constexpr (I == 4 && J == 6) {
      return m_contents[26];
    } else if constexpr (I == 4 && J == 7) {
      return m_contents[27];
    } else if constexpr (I == 5 && J == 0) {
      return m_contents[28];
    } else if constexpr (I == 5 && J == 1) {
      return m_contents[29];
    } else if constexpr (I == 5 && J == 2) {
      return m_contents[30];
    } else if constexpr (I == 5 && J == 3) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 5 && J == 4) {
      return m_contents[31];
    } else if constexpr (I == 5 && J == 5) {
      return m_contents[32];
    } else if constexpr (I == 5 && J == 6) {
      return m_contents[33];
    } else if constexpr (I == 5 && J == 7) {
      return m_contents[34];
    } else if constexpr (I == 6 && J == 0) {
      return m_contents[35];
    } else if constexpr (I == 6 && J == 1) {
      return m_contents[36];
    } else if constexpr (I == 6 && J == 2) {
      return m_contents[37];
    } else if constexpr (I == 6 && J == 3) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 6 && J == 4) {
      return m_contents[38];
    } else if constexpr (I == 6 && J == 5) {
      return m_contents[39];
    } else if constexpr (I == 6 && J == 6) {
      return m_contents[40];
    } else if constexpr (I == 6 && J == 7) {
      return m_contents[41];
    } else if constexpr (I == 7 && J == 0) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 7 && J == 1) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 7 && J == 2) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 7 && J == 3) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 7 && J == 4) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 7 && J == 5) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 7 && J == 6) {
      return static_cast<scalar_type>(0.f);
    } else if constexpr (I == 7 && J == 7) {
      return m_contents[42];
    }
  }
  template <std::size_t I, std::size_t J>
  DETRAY_HOST_DEVICE scalar_type& element()
    requires((I == 0 && J == 0) || (I == 0 && J == 1) || (I == 0 && J == 2) ||
             (I == 0 && J == 4) || (I == 0 && J == 5) || (I == 0 && J == 6) ||
             (I == 0 && J == 7) || (I == 1 && J == 0) || (I == 1 && J == 1) ||
             (I == 1 && J == 2) || (I == 1 && J == 4) || (I == 1 && J == 5) ||
             (I == 1 && J == 6) || (I == 1 && J == 7) || (I == 2 && J == 0) ||
             (I == 2 && J == 1) || (I == 2 && J == 2) || (I == 2 && J == 4) ||
             (I == 2 && J == 5) || (I == 2 && J == 6) || (I == 2 && J == 7) ||
             (I == 4 && J == 0) || (I == 4 && J == 1) || (I == 4 && J == 2) ||
             (I == 4 && J == 4) || (I == 4 && J == 5) || (I == 4 && J == 6) ||
             (I == 4 && J == 7) || (I == 5 && J == 0) || (I == 5 && J == 1) ||
             (I == 5 && J == 2) || (I == 5 && J == 4) || (I == 5 && J == 5) ||
             (I == 5 && J == 6) || (I == 5 && J == 7) || (I == 6 && J == 0) ||
             (I == 6 && J == 1) || (I == 6 && J == 2) || (I == 6 && J == 4) ||
             (I == 6 && J == 5) || (I == 6 && J == 6) || (I == 6 && J == 7) ||
             (I == 7 && J == 7))
  {
    if constexpr (I == 0 && J == 0) {
      return m_contents[0];
    } else if constexpr (I == 0 && J == 1) {
      return m_contents[1];
    } else if constexpr (I == 0 && J == 2) {
      return m_contents[2];
    } else if constexpr (I == 0 && J == 4) {
      return m_contents[3];
    } else if constexpr (I == 0 && J == 5) {
      return m_contents[4];
    } else if constexpr (I == 0 && J == 6) {
      return m_contents[5];
    } else if constexpr (I == 0 && J == 7) {
      return m_contents[6];
    } else if constexpr (I == 1 && J == 0) {
      return m_contents[7];
    } else if constexpr (I == 1 && J == 1) {
      return m_contents[8];
    } else if constexpr (I == 1 && J == 2) {
      return m_contents[9];
    } else if constexpr (I == 1 && J == 4) {
      return m_contents[10];
    } else if constexpr (I == 1 && J == 5) {
      return m_contents[11];
    } else if constexpr (I == 1 && J == 6) {
      return m_contents[12];
    } else if constexpr (I == 1 && J == 7) {
      return m_contents[13];
    } else if constexpr (I == 2 && J == 0) {
      return m_contents[14];
    } else if constexpr (I == 2 && J == 1) {
      return m_contents[15];
    } else if constexpr (I == 2 && J == 2) {
      return m_contents[16];
    } else if constexpr (I == 2 && J == 4) {
      return m_contents[17];
    } else if constexpr (I == 2 && J == 5) {
      return m_contents[18];
    } else if constexpr (I == 2 && J == 6) {
      return m_contents[19];
    } else if constexpr (I == 2 && J == 7) {
      return m_contents[20];
    } else if constexpr (I == 4 && J == 0) {
      return m_contents[21];
    } else if constexpr (I == 4 && J == 1) {
      return m_contents[22];
    } else if constexpr (I == 4 && J == 2) {
      return m_contents[23];
    } else if constexpr (I == 4 && J == 4) {
      return m_contents[24];
    } else if constexpr (I == 4 && J == 5) {
      return m_contents[25];
    } else if constexpr (I == 4 && J == 6) {
      return m_contents[26];
    } else if constexpr (I == 4 && J == 7) {
      return m_contents[27];
    } else if constexpr (I == 5 && J == 0) {
      return m_contents[28];
    } else if constexpr (I == 5 && J == 1) {
      return m_contents[29];
    } else if constexpr (I == 5 && J == 2) {
      return m_contents[30];
    } else if constexpr (I == 5 && J == 4) {
      return m_contents[31];
    } else if constexpr (I == 5 && J == 5) {
      return m_contents[32];
    } else if constexpr (I == 5 && J == 6) {
      return m_contents[33];
    } else if constexpr (I == 5 && J == 7) {
      return m_contents[34];
    } else if constexpr (I == 6 && J == 0) {
      return m_contents[35];
    } else if constexpr (I == 6 && J == 1) {
      return m_contents[36];
    } else if constexpr (I == 6 && J == 2) {
      return m_contents[37];
    } else if constexpr (I == 6 && J == 4) {
      return m_contents[38];
    } else if constexpr (I == 6 && J == 5) {
      return m_contents[39];
    } else if constexpr (I == 6 && J == 6) {
      return m_contents[40];
    } else if constexpr (I == 6 && J == 7) {
      return m_contents[41];
    } else if constexpr (I == 7 && J == 7) {
      return m_contents[42];
    }
  }
  DETRAY_HOST_DEVICE explicit operator dmatrix<algebra_t, 8, 8>() const {
    dmatrix<algebra_t, 8, 8> rv;
    getter::element<0, 0>(rv) = element<0, 0>();
    getter::element<0, 1>(rv) = element<0, 1>();
    getter::element<0, 2>(rv) = element<0, 2>();
    getter::element<0, 3>(rv) = element<0, 3>();
    getter::element<0, 4>(rv) = element<0, 4>();
    getter::element<0, 5>(rv) = element<0, 5>();
    getter::element<0, 6>(rv) = element<0, 6>();
    getter::element<0, 7>(rv) = element<0, 7>();
    getter::element<1, 0>(rv) = element<1, 0>();
    getter::element<1, 1>(rv) = element<1, 1>();
    getter::element<1, 2>(rv) = element<1, 2>();
    getter::element<1, 3>(rv) = element<1, 3>();
    getter::element<1, 4>(rv) = element<1, 4>();
    getter::element<1, 5>(rv) = element<1, 5>();
    getter::element<1, 6>(rv) = element<1, 6>();
    getter::element<1, 7>(rv) = element<1, 7>();
    getter::element<2, 0>(rv) = element<2, 0>();
    getter::element<2, 1>(rv) = element<2, 1>();
    getter::element<2, 2>(rv) = element<2, 2>();
    getter::element<2, 3>(rv) = element<2, 3>();
    getter::element<2, 4>(rv) = element<2, 4>();
    getter::element<2, 5>(rv) = element<2, 5>();
    getter::element<2, 6>(rv) = element<2, 6>();
    getter::element<2, 7>(rv) = element<2, 7>();
    getter::element<3, 0>(rv) = element<3, 0>();
    getter::element<3, 1>(rv) = element<3, 1>();
    getter::element<3, 2>(rv) = element<3, 2>();
    getter::element<3, 3>(rv) = element<3, 3>();
    getter::element<3, 4>(rv) = element<3, 4>();
    getter::element<3, 5>(rv) = element<3, 5>();
    getter::element<3, 6>(rv) = element<3, 6>();
    getter::element<3, 7>(rv) = element<3, 7>();
    getter::element<4, 0>(rv) = element<4, 0>();
    getter::element<4, 1>(rv) = element<4, 1>();
    getter::element<4, 2>(rv) = element<4, 2>();
    getter::element<4, 3>(rv) = element<4, 3>();
    getter::element<4, 4>(rv) = element<4, 4>();
    getter::element<4, 5>(rv) = element<4, 5>();
    getter::element<4, 6>(rv) = element<4, 6>();
    getter::element<4, 7>(rv) = element<4, 7>();
    getter::element<5, 0>(rv) = element<5, 0>();
    getter::element<5, 1>(rv) = element<5, 1>();
    getter::element<5, 2>(rv) = element<5, 2>();
    getter::element<5, 3>(rv) = element<5, 3>();
    getter::element<5, 4>(rv) = element<5, 4>();
    getter::element<5, 5>(rv) = element<5, 5>();
    getter::element<5, 6>(rv) = element<5, 6>();
    getter::element<5, 7>(rv) = element<5, 7>();
    getter::element<6, 0>(rv) = element<6, 0>();
    getter::element<6, 1>(rv) = element<6, 1>();
    getter::element<6, 2>(rv) = element<6, 2>();
    getter::element<6, 3>(rv) = element<6, 3>();
    getter::element<6, 4>(rv) = element<6, 4>();
    getter::element<6, 5>(rv) = element<6, 5>();
    getter::element<6, 6>(rv) = element<6, 6>();
    getter::element<6, 7>(rv) = element<6, 7>();
    getter::element<7, 0>(rv) = element<7, 0>();
    getter::element<7, 1>(rv) = element<7, 1>();
    getter::element<7, 2>(rv) = element<7, 2>();
    getter::element<7, 3>(rv) = element<7, 3>();
    getter::element<7, 4>(rv) = element<7, 4>();
    getter::element<7, 5>(rv) = element<7, 5>();
    getter::element<7, 6>(rv) = element<7, 6>();
    getter::element<7, 7>(rv) = element<7, 7>();
    return rv;
  }
  DETRAY_HOST_DEVICE static constexpr transport_jacobian_matrix_with_gradient
  identity() {
    transport_jacobian_matrix_with_gradient rv;
    rv.element<0, 0>() = static_cast<scalar_type>(1.f);
    rv.element<0, 1>() = static_cast<scalar_type>(0.f);
    rv.element<0, 2>() = static_cast<scalar_type>(0.f);
    rv.element<0, 4>() = static_cast<scalar_type>(0.f);
    rv.element<0, 5>() = static_cast<scalar_type>(0.f);
    rv.element<0, 6>() = static_cast<scalar_type>(0.f);
    rv.element<0, 7>() = static_cast<scalar_type>(0.f);
    rv.element<1, 0>() = static_cast<scalar_type>(0.f);
    rv.element<1, 1>() = static_cast<scalar_type>(1.f);
    rv.element<1, 2>() = static_cast<scalar_type>(0.f);
    rv.element<1, 4>() = static_cast<scalar_type>(0.f);
    rv.element<1, 5>() = static_cast<scalar_type>(0.f);
    rv.element<1, 6>() = static_cast<scalar_type>(0.f);
    rv.element<1, 7>() = static_cast<scalar_type>(0.f);
    rv.element<2, 0>() = static_cast<scalar_type>(0.f);
    rv.element<2, 1>() = static_cast<scalar_type>(0.f);
    rv.element<2, 2>() = static_cast<scalar_type>(1.f);
    rv.element<2, 4>() = static_cast<scalar_type>(0.f);
    rv.element<2, 5>() = static_cast<scalar_type>(0.f);
    rv.element<2, 6>() = static_cast<scalar_type>(0.f);
    rv.element<2, 7>() = static_cast<scalar_type>(0.f);
    rv.element<4, 0>() = static_cast<scalar_type>(0.f);
    rv.element<4, 1>() = static_cast<scalar_type>(0.f);
    rv.element<4, 2>() = static_cast<scalar_type>(0.f);
    rv.element<4, 4>() = static_cast<scalar_type>(1.f);
    rv.element<4, 5>() = static_cast<scalar_type>(0.f);
    rv.element<4, 6>() = static_cast<scalar_type>(0.f);
    rv.element<4, 7>() = static_cast<scalar_type>(0.f);
    rv.element<5, 0>() = static_cast<scalar_type>(0.f);
    rv.element<5, 1>() = static_cast<scalar_type>(0.f);
    rv.element<5, 2>() = static_cast<scalar_type>(0.f);
    rv.element<5, 4>() = static_cast<scalar_type>(0.f);
    rv.element<5, 5>() = static_cast<scalar_type>(1.f);
    rv.element<5, 6>() = static_cast<scalar_type>(0.f);
    rv.element<5, 7>() = static_cast<scalar_type>(0.f);
    rv.element<6, 0>() = static_cast<scalar_type>(0.f);
    rv.element<6, 1>() = static_cast<scalar_type>(0.f);
    rv.element<6, 2>() = static_cast<scalar_type>(0.f);
    rv.element<6, 4>() = static_cast<scalar_type>(0.f);
    rv.element<6, 5>() = static_cast<scalar_type>(0.f);
    rv.element<6, 6>() = static_cast<scalar_type>(1.f);
    rv.element<6, 7>() = static_cast<scalar_type>(0.f);
    rv.element<7, 7>() = static_cast<scalar_type>(1.f);
    return rv;
  }

 private:
  std::array<scalar_type, 43> m_contents;
};
}  // namespace detray::detail
