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
template <typename transport_jacobian_t, typename b2f_dpos_dloc_t,
          typename b2f_ddir_dangle_t, typename b2f_dpos_dangle_t,
          typename path_to_free_derivative_t,
          typename free_to_path_derivative_t, typename f2b_dloc_dpos_t,
          typename f2b_dangle_ddir_t, typename full_jacobian_t>
DETRAY_HOST_DEVICE void inline update_full_jacobian_with_gradient_impl(
    const transport_jacobian_t& transport_jacobian,
    const b2f_dpos_dloc_t& b2f_dpos_dloc,
    const b2f_ddir_dangle_t& b2f_ddir_dangle,
    const b2f_dpos_dangle_t& b2f_dpos_dangle,
    const path_to_free_derivative_t& path_to_free_derivative,
    const free_to_path_derivative_t& free_to_path_derivative,
    const f2b_dloc_dpos_t& f2b_dloc_dpos,
    const f2b_dangle_ddir_t& f2b_dangle_ddir, full_jacobian_t& full_jacobian)
  requires((detray::concepts::square_matrix<transport_jacobian_t> &&
            detray::traits::max_rank<transport_jacobian_t> == 8) &&
           (detray::concepts::matrix<b2f_dpos_dloc_t> &&
            detray::traits::rows<b2f_dpos_dloc_t> == 3 &&
            detray::traits::columns<b2f_dpos_dloc_t> == 2) &&
           (detray::concepts::matrix<b2f_ddir_dangle_t> &&
            detray::traits::rows<b2f_ddir_dangle_t> == 3 &&
            detray::traits::columns<b2f_ddir_dangle_t> == 2) &&
           (detray::concepts::matrix<b2f_dpos_dangle_t> &&
            detray::traits::rows<b2f_dpos_dangle_t> == 3 &&
            detray::traits::columns<b2f_dpos_dangle_t> == 2) &&
           (detray::concepts::matrix<path_to_free_derivative_t> &&
            detray::traits::rows<path_to_free_derivative_t> == 8 &&
            detray::traits::columns<path_to_free_derivative_t> == 1) &&
           (detray::concepts::row_matrix<free_to_path_derivative_t> &&
            detray::traits::columns<free_to_path_derivative_t> == 8) &&
           (detray::concepts::matrix<f2b_dloc_dpos_t> &&
            detray::traits::rows<f2b_dloc_dpos_t> == 2 &&
            detray::traits::columns<f2b_dloc_dpos_t> == 3) &&
           (detray::concepts::matrix<f2b_dangle_ddir_t> &&
            detray::traits::rows<f2b_dangle_ddir_t> == 2 &&
            detray::traits::columns<f2b_dangle_ddir_t> == 3) &&
           (detray::concepts::square_matrix<full_jacobian_t> &&
            detray::traits::max_rank<full_jacobian_t> == 6))
{
  assert((getter::element<0, 3>(transport_jacobian) == 0.f));
  assert((getter::element<1, 3>(transport_jacobian) == 0.f));
  assert((getter::element<2, 3>(transport_jacobian) == 0.f));
  assert((getter::element<3, 0>(transport_jacobian) == 0.f));
  assert((getter::element<3, 1>(transport_jacobian) == 0.f));
  assert((getter::element<3, 2>(transport_jacobian) == 0.f));
  assert((getter::element<3, 3>(transport_jacobian) == 1.f));
  assert((getter::element<3, 4>(transport_jacobian) == 0.f));
  assert((getter::element<3, 5>(transport_jacobian) == 0.f));
  assert((getter::element<3, 6>(transport_jacobian) == 0.f));
  assert((getter::element<3, 7>(transport_jacobian) == 0.f));
  assert((getter::element<4, 3>(transport_jacobian) == 0.f));
  assert((getter::element<5, 3>(transport_jacobian) == 0.f));
  assert((getter::element<6, 3>(transport_jacobian) == 0.f));
  assert((getter::element<7, 0>(transport_jacobian) == 0.f));
  assert((getter::element<7, 1>(transport_jacobian) == 0.f));
  assert((getter::element<7, 2>(transport_jacobian) == 0.f));
  assert((getter::element<7, 3>(transport_jacobian) == 0.f));
  assert((getter::element<7, 4>(transport_jacobian) == 0.f));
  assert((getter::element<7, 5>(transport_jacobian) == 0.f));
  assert((getter::element<7, 6>(transport_jacobian) == 0.f));
  assert((getter::element<2, 0>(b2f_ddir_dangle) == 0.f));
  assert((getter::element<3, 0>(path_to_free_derivative) == 0.f));
  assert((getter::element<0, 3>(free_to_path_derivative) == 0.f));
  assert((getter::element<0, 7>(free_to_path_derivative) == 0.f));
  assert((getter::element<0, 2>(f2b_dangle_ddir) == 0.f));
  const auto x0 = getter::element<0, 0>(f2b_dloc_dpos) *
                  getter::element<0, 0>(path_to_free_derivative);
  const auto x1 = getter::element<0, 1>(f2b_dloc_dpos) *
                  getter::element<1, 0>(path_to_free_derivative);
  const auto x2 = getter::element<0, 2>(f2b_dloc_dpos) *
                  getter::element<2, 0>(path_to_free_derivative);
  const auto x6 = getter::element<0, 0>(free_to_path_derivative) *
                      getter::element<0, 0>(path_to_free_derivative) +
                  1;
  const auto x8 = getter::element<0, 1>(free_to_path_derivative) *
                      getter::element<1, 0>(path_to_free_derivative) +
                  1;
  const auto x10 = getter::element<0, 2>(free_to_path_derivative) *
                       getter::element<2, 0>(path_to_free_derivative) +
                   1;
  const auto x17 = getter::element<1, 0>(f2b_dloc_dpos) *
                   getter::element<0, 0>(path_to_free_derivative);
  const auto x18 = getter::element<1, 1>(f2b_dloc_dpos) *
                   getter::element<1, 0>(path_to_free_derivative);
  const auto x19 = getter::element<1, 2>(f2b_dloc_dpos) *
                   getter::element<2, 0>(path_to_free_derivative);
  const auto x31 = getter::element<0, 0>(f2b_dangle_ddir) *
                   getter::element<4, 0>(path_to_free_derivative);
  const auto x32 = getter::element<0, 1>(f2b_dangle_ddir) *
                   getter::element<5, 0>(path_to_free_derivative);
  const auto x37 = getter::element<0, 4>(free_to_path_derivative) *
                       getter::element<4, 0>(path_to_free_derivative) +
                   1;
  const auto x39 = getter::element<0, 5>(free_to_path_derivative) *
                       getter::element<5, 0>(path_to_free_derivative) +
                   1;
  const auto x46 = getter::element<1, 0>(f2b_dangle_ddir) *
                   getter::element<4, 0>(path_to_free_derivative);
  const auto x47 = getter::element<1, 1>(f2b_dangle_ddir) *
                   getter::element<5, 0>(path_to_free_derivative);
  const auto x48 = getter::element<1, 2>(f2b_dangle_ddir) *
                   getter::element<6, 0>(path_to_free_derivative);
  const auto x60 = getter::element<0, 0>(free_to_path_derivative) *
                   getter::element<7, 0>(path_to_free_derivative);
  const auto x61 = getter::element<0, 1>(free_to_path_derivative) *
                   getter::element<7, 0>(path_to_free_derivative);
  const auto x62 = getter::element<0, 2>(free_to_path_derivative) *
                   getter::element<7, 0>(path_to_free_derivative);
  const auto x63 = getter::element<0, 4>(free_to_path_derivative) *
                   getter::element<7, 0>(path_to_free_derivative);
  const auto x64 = getter::element<0, 5>(free_to_path_derivative) *
                   getter::element<7, 0>(path_to_free_derivative);
  const auto x65 = getter::element<0, 6>(free_to_path_derivative) *
                   getter::element<7, 0>(path_to_free_derivative);
  const auto x33 = (x31 + x32) * getter::element<0, 0>(free_to_path_derivative);
  const auto x34 = (x31 + x32) * getter::element<0, 1>(free_to_path_derivative);
  const auto x35 = (x31 + x32) * getter::element<0, 2>(free_to_path_derivative);
  const auto x36 = (x31 + x32) * getter::element<0, 6>(free_to_path_derivative);
  const auto x3 =
      (x0 + x1 + x2) * getter::element<0, 4>(free_to_path_derivative);
  const auto x4 =
      (x0 + x1 + x2) * getter::element<0, 5>(free_to_path_derivative);
  const auto x5 =
      (x0 + x1 + x2) * getter::element<0, 6>(free_to_path_derivative);
  const auto x20 =
      (x17 + x18 + x19) * getter::element<0, 4>(free_to_path_derivative);
  const auto x21 =
      (x17 + x18 + x19) * getter::element<0, 5>(free_to_path_derivative);
  const auto x22 =
      (x17 + x18 + x19) * getter::element<0, 6>(free_to_path_derivative);
  const auto x38 = x32 * getter::element<0, 4>(free_to_path_derivative) +
                   x37 * getter::element<0, 0>(f2b_dangle_ddir);
  const auto x40 = x31 * getter::element<0, 5>(free_to_path_derivative) +
                   x39 * getter::element<0, 1>(f2b_dangle_ddir);
  const auto x49 =
      (x46 + x47 + x48) * getter::element<0, 0>(free_to_path_derivative);
  const auto x50 =
      (x46 + x47 + x48) * getter::element<0, 1>(free_to_path_derivative);
  const auto x51 =
      (x46 + x47 + x48) * getter::element<0, 2>(free_to_path_derivative);
  const auto x7 = x1 * getter::element<0, 0>(free_to_path_derivative) +
                  x2 * getter::element<0, 0>(free_to_path_derivative) +
                  x6 * getter::element<0, 0>(f2b_dloc_dpos);
  const auto x9 = x0 * getter::element<0, 1>(free_to_path_derivative) +
                  x2 * getter::element<0, 1>(free_to_path_derivative) +
                  x8 * getter::element<0, 1>(f2b_dloc_dpos);
  const auto x11 = x0 * getter::element<0, 2>(free_to_path_derivative) +
                   x1 * getter::element<0, 2>(free_to_path_derivative) +
                   x10 * getter::element<0, 2>(f2b_dloc_dpos);
  const auto x23 = x18 * getter::element<0, 0>(free_to_path_derivative) +
                   x19 * getter::element<0, 0>(free_to_path_derivative) +
                   x6 * getter::element<1, 0>(f2b_dloc_dpos);
  const auto x24 = x17 * getter::element<0, 1>(free_to_path_derivative) +
                   x19 * getter::element<0, 1>(free_to_path_derivative) +
                   x8 * getter::element<1, 1>(f2b_dloc_dpos);
  const auto x25 = x10 * getter::element<1, 2>(f2b_dloc_dpos) +
                   x17 * getter::element<0, 2>(free_to_path_derivative) +
                   x18 * getter::element<0, 2>(free_to_path_derivative);
  const auto x52 = x37 * getter::element<1, 0>(f2b_dangle_ddir) +
                   x47 * getter::element<0, 4>(free_to_path_derivative) +
                   x48 * getter::element<0, 4>(free_to_path_derivative);
  const auto x53 = x39 * getter::element<1, 1>(f2b_dangle_ddir) +
                   x46 * getter::element<0, 5>(free_to_path_derivative) +
                   x48 * getter::element<0, 5>(free_to_path_derivative);
  const auto x54 = x46 * getter::element<0, 6>(free_to_path_derivative) +
                   x47 * getter::element<0, 6>(free_to_path_derivative) +
                   (getter::element<0, 6>(free_to_path_derivative) *
                        getter::element<6, 0>(path_to_free_derivative) +
                    1) *
                       getter::element<1, 2>(f2b_dangle_ddir);
  const auto x66 = x60 * getter::element<0, 0>(transport_jacobian) +
                   x61 * getter::element<1, 0>(transport_jacobian) +
                   x62 * getter::element<2, 0>(transport_jacobian) +
                   x63 * getter::element<4, 0>(transport_jacobian) +
                   x64 * getter::element<5, 0>(transport_jacobian) +
                   x65 * getter::element<6, 0>(transport_jacobian);
  const auto x67 = x60 * getter::element<0, 1>(transport_jacobian) +
                   x61 * getter::element<1, 1>(transport_jacobian) +
                   x62 * getter::element<2, 1>(transport_jacobian) +
                   x63 * getter::element<4, 1>(transport_jacobian) +
                   x64 * getter::element<5, 1>(transport_jacobian) +
                   x65 * getter::element<6, 1>(transport_jacobian);
  const auto x68 = x60 * getter::element<0, 2>(transport_jacobian) +
                   x61 * getter::element<1, 2>(transport_jacobian) +
                   x62 * getter::element<2, 2>(transport_jacobian) +
                   x63 * getter::element<4, 2>(transport_jacobian) +
                   x64 * getter::element<5, 2>(transport_jacobian) +
                   x65 * getter::element<6, 2>(transport_jacobian);
  const auto x69 = x60 * getter::element<0, 4>(transport_jacobian) +
                   x61 * getter::element<1, 4>(transport_jacobian) +
                   x62 * getter::element<2, 4>(transport_jacobian) +
                   x63 * getter::element<4, 4>(transport_jacobian) +
                   x64 * getter::element<5, 4>(transport_jacobian) +
                   x65 * getter::element<6, 4>(transport_jacobian);
  const auto x70 = x60 * getter::element<0, 5>(transport_jacobian) +
                   x61 * getter::element<1, 5>(transport_jacobian) +
                   x62 * getter::element<2, 5>(transport_jacobian) +
                   x63 * getter::element<4, 5>(transport_jacobian) +
                   x64 * getter::element<5, 5>(transport_jacobian) +
                   x65 * getter::element<6, 5>(transport_jacobian);
  const auto x12 = x11 * getter::element<2, 0>(transport_jacobian) +
                   x3 * getter::element<4, 0>(transport_jacobian) +
                   x4 * getter::element<5, 0>(transport_jacobian) +
                   x5 * getter::element<6, 0>(transport_jacobian) +
                   x7 * getter::element<0, 0>(transport_jacobian) +
                   x9 * getter::element<1, 0>(transport_jacobian);
  const auto x13 = x11 * getter::element<2, 1>(transport_jacobian) +
                   x3 * getter::element<4, 1>(transport_jacobian) +
                   x4 * getter::element<5, 1>(transport_jacobian) +
                   x5 * getter::element<6, 1>(transport_jacobian) +
                   x7 * getter::element<0, 1>(transport_jacobian) +
                   x9 * getter::element<1, 1>(transport_jacobian);
  const auto x14 = x11 * getter::element<2, 2>(transport_jacobian) +
                   x3 * getter::element<4, 2>(transport_jacobian) +
                   x4 * getter::element<5, 2>(transport_jacobian) +
                   x5 * getter::element<6, 2>(transport_jacobian) +
                   x7 * getter::element<0, 2>(transport_jacobian) +
                   x9 * getter::element<1, 2>(transport_jacobian);
  const auto x15 = x11 * getter::element<2, 4>(transport_jacobian) +
                   x3 * getter::element<4, 4>(transport_jacobian) +
                   x4 * getter::element<5, 4>(transport_jacobian) +
                   x5 * getter::element<6, 4>(transport_jacobian) +
                   x7 * getter::element<0, 4>(transport_jacobian) +
                   x9 * getter::element<1, 4>(transport_jacobian);
  const auto x16 = x11 * getter::element<2, 5>(transport_jacobian) +
                   x3 * getter::element<4, 5>(transport_jacobian) +
                   x4 * getter::element<5, 5>(transport_jacobian) +
                   x5 * getter::element<6, 5>(transport_jacobian) +
                   x7 * getter::element<0, 5>(transport_jacobian) +
                   x9 * getter::element<1, 5>(transport_jacobian);
  const auto x26 = x20 * getter::element<4, 0>(transport_jacobian) +
                   x21 * getter::element<5, 0>(transport_jacobian) +
                   x22 * getter::element<6, 0>(transport_jacobian) +
                   x23 * getter::element<0, 0>(transport_jacobian) +
                   x24 * getter::element<1, 0>(transport_jacobian) +
                   x25 * getter::element<2, 0>(transport_jacobian);
  const auto x27 = x20 * getter::element<4, 1>(transport_jacobian) +
                   x21 * getter::element<5, 1>(transport_jacobian) +
                   x22 * getter::element<6, 1>(transport_jacobian) +
                   x23 * getter::element<0, 1>(transport_jacobian) +
                   x24 * getter::element<1, 1>(transport_jacobian) +
                   x25 * getter::element<2, 1>(transport_jacobian);
  const auto x28 = x20 * getter::element<4, 2>(transport_jacobian) +
                   x21 * getter::element<5, 2>(transport_jacobian) +
                   x22 * getter::element<6, 2>(transport_jacobian) +
                   x23 * getter::element<0, 2>(transport_jacobian) +
                   x24 * getter::element<1, 2>(transport_jacobian) +
                   x25 * getter::element<2, 2>(transport_jacobian);
  const auto x29 = x20 * getter::element<4, 4>(transport_jacobian) +
                   x21 * getter::element<5, 4>(transport_jacobian) +
                   x22 * getter::element<6, 4>(transport_jacobian) +
                   x23 * getter::element<0, 4>(transport_jacobian) +
                   x24 * getter::element<1, 4>(transport_jacobian) +
                   x25 * getter::element<2, 4>(transport_jacobian);
  const auto x30 = x20 * getter::element<4, 5>(transport_jacobian) +
                   x21 * getter::element<5, 5>(transport_jacobian) +
                   x22 * getter::element<6, 5>(transport_jacobian) +
                   x23 * getter::element<0, 5>(transport_jacobian) +
                   x24 * getter::element<1, 5>(transport_jacobian) +
                   x25 * getter::element<2, 5>(transport_jacobian);
  const auto x41 = x33 * getter::element<0, 0>(transport_jacobian) +
                   x34 * getter::element<1, 0>(transport_jacobian) +
                   x35 * getter::element<2, 0>(transport_jacobian) +
                   x36 * getter::element<6, 0>(transport_jacobian) +
                   x38 * getter::element<4, 0>(transport_jacobian) +
                   x40 * getter::element<5, 0>(transport_jacobian);
  const auto x42 = x33 * getter::element<0, 1>(transport_jacobian) +
                   x34 * getter::element<1, 1>(transport_jacobian) +
                   x35 * getter::element<2, 1>(transport_jacobian) +
                   x36 * getter::element<6, 1>(transport_jacobian) +
                   x38 * getter::element<4, 1>(transport_jacobian) +
                   x40 * getter::element<5, 1>(transport_jacobian);
  const auto x43 = x33 * getter::element<0, 2>(transport_jacobian) +
                   x34 * getter::element<1, 2>(transport_jacobian) +
                   x35 * getter::element<2, 2>(transport_jacobian) +
                   x36 * getter::element<6, 2>(transport_jacobian) +
                   x38 * getter::element<4, 2>(transport_jacobian) +
                   x40 * getter::element<5, 2>(transport_jacobian);
  const auto x44 = x33 * getter::element<0, 4>(transport_jacobian) +
                   x34 * getter::element<1, 4>(transport_jacobian) +
                   x35 * getter::element<2, 4>(transport_jacobian) +
                   x36 * getter::element<6, 4>(transport_jacobian) +
                   x38 * getter::element<4, 4>(transport_jacobian) +
                   x40 * getter::element<5, 4>(transport_jacobian);
  const auto x45 = x33 * getter::element<0, 5>(transport_jacobian) +
                   x34 * getter::element<1, 5>(transport_jacobian) +
                   x35 * getter::element<2, 5>(transport_jacobian) +
                   x36 * getter::element<6, 5>(transport_jacobian) +
                   x38 * getter::element<4, 5>(transport_jacobian) +
                   x40 * getter::element<5, 5>(transport_jacobian);
  const auto x55 = x49 * getter::element<0, 0>(transport_jacobian) +
                   x50 * getter::element<1, 0>(transport_jacobian) +
                   x51 * getter::element<2, 0>(transport_jacobian) +
                   x52 * getter::element<4, 0>(transport_jacobian) +
                   x53 * getter::element<5, 0>(transport_jacobian) +
                   x54 * getter::element<6, 0>(transport_jacobian);
  const auto x56 = x49 * getter::element<0, 1>(transport_jacobian) +
                   x50 * getter::element<1, 1>(transport_jacobian) +
                   x51 * getter::element<2, 1>(transport_jacobian) +
                   x52 * getter::element<4, 1>(transport_jacobian) +
                   x53 * getter::element<5, 1>(transport_jacobian) +
                   x54 * getter::element<6, 1>(transport_jacobian);
  const auto x57 = x49 * getter::element<0, 2>(transport_jacobian) +
                   x50 * getter::element<1, 2>(transport_jacobian) +
                   x51 * getter::element<2, 2>(transport_jacobian) +
                   x52 * getter::element<4, 2>(transport_jacobian) +
                   x53 * getter::element<5, 2>(transport_jacobian) +
                   x54 * getter::element<6, 2>(transport_jacobian);
  const auto x58 = x49 * getter::element<0, 4>(transport_jacobian) +
                   x50 * getter::element<1, 4>(transport_jacobian) +
                   x51 * getter::element<2, 4>(transport_jacobian) +
                   x52 * getter::element<4, 4>(transport_jacobian) +
                   x53 * getter::element<5, 4>(transport_jacobian) +
                   x54 * getter::element<6, 4>(transport_jacobian);
  const auto x59 = x49 * getter::element<0, 5>(transport_jacobian) +
                   x50 * getter::element<1, 5>(transport_jacobian) +
                   x51 * getter::element<2, 5>(transport_jacobian) +
                   x52 * getter::element<4, 5>(transport_jacobian) +
                   x53 * getter::element<5, 5>(transport_jacobian) +
                   x54 * getter::element<6, 5>(transport_jacobian);
  getter::element<0, 0>(full_jacobian) =
      x12 * getter::element<0, 0>(b2f_dpos_dloc) +
      x13 * getter::element<1, 0>(b2f_dpos_dloc) +
      x14 * getter::element<2, 0>(b2f_dpos_dloc);
  getter::element<1, 0>(full_jacobian) =
      x26 * getter::element<0, 0>(b2f_dpos_dloc) +
      x27 * getter::element<1, 0>(b2f_dpos_dloc) +
      x28 * getter::element<2, 0>(b2f_dpos_dloc);
  getter::element<2, 0>(full_jacobian) =
      x41 * getter::element<0, 0>(b2f_dpos_dloc) +
      x42 * getter::element<1, 0>(b2f_dpos_dloc) +
      x43 * getter::element<2, 0>(b2f_dpos_dloc);
  getter::element<3, 0>(full_jacobian) =
      x55 * getter::element<0, 0>(b2f_dpos_dloc) +
      x56 * getter::element<1, 0>(b2f_dpos_dloc) +
      x57 * getter::element<2, 0>(b2f_dpos_dloc);
  getter::element<4, 0>(full_jacobian) =
      x66 * getter::element<0, 0>(b2f_dpos_dloc) +
      x67 * getter::element<1, 0>(b2f_dpos_dloc) +
      x68 * getter::element<2, 0>(b2f_dpos_dloc);
  getter::element<5, 0>(full_jacobian) = 0;
  getter::element<0, 1>(full_jacobian) =
      x12 * getter::element<0, 1>(b2f_dpos_dloc) +
      x13 * getter::element<1, 1>(b2f_dpos_dloc) +
      x14 * getter::element<2, 1>(b2f_dpos_dloc);
  getter::element<1, 1>(full_jacobian) =
      x26 * getter::element<0, 1>(b2f_dpos_dloc) +
      x27 * getter::element<1, 1>(b2f_dpos_dloc) +
      x28 * getter::element<2, 1>(b2f_dpos_dloc);
  getter::element<2, 1>(full_jacobian) =
      x41 * getter::element<0, 1>(b2f_dpos_dloc) +
      x42 * getter::element<1, 1>(b2f_dpos_dloc) +
      x43 * getter::element<2, 1>(b2f_dpos_dloc);
  getter::element<3, 1>(full_jacobian) =
      x55 * getter::element<0, 1>(b2f_dpos_dloc) +
      x56 * getter::element<1, 1>(b2f_dpos_dloc) +
      x57 * getter::element<2, 1>(b2f_dpos_dloc);
  getter::element<4, 1>(full_jacobian) =
      x66 * getter::element<0, 1>(b2f_dpos_dloc) +
      x67 * getter::element<1, 1>(b2f_dpos_dloc) +
      x68 * getter::element<2, 1>(b2f_dpos_dloc);
  getter::element<5, 1>(full_jacobian) = 0;
  getter::element<0, 2>(full_jacobian) =
      x12 * getter::element<0, 0>(b2f_dpos_dangle) +
      x13 * getter::element<1, 0>(b2f_dpos_dangle) +
      x14 * getter::element<2, 0>(b2f_dpos_dangle) +
      x15 * getter::element<0, 0>(b2f_ddir_dangle) +
      x16 * getter::element<1, 0>(b2f_ddir_dangle);
  getter::element<1, 2>(full_jacobian) =
      x26 * getter::element<0, 0>(b2f_dpos_dangle) +
      x27 * getter::element<1, 0>(b2f_dpos_dangle) +
      x28 * getter::element<2, 0>(b2f_dpos_dangle) +
      x29 * getter::element<0, 0>(b2f_ddir_dangle) +
      x30 * getter::element<1, 0>(b2f_ddir_dangle);
  getter::element<2, 2>(full_jacobian) =
      x41 * getter::element<0, 0>(b2f_dpos_dangle) +
      x42 * getter::element<1, 0>(b2f_dpos_dangle) +
      x43 * getter::element<2, 0>(b2f_dpos_dangle) +
      x44 * getter::element<0, 0>(b2f_ddir_dangle) +
      x45 * getter::element<1, 0>(b2f_ddir_dangle);
  getter::element<3, 2>(full_jacobian) =
      x55 * getter::element<0, 0>(b2f_dpos_dangle) +
      x56 * getter::element<1, 0>(b2f_dpos_dangle) +
      x57 * getter::element<2, 0>(b2f_dpos_dangle) +
      x58 * getter::element<0, 0>(b2f_ddir_dangle) +
      x59 * getter::element<1, 0>(b2f_ddir_dangle);
  getter::element<4, 2>(full_jacobian) =
      x66 * getter::element<0, 0>(b2f_dpos_dangle) +
      x67 * getter::element<1, 0>(b2f_dpos_dangle) +
      x68 * getter::element<2, 0>(b2f_dpos_dangle) +
      x69 * getter::element<0, 0>(b2f_ddir_dangle) +
      x70 * getter::element<1, 0>(b2f_ddir_dangle);
  getter::element<5, 2>(full_jacobian) = 0;
  getter::element<0, 3>(full_jacobian) =
      x12 * getter::element<0, 1>(b2f_dpos_dangle) +
      x13 * getter::element<1, 1>(b2f_dpos_dangle) +
      x14 * getter::element<2, 1>(b2f_dpos_dangle) +
      x15 * getter::element<0, 1>(b2f_ddir_dangle) +
      x16 * getter::element<1, 1>(b2f_ddir_dangle) +
      (x11 * getter::element<2, 6>(transport_jacobian) +
       x3 * getter::element<4, 6>(transport_jacobian) +
       x4 * getter::element<5, 6>(transport_jacobian) +
       x5 * getter::element<6, 6>(transport_jacobian) +
       x7 * getter::element<0, 6>(transport_jacobian) +
       x9 * getter::element<1, 6>(transport_jacobian)) *
          getter::element<2, 1>(b2f_ddir_dangle);
  getter::element<1, 3>(full_jacobian) =
      x26 * getter::element<0, 1>(b2f_dpos_dangle) +
      x27 * getter::element<1, 1>(b2f_dpos_dangle) +
      x28 * getter::element<2, 1>(b2f_dpos_dangle) +
      x29 * getter::element<0, 1>(b2f_ddir_dangle) +
      x30 * getter::element<1, 1>(b2f_ddir_dangle) +
      (x20 * getter::element<4, 6>(transport_jacobian) +
       x21 * getter::element<5, 6>(transport_jacobian) +
       x22 * getter::element<6, 6>(transport_jacobian) +
       x23 * getter::element<0, 6>(transport_jacobian) +
       x24 * getter::element<1, 6>(transport_jacobian) +
       x25 * getter::element<2, 6>(transport_jacobian)) *
          getter::element<2, 1>(b2f_ddir_dangle);
  getter::element<2, 3>(full_jacobian) =
      x41 * getter::element<0, 1>(b2f_dpos_dangle) +
      x42 * getter::element<1, 1>(b2f_dpos_dangle) +
      x43 * getter::element<2, 1>(b2f_dpos_dangle) +
      x44 * getter::element<0, 1>(b2f_ddir_dangle) +
      x45 * getter::element<1, 1>(b2f_ddir_dangle) +
      (x33 * getter::element<0, 6>(transport_jacobian) +
       x34 * getter::element<1, 6>(transport_jacobian) +
       x35 * getter::element<2, 6>(transport_jacobian) +
       x36 * getter::element<6, 6>(transport_jacobian) +
       x38 * getter::element<4, 6>(transport_jacobian) +
       x40 * getter::element<5, 6>(transport_jacobian)) *
          getter::element<2, 1>(b2f_ddir_dangle);
  getter::element<3, 3>(full_jacobian) =
      x55 * getter::element<0, 1>(b2f_dpos_dangle) +
      x56 * getter::element<1, 1>(b2f_dpos_dangle) +
      x57 * getter::element<2, 1>(b2f_dpos_dangle) +
      x58 * getter::element<0, 1>(b2f_ddir_dangle) +
      x59 * getter::element<1, 1>(b2f_ddir_dangle) +
      (x49 * getter::element<0, 6>(transport_jacobian) +
       x50 * getter::element<1, 6>(transport_jacobian) +
       x51 * getter::element<2, 6>(transport_jacobian) +
       x52 * getter::element<4, 6>(transport_jacobian) +
       x53 * getter::element<5, 6>(transport_jacobian) +
       x54 * getter::element<6, 6>(transport_jacobian)) *
          getter::element<2, 1>(b2f_ddir_dangle);
  getter::element<4, 3>(full_jacobian) =
      x66 * getter::element<0, 1>(b2f_dpos_dangle) +
      x67 * getter::element<1, 1>(b2f_dpos_dangle) +
      x68 * getter::element<2, 1>(b2f_dpos_dangle) +
      x69 * getter::element<0, 1>(b2f_ddir_dangle) +
      x70 * getter::element<1, 1>(b2f_ddir_dangle) +
      (x60 * getter::element<0, 6>(transport_jacobian) +
       x61 * getter::element<1, 6>(transport_jacobian) +
       x62 * getter::element<2, 6>(transport_jacobian) +
       x63 * getter::element<4, 6>(transport_jacobian) +
       x64 * getter::element<5, 6>(transport_jacobian) +
       x65 * getter::element<6, 6>(transport_jacobian)) *
          getter::element<2, 1>(b2f_ddir_dangle);
  getter::element<5, 3>(full_jacobian) = 0;
  getter::element<0, 4>(full_jacobian) =
      x11 * getter::element<2, 7>(transport_jacobian) +
      x3 * getter::element<4, 7>(transport_jacobian) +
      x4 * getter::element<5, 7>(transport_jacobian) +
      x5 * getter::element<6, 7>(transport_jacobian) +
      x7 * getter::element<0, 7>(transport_jacobian) +
      x9 * getter::element<1, 7>(transport_jacobian);
  getter::element<1, 4>(full_jacobian) =
      x20 * getter::element<4, 7>(transport_jacobian) +
      x21 * getter::element<5, 7>(transport_jacobian) +
      x22 * getter::element<6, 7>(transport_jacobian) +
      x23 * getter::element<0, 7>(transport_jacobian) +
      x24 * getter::element<1, 7>(transport_jacobian) +
      x25 * getter::element<2, 7>(transport_jacobian);
  getter::element<2, 4>(full_jacobian) =
      x33 * getter::element<0, 7>(transport_jacobian) +
      x34 * getter::element<1, 7>(transport_jacobian) +
      x35 * getter::element<2, 7>(transport_jacobian) +
      x36 * getter::element<6, 7>(transport_jacobian) +
      x38 * getter::element<4, 7>(transport_jacobian) +
      x40 * getter::element<5, 7>(transport_jacobian);
  getter::element<3, 4>(full_jacobian) =
      x49 * getter::element<0, 7>(transport_jacobian) +
      x50 * getter::element<1, 7>(transport_jacobian) +
      x51 * getter::element<2, 7>(transport_jacobian) +
      x52 * getter::element<4, 7>(transport_jacobian) +
      x53 * getter::element<5, 7>(transport_jacobian) +
      x54 * getter::element<6, 7>(transport_jacobian);
  getter::element<4, 4>(full_jacobian) =
      x60 * getter::element<0, 7>(transport_jacobian) +
      x61 * getter::element<1, 7>(transport_jacobian) +
      x62 * getter::element<2, 7>(transport_jacobian) +
      x63 * getter::element<4, 7>(transport_jacobian) +
      x64 * getter::element<5, 7>(transport_jacobian) +
      x65 * getter::element<6, 7>(transport_jacobian) +
      getter::element<7, 7>(transport_jacobian);
  getter::element<5, 4>(full_jacobian) = 0;
  getter::element<0, 5>(full_jacobian) = 0;
  getter::element<1, 5>(full_jacobian) = 0;
  getter::element<2, 5>(full_jacobian) = 0;
  getter::element<3, 5>(full_jacobian) = 0;
  getter::element<4, 5>(full_jacobian) = 0;
  getter::element<5, 5>(full_jacobian) = 1;
}
template <typename transport_jacobian_t, typename b2f_dpos_dloc_t,
          typename b2f_ddir_dangle_t, typename b2f_dpos_dangle_t,
          typename path_to_free_derivative_t,
          typename free_to_path_derivative_t, typename f2b_dloc_dpos_t,
          typename f2b_dangle_ddir_t, typename full_jacobian_t>
DETRAY_HOST_DEVICE void inline update_full_jacobian_without_gradient_impl(
    const transport_jacobian_t& transport_jacobian,
    const b2f_dpos_dloc_t& b2f_dpos_dloc,
    const b2f_ddir_dangle_t& b2f_ddir_dangle,
    const b2f_dpos_dangle_t& b2f_dpos_dangle,
    const path_to_free_derivative_t& path_to_free_derivative,
    const free_to_path_derivative_t& free_to_path_derivative,
    const f2b_dloc_dpos_t& f2b_dloc_dpos,
    const f2b_dangle_ddir_t& f2b_dangle_ddir, full_jacobian_t& full_jacobian)
  requires((detray::concepts::square_matrix<transport_jacobian_t> &&
            detray::traits::max_rank<transport_jacobian_t> == 8) &&
           (detray::concepts::matrix<b2f_dpos_dloc_t> &&
            detray::traits::rows<b2f_dpos_dloc_t> == 3 &&
            detray::traits::columns<b2f_dpos_dloc_t> == 2) &&
           (detray::concepts::matrix<b2f_ddir_dangle_t> &&
            detray::traits::rows<b2f_ddir_dangle_t> == 3 &&
            detray::traits::columns<b2f_ddir_dangle_t> == 2) &&
           (detray::concepts::matrix<b2f_dpos_dangle_t> &&
            detray::traits::rows<b2f_dpos_dangle_t> == 3 &&
            detray::traits::columns<b2f_dpos_dangle_t> == 2) &&
           (detray::concepts::matrix<path_to_free_derivative_t> &&
            detray::traits::rows<path_to_free_derivative_t> == 8 &&
            detray::traits::columns<path_to_free_derivative_t> == 1) &&
           (detray::concepts::row_matrix<free_to_path_derivative_t> &&
            detray::traits::columns<free_to_path_derivative_t> == 8) &&
           (detray::concepts::matrix<f2b_dloc_dpos_t> &&
            detray::traits::rows<f2b_dloc_dpos_t> == 2 &&
            detray::traits::columns<f2b_dloc_dpos_t> == 3) &&
           (detray::concepts::matrix<f2b_dangle_ddir_t> &&
            detray::traits::rows<f2b_dangle_ddir_t> == 2 &&
            detray::traits::columns<f2b_dangle_ddir_t> == 3) &&
           (detray::concepts::square_matrix<full_jacobian_t> &&
            detray::traits::max_rank<full_jacobian_t> == 6))
{
  assert((getter::element<0, 0>(transport_jacobian) == 1.f));
  assert((getter::element<0, 1>(transport_jacobian) == 0.f));
  assert((getter::element<0, 2>(transport_jacobian) == 0.f));
  assert((getter::element<0, 3>(transport_jacobian) == 0.f));
  assert((getter::element<1, 0>(transport_jacobian) == 0.f));
  assert((getter::element<1, 1>(transport_jacobian) == 1.f));
  assert((getter::element<1, 2>(transport_jacobian) == 0.f));
  assert((getter::element<1, 3>(transport_jacobian) == 0.f));
  assert((getter::element<2, 0>(transport_jacobian) == 0.f));
  assert((getter::element<2, 1>(transport_jacobian) == 0.f));
  assert((getter::element<2, 2>(transport_jacobian) == 1.f));
  assert((getter::element<2, 3>(transport_jacobian) == 0.f));
  assert((getter::element<3, 0>(transport_jacobian) == 0.f));
  assert((getter::element<3, 1>(transport_jacobian) == 0.f));
  assert((getter::element<3, 2>(transport_jacobian) == 0.f));
  assert((getter::element<3, 3>(transport_jacobian) == 1.f));
  assert((getter::element<3, 4>(transport_jacobian) == 0.f));
  assert((getter::element<3, 5>(transport_jacobian) == 0.f));
  assert((getter::element<3, 6>(transport_jacobian) == 0.f));
  assert((getter::element<3, 7>(transport_jacobian) == 0.f));
  assert((getter::element<4, 0>(transport_jacobian) == 0.f));
  assert((getter::element<4, 1>(transport_jacobian) == 0.f));
  assert((getter::element<4, 2>(transport_jacobian) == 0.f));
  assert((getter::element<4, 3>(transport_jacobian) == 0.f));
  assert((getter::element<5, 0>(transport_jacobian) == 0.f));
  assert((getter::element<5, 1>(transport_jacobian) == 0.f));
  assert((getter::element<5, 2>(transport_jacobian) == 0.f));
  assert((getter::element<5, 3>(transport_jacobian) == 0.f));
  assert((getter::element<6, 0>(transport_jacobian) == 0.f));
  assert((getter::element<6, 1>(transport_jacobian) == 0.f));
  assert((getter::element<6, 2>(transport_jacobian) == 0.f));
  assert((getter::element<6, 3>(transport_jacobian) == 0.f));
  assert((getter::element<7, 0>(transport_jacobian) == 0.f));
  assert((getter::element<7, 1>(transport_jacobian) == 0.f));
  assert((getter::element<7, 2>(transport_jacobian) == 0.f));
  assert((getter::element<7, 3>(transport_jacobian) == 0.f));
  assert((getter::element<7, 4>(transport_jacobian) == 0.f));
  assert((getter::element<7, 5>(transport_jacobian) == 0.f));
  assert((getter::element<7, 6>(transport_jacobian) == 0.f));
  assert((getter::element<2, 0>(b2f_ddir_dangle) == 0.f));
  assert((getter::element<3, 0>(path_to_free_derivative) == 0.f));
  assert((getter::element<0, 3>(free_to_path_derivative) == 0.f));
  assert((getter::element<0, 7>(free_to_path_derivative) == 0.f));
  assert((getter::element<0, 2>(f2b_dangle_ddir) == 0.f));
  const auto x0 = getter::element<0, 1>(f2b_dloc_dpos) *
                  getter::element<1, 0>(path_to_free_derivative);
  const auto x1 = getter::element<0, 2>(f2b_dloc_dpos) *
                  getter::element<2, 0>(path_to_free_derivative);
  const auto x2 = getter::element<0, 0>(free_to_path_derivative) *
                      getter::element<0, 0>(path_to_free_derivative) +
                  1;
  const auto x4 = getter::element<0, 0>(f2b_dloc_dpos) *
                  getter::element<0, 0>(path_to_free_derivative);
  const auto x5 = getter::element<0, 1>(free_to_path_derivative) *
                      getter::element<1, 0>(path_to_free_derivative) +
                  1;
  const auto x7 = getter::element<0, 2>(free_to_path_derivative) *
                      getter::element<2, 0>(path_to_free_derivative) +
                  1;
  const auto x14 = getter::element<1, 1>(f2b_dloc_dpos) *
                   getter::element<1, 0>(path_to_free_derivative);
  const auto x15 = getter::element<1, 2>(f2b_dloc_dpos) *
                   getter::element<2, 0>(path_to_free_derivative);
  const auto x17 = getter::element<1, 0>(f2b_dloc_dpos) *
                   getter::element<0, 0>(path_to_free_derivative);
  const auto x25 = getter::element<0, 0>(f2b_dangle_ddir) *
                   getter::element<4, 0>(path_to_free_derivative);
  const auto x26 = getter::element<0, 1>(f2b_dangle_ddir) *
                   getter::element<5, 0>(path_to_free_derivative);
  const auto x31 = getter::element<0, 4>(free_to_path_derivative) *
                       getter::element<4, 0>(path_to_free_derivative) +
                   1;
  const auto x33 = getter::element<0, 5>(free_to_path_derivative) *
                       getter::element<5, 0>(path_to_free_derivative) +
                   1;
  const auto x37 = getter::element<1, 0>(f2b_dangle_ddir) *
                   getter::element<4, 0>(path_to_free_derivative);
  const auto x38 = getter::element<1, 1>(f2b_dangle_ddir) *
                   getter::element<5, 0>(path_to_free_derivative);
  const auto x39 = getter::element<1, 2>(f2b_dangle_ddir) *
                   getter::element<6, 0>(path_to_free_derivative);
  const auto x48 = getter::element<0, 0>(free_to_path_derivative) *
                   getter::element<7, 0>(path_to_free_derivative);
  const auto x49 = getter::element<0, 1>(free_to_path_derivative) *
                   getter::element<7, 0>(path_to_free_derivative);
  const auto x50 = getter::element<0, 2>(free_to_path_derivative) *
                   getter::element<7, 0>(path_to_free_derivative);
  const auto x51 = getter::element<0, 4>(free_to_path_derivative) *
                   getter::element<7, 0>(path_to_free_derivative);
  const auto x52 = getter::element<0, 5>(free_to_path_derivative) *
                   getter::element<7, 0>(path_to_free_derivative);
  const auto x53 = getter::element<0, 6>(free_to_path_derivative) *
                   getter::element<7, 0>(path_to_free_derivative);
  const auto x27 = (x25 + x26) * getter::element<0, 0>(free_to_path_derivative);
  const auto x28 = (x25 + x26) * getter::element<0, 1>(free_to_path_derivative);
  const auto x29 = (x25 + x26) * getter::element<0, 2>(free_to_path_derivative);
  const auto x30 = (x25 + x26) * getter::element<0, 6>(free_to_path_derivative);
  const auto x9 =
      (x0 + x1 + x4) * getter::element<0, 4>(free_to_path_derivative);
  const auto x10 =
      (x0 + x1 + x4) * getter::element<0, 5>(free_to_path_derivative);
  const auto x11 =
      (x0 + x1 + x4) * getter::element<0, 6>(free_to_path_derivative);
  const auto x20 =
      (x14 + x15 + x17) * getter::element<0, 4>(free_to_path_derivative);
  const auto x21 =
      (x14 + x15 + x17) * getter::element<0, 5>(free_to_path_derivative);
  const auto x22 =
      (x14 + x15 + x17) * getter::element<0, 6>(free_to_path_derivative);
  const auto x32 = x26 * getter::element<0, 4>(free_to_path_derivative) +
                   x31 * getter::element<0, 0>(f2b_dangle_ddir);
  const auto x34 = x25 * getter::element<0, 5>(free_to_path_derivative) +
                   x33 * getter::element<0, 1>(f2b_dangle_ddir);
  const auto x40 =
      (x37 + x38 + x39) * getter::element<0, 0>(free_to_path_derivative);
  const auto x41 =
      (x37 + x38 + x39) * getter::element<0, 1>(free_to_path_derivative);
  const auto x42 =
      (x37 + x38 + x39) * getter::element<0, 2>(free_to_path_derivative);
  const auto x3 = x0 * getter::element<0, 0>(free_to_path_derivative) +
                  x1 * getter::element<0, 0>(free_to_path_derivative) +
                  x2 * getter::element<0, 0>(f2b_dloc_dpos);
  const auto x6 = x1 * getter::element<0, 1>(free_to_path_derivative) +
                  x4 * getter::element<0, 1>(free_to_path_derivative) +
                  x5 * getter::element<0, 1>(f2b_dloc_dpos);
  const auto x8 = x0 * getter::element<0, 2>(free_to_path_derivative) +
                  x4 * getter::element<0, 2>(free_to_path_derivative) +
                  x7 * getter::element<0, 2>(f2b_dloc_dpos);
  const auto x16 = x14 * getter::element<0, 0>(free_to_path_derivative) +
                   x15 * getter::element<0, 0>(free_to_path_derivative) +
                   x2 * getter::element<1, 0>(f2b_dloc_dpos);
  const auto x18 = x15 * getter::element<0, 1>(free_to_path_derivative) +
                   x17 * getter::element<0, 1>(free_to_path_derivative) +
                   x5 * getter::element<1, 1>(f2b_dloc_dpos);
  const auto x19 = x14 * getter::element<0, 2>(free_to_path_derivative) +
                   x17 * getter::element<0, 2>(free_to_path_derivative) +
                   x7 * getter::element<1, 2>(f2b_dloc_dpos);
  const auto x43 = x31 * getter::element<1, 0>(f2b_dangle_ddir) +
                   x38 * getter::element<0, 4>(free_to_path_derivative) +
                   x39 * getter::element<0, 4>(free_to_path_derivative);
  const auto x44 = x33 * getter::element<1, 1>(f2b_dangle_ddir) +
                   x37 * getter::element<0, 5>(free_to_path_derivative) +
                   x39 * getter::element<0, 5>(free_to_path_derivative);
  const auto x45 = x37 * getter::element<0, 6>(free_to_path_derivative) +
                   x38 * getter::element<0, 6>(free_to_path_derivative) +
                   (getter::element<0, 6>(free_to_path_derivative) *
                        getter::element<6, 0>(path_to_free_derivative) +
                    1) *
                       getter::element<1, 2>(f2b_dangle_ddir);
  const auto x54 = x48 * getter::element<0, 4>(transport_jacobian) +
                   x49 * getter::element<1, 4>(transport_jacobian) +
                   x50 * getter::element<2, 4>(transport_jacobian) +
                   x51 * getter::element<4, 4>(transport_jacobian) +
                   x52 * getter::element<5, 4>(transport_jacobian) +
                   x53 * getter::element<6, 4>(transport_jacobian);
  const auto x55 = x48 * getter::element<0, 5>(transport_jacobian) +
                   x49 * getter::element<1, 5>(transport_jacobian) +
                   x50 * getter::element<2, 5>(transport_jacobian) +
                   x51 * getter::element<4, 5>(transport_jacobian) +
                   x52 * getter::element<5, 5>(transport_jacobian) +
                   x53 * getter::element<6, 5>(transport_jacobian);
  const auto x12 = x10 * getter::element<5, 4>(transport_jacobian) +
                   x11 * getter::element<6, 4>(transport_jacobian) +
                   x3 * getter::element<0, 4>(transport_jacobian) +
                   x6 * getter::element<1, 4>(transport_jacobian) +
                   x8 * getter::element<2, 4>(transport_jacobian) +
                   x9 * getter::element<4, 4>(transport_jacobian);
  const auto x13 = x10 * getter::element<5, 5>(transport_jacobian) +
                   x11 * getter::element<6, 5>(transport_jacobian) +
                   x3 * getter::element<0, 5>(transport_jacobian) +
                   x6 * getter::element<1, 5>(transport_jacobian) +
                   x8 * getter::element<2, 5>(transport_jacobian) +
                   x9 * getter::element<4, 5>(transport_jacobian);
  const auto x23 = x16 * getter::element<0, 4>(transport_jacobian) +
                   x18 * getter::element<1, 4>(transport_jacobian) +
                   x19 * getter::element<2, 4>(transport_jacobian) +
                   x20 * getter::element<4, 4>(transport_jacobian) +
                   x21 * getter::element<5, 4>(transport_jacobian) +
                   x22 * getter::element<6, 4>(transport_jacobian);
  const auto x24 = x16 * getter::element<0, 5>(transport_jacobian) +
                   x18 * getter::element<1, 5>(transport_jacobian) +
                   x19 * getter::element<2, 5>(transport_jacobian) +
                   x20 * getter::element<4, 5>(transport_jacobian) +
                   x21 * getter::element<5, 5>(transport_jacobian) +
                   x22 * getter::element<6, 5>(transport_jacobian);
  const auto x35 = x27 * getter::element<0, 4>(transport_jacobian) +
                   x28 * getter::element<1, 4>(transport_jacobian) +
                   x29 * getter::element<2, 4>(transport_jacobian) +
                   x30 * getter::element<6, 4>(transport_jacobian) +
                   x32 * getter::element<4, 4>(transport_jacobian) +
                   x34 * getter::element<5, 4>(transport_jacobian);
  const auto x36 = x27 * getter::element<0, 5>(transport_jacobian) +
                   x28 * getter::element<1, 5>(transport_jacobian) +
                   x29 * getter::element<2, 5>(transport_jacobian) +
                   x30 * getter::element<6, 5>(transport_jacobian) +
                   x32 * getter::element<4, 5>(transport_jacobian) +
                   x34 * getter::element<5, 5>(transport_jacobian);
  const auto x46 = x40 * getter::element<0, 4>(transport_jacobian) +
                   x41 * getter::element<1, 4>(transport_jacobian) +
                   x42 * getter::element<2, 4>(transport_jacobian) +
                   x43 * getter::element<4, 4>(transport_jacobian) +
                   x44 * getter::element<5, 4>(transport_jacobian) +
                   x45 * getter::element<6, 4>(transport_jacobian);
  const auto x47 = x40 * getter::element<0, 5>(transport_jacobian) +
                   x41 * getter::element<1, 5>(transport_jacobian) +
                   x42 * getter::element<2, 5>(transport_jacobian) +
                   x43 * getter::element<4, 5>(transport_jacobian) +
                   x44 * getter::element<5, 5>(transport_jacobian) +
                   x45 * getter::element<6, 5>(transport_jacobian);
  getter::element<0, 0>(full_jacobian) =
      x3 * getter::element<0, 0>(b2f_dpos_dloc) +
      x6 * getter::element<1, 0>(b2f_dpos_dloc) +
      x8 * getter::element<2, 0>(b2f_dpos_dloc);
  getter::element<1, 0>(full_jacobian) =
      x16 * getter::element<0, 0>(b2f_dpos_dloc) +
      x18 * getter::element<1, 0>(b2f_dpos_dloc) +
      x19 * getter::element<2, 0>(b2f_dpos_dloc);
  getter::element<2, 0>(full_jacobian) =
      x27 * getter::element<0, 0>(b2f_dpos_dloc) +
      x28 * getter::element<1, 0>(b2f_dpos_dloc) +
      x29 * getter::element<2, 0>(b2f_dpos_dloc);
  getter::element<3, 0>(full_jacobian) =
      x40 * getter::element<0, 0>(b2f_dpos_dloc) +
      x41 * getter::element<1, 0>(b2f_dpos_dloc) +
      x42 * getter::element<2, 0>(b2f_dpos_dloc);
  getter::element<4, 0>(full_jacobian) =
      x48 * getter::element<0, 0>(b2f_dpos_dloc) +
      x49 * getter::element<1, 0>(b2f_dpos_dloc) +
      x50 * getter::element<2, 0>(b2f_dpos_dloc);
  getter::element<5, 0>(full_jacobian) = 0;
  getter::element<0, 1>(full_jacobian) =
      x3 * getter::element<0, 1>(b2f_dpos_dloc) +
      x6 * getter::element<1, 1>(b2f_dpos_dloc) +
      x8 * getter::element<2, 1>(b2f_dpos_dloc);
  getter::element<1, 1>(full_jacobian) =
      x16 * getter::element<0, 1>(b2f_dpos_dloc) +
      x18 * getter::element<1, 1>(b2f_dpos_dloc) +
      x19 * getter::element<2, 1>(b2f_dpos_dloc);
  getter::element<2, 1>(full_jacobian) =
      x27 * getter::element<0, 1>(b2f_dpos_dloc) +
      x28 * getter::element<1, 1>(b2f_dpos_dloc) +
      x29 * getter::element<2, 1>(b2f_dpos_dloc);
  getter::element<3, 1>(full_jacobian) =
      x40 * getter::element<0, 1>(b2f_dpos_dloc) +
      x41 * getter::element<1, 1>(b2f_dpos_dloc) +
      x42 * getter::element<2, 1>(b2f_dpos_dloc);
  getter::element<4, 1>(full_jacobian) =
      x48 * getter::element<0, 1>(b2f_dpos_dloc) +
      x49 * getter::element<1, 1>(b2f_dpos_dloc) +
      x50 * getter::element<2, 1>(b2f_dpos_dloc);
  getter::element<5, 1>(full_jacobian) = 0;
  getter::element<0, 2>(full_jacobian) =
      x12 * getter::element<0, 0>(b2f_ddir_dangle) +
      x13 * getter::element<1, 0>(b2f_ddir_dangle) +
      x3 * getter::element<0, 0>(b2f_dpos_dangle) +
      x6 * getter::element<1, 0>(b2f_dpos_dangle) +
      x8 * getter::element<2, 0>(b2f_dpos_dangle);
  getter::element<1, 2>(full_jacobian) =
      x16 * getter::element<0, 0>(b2f_dpos_dangle) +
      x18 * getter::element<1, 0>(b2f_dpos_dangle) +
      x19 * getter::element<2, 0>(b2f_dpos_dangle) +
      x23 * getter::element<0, 0>(b2f_ddir_dangle) +
      x24 * getter::element<1, 0>(b2f_ddir_dangle);
  getter::element<2, 2>(full_jacobian) =
      x27 * getter::element<0, 0>(b2f_dpos_dangle) +
      x28 * getter::element<1, 0>(b2f_dpos_dangle) +
      x29 * getter::element<2, 0>(b2f_dpos_dangle) +
      x35 * getter::element<0, 0>(b2f_ddir_dangle) +
      x36 * getter::element<1, 0>(b2f_ddir_dangle);
  getter::element<3, 2>(full_jacobian) =
      x40 * getter::element<0, 0>(b2f_dpos_dangle) +
      x41 * getter::element<1, 0>(b2f_dpos_dangle) +
      x42 * getter::element<2, 0>(b2f_dpos_dangle) +
      x46 * getter::element<0, 0>(b2f_ddir_dangle) +
      x47 * getter::element<1, 0>(b2f_ddir_dangle);
  getter::element<4, 2>(full_jacobian) =
      x48 * getter::element<0, 0>(b2f_dpos_dangle) +
      x49 * getter::element<1, 0>(b2f_dpos_dangle) +
      x50 * getter::element<2, 0>(b2f_dpos_dangle) +
      x54 * getter::element<0, 0>(b2f_ddir_dangle) +
      x55 * getter::element<1, 0>(b2f_ddir_dangle);
  getter::element<5, 2>(full_jacobian) = 0;
  getter::element<0, 3>(full_jacobian) =
      x12 * getter::element<0, 1>(b2f_ddir_dangle) +
      x13 * getter::element<1, 1>(b2f_ddir_dangle) +
      x3 * getter::element<0, 1>(b2f_dpos_dangle) +
      x6 * getter::element<1, 1>(b2f_dpos_dangle) +
      x8 * getter::element<2, 1>(b2f_dpos_dangle) +
      (x10 * getter::element<5, 6>(transport_jacobian) +
       x11 * getter::element<6, 6>(transport_jacobian) +
       x3 * getter::element<0, 6>(transport_jacobian) +
       x6 * getter::element<1, 6>(transport_jacobian) +
       x8 * getter::element<2, 6>(transport_jacobian) +
       x9 * getter::element<4, 6>(transport_jacobian)) *
          getter::element<2, 1>(b2f_ddir_dangle);
  getter::element<1, 3>(full_jacobian) =
      x16 * getter::element<0, 1>(b2f_dpos_dangle) +
      x18 * getter::element<1, 1>(b2f_dpos_dangle) +
      x19 * getter::element<2, 1>(b2f_dpos_dangle) +
      x23 * getter::element<0, 1>(b2f_ddir_dangle) +
      x24 * getter::element<1, 1>(b2f_ddir_dangle) +
      (x16 * getter::element<0, 6>(transport_jacobian) +
       x18 * getter::element<1, 6>(transport_jacobian) +
       x19 * getter::element<2, 6>(transport_jacobian) +
       x20 * getter::element<4, 6>(transport_jacobian) +
       x21 * getter::element<5, 6>(transport_jacobian) +
       x22 * getter::element<6, 6>(transport_jacobian)) *
          getter::element<2, 1>(b2f_ddir_dangle);
  getter::element<2, 3>(full_jacobian) =
      x27 * getter::element<0, 1>(b2f_dpos_dangle) +
      x28 * getter::element<1, 1>(b2f_dpos_dangle) +
      x29 * getter::element<2, 1>(b2f_dpos_dangle) +
      x35 * getter::element<0, 1>(b2f_ddir_dangle) +
      x36 * getter::element<1, 1>(b2f_ddir_dangle) +
      (x27 * getter::element<0, 6>(transport_jacobian) +
       x28 * getter::element<1, 6>(transport_jacobian) +
       x29 * getter::element<2, 6>(transport_jacobian) +
       x30 * getter::element<6, 6>(transport_jacobian) +
       x32 * getter::element<4, 6>(transport_jacobian) +
       x34 * getter::element<5, 6>(transport_jacobian)) *
          getter::element<2, 1>(b2f_ddir_dangle);
  getter::element<3, 3>(full_jacobian) =
      x40 * getter::element<0, 1>(b2f_dpos_dangle) +
      x41 * getter::element<1, 1>(b2f_dpos_dangle) +
      x42 * getter::element<2, 1>(b2f_dpos_dangle) +
      x46 * getter::element<0, 1>(b2f_ddir_dangle) +
      x47 * getter::element<1, 1>(b2f_ddir_dangle) +
      (x40 * getter::element<0, 6>(transport_jacobian) +
       x41 * getter::element<1, 6>(transport_jacobian) +
       x42 * getter::element<2, 6>(transport_jacobian) +
       x43 * getter::element<4, 6>(transport_jacobian) +
       x44 * getter::element<5, 6>(transport_jacobian) +
       x45 * getter::element<6, 6>(transport_jacobian)) *
          getter::element<2, 1>(b2f_ddir_dangle);
  getter::element<4, 3>(full_jacobian) =
      x48 * getter::element<0, 1>(b2f_dpos_dangle) +
      x49 * getter::element<1, 1>(b2f_dpos_dangle) +
      x50 * getter::element<2, 1>(b2f_dpos_dangle) +
      x54 * getter::element<0, 1>(b2f_ddir_dangle) +
      x55 * getter::element<1, 1>(b2f_ddir_dangle) +
      (x48 * getter::element<0, 6>(transport_jacobian) +
       x49 * getter::element<1, 6>(transport_jacobian) +
       x50 * getter::element<2, 6>(transport_jacobian) +
       x51 * getter::element<4, 6>(transport_jacobian) +
       x52 * getter::element<5, 6>(transport_jacobian) +
       x53 * getter::element<6, 6>(transport_jacobian)) *
          getter::element<2, 1>(b2f_ddir_dangle);
  getter::element<5, 3>(full_jacobian) = 0;
  getter::element<0, 4>(full_jacobian) =
      x10 * getter::element<5, 7>(transport_jacobian) +
      x11 * getter::element<6, 7>(transport_jacobian) +
      x3 * getter::element<0, 7>(transport_jacobian) +
      x6 * getter::element<1, 7>(transport_jacobian) +
      x8 * getter::element<2, 7>(transport_jacobian) +
      x9 * getter::element<4, 7>(transport_jacobian);
  getter::element<1, 4>(full_jacobian) =
      x16 * getter::element<0, 7>(transport_jacobian) +
      x18 * getter::element<1, 7>(transport_jacobian) +
      x19 * getter::element<2, 7>(transport_jacobian) +
      x20 * getter::element<4, 7>(transport_jacobian) +
      x21 * getter::element<5, 7>(transport_jacobian) +
      x22 * getter::element<6, 7>(transport_jacobian);
  getter::element<2, 4>(full_jacobian) =
      x27 * getter::element<0, 7>(transport_jacobian) +
      x28 * getter::element<1, 7>(transport_jacobian) +
      x29 * getter::element<2, 7>(transport_jacobian) +
      x30 * getter::element<6, 7>(transport_jacobian) +
      x32 * getter::element<4, 7>(transport_jacobian) +
      x34 * getter::element<5, 7>(transport_jacobian);
  getter::element<3, 4>(full_jacobian) =
      x40 * getter::element<0, 7>(transport_jacobian) +
      x41 * getter::element<1, 7>(transport_jacobian) +
      x42 * getter::element<2, 7>(transport_jacobian) +
      x43 * getter::element<4, 7>(transport_jacobian) +
      x44 * getter::element<5, 7>(transport_jacobian) +
      x45 * getter::element<6, 7>(transport_jacobian);
  getter::element<4, 4>(full_jacobian) =
      x48 * getter::element<0, 7>(transport_jacobian) +
      x49 * getter::element<1, 7>(transport_jacobian) +
      x50 * getter::element<2, 7>(transport_jacobian) +
      x51 * getter::element<4, 7>(transport_jacobian) +
      x52 * getter::element<5, 7>(transport_jacobian) +
      x53 * getter::element<6, 7>(transport_jacobian) +
      getter::element<7, 7>(transport_jacobian);
  getter::element<5, 4>(full_jacobian) = 0;
  getter::element<0, 5>(full_jacobian) = 0;
  getter::element<1, 5>(full_jacobian) = 0;
  getter::element<2, 5>(full_jacobian) = 0;
  getter::element<3, 5>(full_jacobian) = 0;
  getter::element<4, 5>(full_jacobian) = 0;
  getter::element<5, 5>(full_jacobian) = 1;
}
}  // namespace detray::detail
