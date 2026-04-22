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
template <typename J_transport_t, typename dFdt_t, typename dGdt_t,
          typename dFdqop_t, typename dGdqop_t, typename dqopqop_t,
          typename new_J_t>
DETRAY_HOST_DEVICE void inline update_transport_jacobian_without_gradient_impl(
    const J_transport_t& J_transport, const dFdt_t& dFdt, const dGdt_t& dGdt,
    const dFdqop_t& dFdqop, const dGdqop_t& dGdqop, const dqopqop_t& dqopqop,
    new_J_t& new_J)
  requires((detray::concepts::square_matrix<J_transport_t> &&
            detray::traits::max_rank<J_transport_t> == 8) &&
           (detray::concepts::square_matrix<dFdt_t> &&
            detray::traits::max_rank<dFdt_t> == 3) &&
           (detray::concepts::square_matrix<dGdt_t> &&
            detray::traits::max_rank<dGdt_t> == 3) &&
           (detray::concepts::vector3D<dFdqop_t>) &&
           (detray::concepts::vector3D<dGdqop_t>) &&
           (detray::concepts::scalar<dqopqop_t>) &&
           (detray::concepts::square_matrix<new_J_t> &&
            detray::traits::max_rank<new_J_t> == 8))
{
  assert((getter::element<0, 0>(J_transport) == 1.f));
  assert((getter::element<0, 1>(J_transport) == 0.f));
  assert((getter::element<0, 2>(J_transport) == 0.f));
  assert((getter::element<0, 3>(J_transport) == 0.f));
  assert((getter::element<1, 0>(J_transport) == 0.f));
  assert((getter::element<1, 1>(J_transport) == 1.f));
  assert((getter::element<1, 2>(J_transport) == 0.f));
  assert((getter::element<1, 3>(J_transport) == 0.f));
  assert((getter::element<2, 0>(J_transport) == 0.f));
  assert((getter::element<2, 1>(J_transport) == 0.f));
  assert((getter::element<2, 2>(J_transport) == 1.f));
  assert((getter::element<2, 3>(J_transport) == 0.f));
  assert((getter::element<3, 0>(J_transport) == 0.f));
  assert((getter::element<3, 1>(J_transport) == 0.f));
  assert((getter::element<3, 2>(J_transport) == 0.f));
  assert((getter::element<3, 3>(J_transport) == 1.f));
  assert((getter::element<3, 4>(J_transport) == 0.f));
  assert((getter::element<3, 5>(J_transport) == 0.f));
  assert((getter::element<3, 6>(J_transport) == 0.f));
  assert((getter::element<3, 7>(J_transport) == 0.f));
  assert((getter::element<4, 0>(J_transport) == 0.f));
  assert((getter::element<4, 1>(J_transport) == 0.f));
  assert((getter::element<4, 2>(J_transport) == 0.f));
  assert((getter::element<4, 3>(J_transport) == 0.f));
  assert((getter::element<5, 0>(J_transport) == 0.f));
  assert((getter::element<5, 1>(J_transport) == 0.f));
  assert((getter::element<5, 2>(J_transport) == 0.f));
  assert((getter::element<5, 3>(J_transport) == 0.f));
  assert((getter::element<6, 0>(J_transport) == 0.f));
  assert((getter::element<6, 1>(J_transport) == 0.f));
  assert((getter::element<6, 2>(J_transport) == 0.f));
  assert((getter::element<6, 3>(J_transport) == 0.f));
  assert((getter::element<7, 0>(J_transport) == 0.f));
  assert((getter::element<7, 1>(J_transport) == 0.f));
  assert((getter::element<7, 2>(J_transport) == 0.f));
  assert((getter::element<7, 3>(J_transport) == 0.f));
  assert((getter::element<7, 4>(J_transport) == 0.f));
  assert((getter::element<7, 5>(J_transport) == 0.f));
  assert((getter::element<7, 6>(J_transport) == 0.f));

  getter::element<0, 4>(new_J) =
      getter::element<0, 4>(J_transport) +
      getter::element<4, 4>(J_transport) * getter::element<0, 0>(dFdt) +
      getter::element<5, 4>(J_transport) * getter::element<0, 1>(dFdt) +
      getter::element<6, 4>(J_transport) * getter::element<0, 2>(dFdt);
  getter::element<1, 4>(new_J) =
      getter::element<1, 4>(J_transport) +
      getter::element<4, 4>(J_transport) * getter::element<1, 0>(dFdt) +
      getter::element<5, 4>(J_transport) * getter::element<1, 1>(dFdt) +
      getter::element<6, 4>(J_transport) * getter::element<1, 2>(dFdt);
  getter::element<2, 4>(new_J) =
      getter::element<2, 4>(J_transport) +
      getter::element<4, 4>(J_transport) * getter::element<2, 0>(dFdt) +
      getter::element<5, 4>(J_transport) * getter::element<2, 1>(dFdt) +
      getter::element<6, 4>(J_transport) * getter::element<2, 2>(dFdt);

  getter::element<4, 4>(new_J) =
      getter::element<4, 4>(J_transport) * getter::element<0, 0>(dGdt) +
      getter::element<5, 4>(J_transport) * getter::element<0, 1>(dGdt) +
      getter::element<6, 4>(J_transport) * getter::element<0, 2>(dGdt);
  getter::element<5, 4>(new_J) =
      getter::element<4, 4>(J_transport) * getter::element<1, 0>(dGdt) +
      getter::element<5, 4>(J_transport) * getter::element<1, 1>(dGdt) +
      getter::element<6, 4>(J_transport) * getter::element<1, 2>(dGdt);
  getter::element<6, 4>(new_J) =
      getter::element<4, 4>(J_transport) * getter::element<2, 0>(dGdt) +
      getter::element<5, 4>(J_transport) * getter::element<2, 1>(dGdt) +
      getter::element<6, 4>(J_transport) * getter::element<2, 2>(dGdt);

  getter::element<0, 5>(new_J) =
      getter::element<0, 5>(J_transport) +
      getter::element<4, 5>(J_transport) * getter::element<0, 0>(dFdt) +
      getter::element<5, 5>(J_transport) * getter::element<0, 1>(dFdt) +
      getter::element<6, 5>(J_transport) * getter::element<0, 2>(dFdt);
  getter::element<1, 5>(new_J) =
      getter::element<1, 5>(J_transport) +
      getter::element<4, 5>(J_transport) * getter::element<1, 0>(dFdt) +
      getter::element<5, 5>(J_transport) * getter::element<1, 1>(dFdt) +
      getter::element<6, 5>(J_transport) * getter::element<1, 2>(dFdt);
  getter::element<2, 5>(new_J) =
      getter::element<2, 5>(J_transport) +
      getter::element<4, 5>(J_transport) * getter::element<2, 0>(dFdt) +
      getter::element<5, 5>(J_transport) * getter::element<2, 1>(dFdt) +
      getter::element<6, 5>(J_transport) * getter::element<2, 2>(dFdt);

  getter::element<4, 5>(new_J) =
      getter::element<4, 5>(J_transport) * getter::element<0, 0>(dGdt) +
      getter::element<5, 5>(J_transport) * getter::element<0, 1>(dGdt) +
      getter::element<6, 5>(J_transport) * getter::element<0, 2>(dGdt);
  getter::element<5, 5>(new_J) =
      getter::element<4, 5>(J_transport) * getter::element<1, 0>(dGdt) +
      getter::element<5, 5>(J_transport) * getter::element<1, 1>(dGdt) +
      getter::element<6, 5>(J_transport) * getter::element<1, 2>(dGdt);
  getter::element<6, 5>(new_J) =
      getter::element<4, 5>(J_transport) * getter::element<2, 0>(dGdt) +
      getter::element<5, 5>(J_transport) * getter::element<2, 1>(dGdt) +
      getter::element<6, 5>(J_transport) * getter::element<2, 2>(dGdt);

  getter::element<0, 6>(new_J) =
      getter::element<0, 6>(J_transport) +
      getter::element<4, 6>(J_transport) * getter::element<0, 0>(dFdt) +
      getter::element<5, 6>(J_transport) * getter::element<0, 1>(dFdt) +
      getter::element<6, 6>(J_transport) * getter::element<0, 2>(dFdt);
  getter::element<1, 6>(new_J) =
      getter::element<1, 6>(J_transport) +
      getter::element<4, 6>(J_transport) * getter::element<1, 0>(dFdt) +
      getter::element<5, 6>(J_transport) * getter::element<1, 1>(dFdt) +
      getter::element<6, 6>(J_transport) * getter::element<1, 2>(dFdt);
  getter::element<2, 6>(new_J) =
      getter::element<2, 6>(J_transport) +
      getter::element<4, 6>(J_transport) * getter::element<2, 0>(dFdt) +
      getter::element<5, 6>(J_transport) * getter::element<2, 1>(dFdt) +
      getter::element<6, 6>(J_transport) * getter::element<2, 2>(dFdt);

  getter::element<4, 6>(new_J) =
      getter::element<4, 6>(J_transport) * getter::element<0, 0>(dGdt) +
      getter::element<5, 6>(J_transport) * getter::element<0, 1>(dGdt) +
      getter::element<6, 6>(J_transport) * getter::element<0, 2>(dGdt);
  getter::element<5, 6>(new_J) =
      getter::element<4, 6>(J_transport) * getter::element<1, 0>(dGdt) +
      getter::element<5, 6>(J_transport) * getter::element<1, 1>(dGdt) +
      getter::element<6, 6>(J_transport) * getter::element<1, 2>(dGdt);
  getter::element<6, 6>(new_J) =
      getter::element<4, 6>(J_transport) * getter::element<2, 0>(dGdt) +
      getter::element<5, 6>(J_transport) * getter::element<2, 1>(dGdt) +
      getter::element<6, 6>(J_transport) * getter::element<2, 2>(dGdt);

  getter::element<0, 7>(new_J) =
      getter::element<0, 7>(J_transport) +
      getter::element<4, 7>(J_transport) * getter::element<0, 0>(dFdt) +
      getter::element<5, 7>(J_transport) * getter::element<0, 1>(dFdt) +
      getter::element<6, 7>(J_transport) * getter::element<0, 2>(dFdt) +
      getter::element<7, 7>(J_transport) * getter::element<0>(dFdqop);
  getter::element<1, 7>(new_J) =
      getter::element<1, 7>(J_transport) +
      getter::element<4, 7>(J_transport) * getter::element<1, 0>(dFdt) +
      getter::element<5, 7>(J_transport) * getter::element<1, 1>(dFdt) +
      getter::element<6, 7>(J_transport) * getter::element<1, 2>(dFdt) +
      getter::element<7, 7>(J_transport) * getter::element<1>(dFdqop);
  getter::element<2, 7>(new_J) =
      getter::element<2, 7>(J_transport) +
      getter::element<4, 7>(J_transport) * getter::element<2, 0>(dFdt) +
      getter::element<5, 7>(J_transport) * getter::element<2, 1>(dFdt) +
      getter::element<6, 7>(J_transport) * getter::element<2, 2>(dFdt) +
      getter::element<7, 7>(J_transport) * getter::element<2>(dFdqop);

  getter::element<4, 7>(new_J) =
      getter::element<4, 7>(J_transport) * getter::element<0, 0>(dGdt) +
      getter::element<5, 7>(J_transport) * getter::element<0, 1>(dGdt) +
      getter::element<6, 7>(J_transport) * getter::element<0, 2>(dGdt) +
      getter::element<7, 7>(J_transport) * getter::element<0>(dGdqop);
  getter::element<5, 7>(new_J) =
      getter::element<4, 7>(J_transport) * getter::element<1, 0>(dGdt) +
      getter::element<5, 7>(J_transport) * getter::element<1, 1>(dGdt) +
      getter::element<6, 7>(J_transport) * getter::element<1, 2>(dGdt) +
      getter::element<7, 7>(J_transport) * getter::element<1>(dGdqop);
  getter::element<6, 7>(new_J) =
      getter::element<4, 7>(J_transport) * getter::element<2, 0>(dGdt) +
      getter::element<5, 7>(J_transport) * getter::element<2, 1>(dGdt) +
      getter::element<6, 7>(J_transport) * getter::element<2, 2>(dGdt) +
      getter::element<7, 7>(J_transport) * getter::element<2>(dGdqop);
  getter::element<7, 7>(new_J) = dqopqop * getter::element<7, 7>(J_transport);
}
template <typename J_transport_t, typename dFdt_t, typename dGdt_t,
          typename dFdr_t, typename dGdr_t, typename dFdqop_t,
          typename dGdqop_t, typename dqopqop_t, typename new_J_t>
DETRAY_HOST_DEVICE void inline update_transport_jacobian_with_gradient_impl(
    const J_transport_t& J_transport, const dFdt_t& dFdt, const dGdt_t& dGdt,
    const dFdr_t& dFdr, const dGdr_t& dGdr, const dFdqop_t& dFdqop,
    const dGdqop_t& dGdqop, const dqopqop_t& dqopqop, new_J_t& new_J)
  requires((detray::concepts::square_matrix<J_transport_t> &&
            detray::traits::max_rank<J_transport_t> == 8) &&
           (detray::concepts::square_matrix<dFdt_t> &&
            detray::traits::max_rank<dFdt_t> == 3) &&
           (detray::concepts::square_matrix<dGdt_t> &&
            detray::traits::max_rank<dGdt_t> == 3) &&
           (detray::concepts::square_matrix<dFdr_t> &&
            detray::traits::max_rank<dFdr_t> == 3) &&
           (detray::concepts::square_matrix<dGdr_t> &&
            detray::traits::max_rank<dGdr_t> == 3) &&
           (detray::concepts::vector3D<dFdqop_t>) &&
           (detray::concepts::vector3D<dGdqop_t>) &&
           (detray::concepts::scalar<dqopqop_t>) &&
           (detray::concepts::square_matrix<new_J_t> &&
            detray::traits::max_rank<new_J_t> == 8))
{
  assert((getter::element<0, 3>(J_transport) == 0.f));
  assert((getter::element<1, 3>(J_transport) == 0.f));
  assert((getter::element<2, 3>(J_transport) == 0.f));
  assert((getter::element<3, 0>(J_transport) == 0.f));
  assert((getter::element<3, 1>(J_transport) == 0.f));
  assert((getter::element<3, 2>(J_transport) == 0.f));
  assert((getter::element<3, 3>(J_transport) == 1.f));
  assert((getter::element<3, 4>(J_transport) == 0.f));
  assert((getter::element<3, 5>(J_transport) == 0.f));
  assert((getter::element<3, 6>(J_transport) == 0.f));
  assert((getter::element<3, 7>(J_transport) == 0.f));
  assert((getter::element<4, 3>(J_transport) == 0.f));
  assert((getter::element<5, 3>(J_transport) == 0.f));
  assert((getter::element<6, 3>(J_transport) == 0.f));
  assert((getter::element<7, 0>(J_transport) == 0.f));
  assert((getter::element<7, 1>(J_transport) == 0.f));
  assert((getter::element<7, 2>(J_transport) == 0.f));
  assert((getter::element<7, 3>(J_transport) == 0.f));
  assert((getter::element<7, 4>(J_transport) == 0.f));
  assert((getter::element<7, 5>(J_transport) == 0.f));
  assert((getter::element<7, 6>(J_transport) == 0.f));
  getter::element<0, 0>(new_J) =
      getter::element<0, 0>(J_transport) * getter::element<0, 0>(dFdr) +
      getter::element<1, 0>(J_transport) * getter::element<0, 1>(dFdr) +
      getter::element<2, 0>(J_transport) * getter::element<0, 2>(dFdr) +
      getter::element<4, 0>(J_transport) * getter::element<0, 0>(dFdt) +
      getter::element<5, 0>(J_transport) * getter::element<0, 1>(dFdt) +
      getter::element<6, 0>(J_transport) * getter::element<0, 2>(dFdt);
  getter::element<1, 0>(new_J) =
      getter::element<0, 0>(J_transport) * getter::element<1, 0>(dFdr) +
      getter::element<1, 0>(J_transport) * getter::element<1, 1>(dFdr) +
      getter::element<2, 0>(J_transport) * getter::element<1, 2>(dFdr) +
      getter::element<4, 0>(J_transport) * getter::element<1, 0>(dFdt) +
      getter::element<5, 0>(J_transport) * getter::element<1, 1>(dFdt) +
      getter::element<6, 0>(J_transport) * getter::element<1, 2>(dFdt);
  getter::element<2, 0>(new_J) =
      getter::element<0, 0>(J_transport) * getter::element<2, 0>(dFdr) +
      getter::element<1, 0>(J_transport) * getter::element<2, 1>(dFdr) +
      getter::element<2, 0>(J_transport) * getter::element<2, 2>(dFdr) +
      getter::element<4, 0>(J_transport) * getter::element<2, 0>(dFdt) +
      getter::element<5, 0>(J_transport) * getter::element<2, 1>(dFdt) +
      getter::element<6, 0>(J_transport) * getter::element<2, 2>(dFdt);

  getter::element<4, 0>(new_J) =
      getter::element<0, 0>(J_transport) * getter::element<0, 0>(dGdr) +
      getter::element<1, 0>(J_transport) * getter::element<0, 1>(dGdr) +
      getter::element<2, 0>(J_transport) * getter::element<0, 2>(dGdr) +
      getter::element<4, 0>(J_transport) * getter::element<0, 0>(dGdt) +
      getter::element<5, 0>(J_transport) * getter::element<0, 1>(dGdt) +
      getter::element<6, 0>(J_transport) * getter::element<0, 2>(dGdt);
  getter::element<5, 0>(new_J) =
      getter::element<0, 0>(J_transport) * getter::element<1, 0>(dGdr) +
      getter::element<1, 0>(J_transport) * getter::element<1, 1>(dGdr) +
      getter::element<2, 0>(J_transport) * getter::element<1, 2>(dGdr) +
      getter::element<4, 0>(J_transport) * getter::element<1, 0>(dGdt) +
      getter::element<5, 0>(J_transport) * getter::element<1, 1>(dGdt) +
      getter::element<6, 0>(J_transport) * getter::element<1, 2>(dGdt);
  getter::element<6, 0>(new_J) =
      getter::element<0, 0>(J_transport) * getter::element<2, 0>(dGdr) +
      getter::element<1, 0>(J_transport) * getter::element<2, 1>(dGdr) +
      getter::element<2, 0>(J_transport) * getter::element<2, 2>(dGdr) +
      getter::element<4, 0>(J_transport) * getter::element<2, 0>(dGdt) +
      getter::element<5, 0>(J_transport) * getter::element<2, 1>(dGdt) +
      getter::element<6, 0>(J_transport) * getter::element<2, 2>(dGdt);

  getter::element<0, 1>(new_J) =
      getter::element<0, 1>(J_transport) * getter::element<0, 0>(dFdr) +
      getter::element<1, 1>(J_transport) * getter::element<0, 1>(dFdr) +
      getter::element<2, 1>(J_transport) * getter::element<0, 2>(dFdr) +
      getter::element<4, 1>(J_transport) * getter::element<0, 0>(dFdt) +
      getter::element<5, 1>(J_transport) * getter::element<0, 1>(dFdt) +
      getter::element<6, 1>(J_transport) * getter::element<0, 2>(dFdt);
  getter::element<1, 1>(new_J) =
      getter::element<0, 1>(J_transport) * getter::element<1, 0>(dFdr) +
      getter::element<1, 1>(J_transport) * getter::element<1, 1>(dFdr) +
      getter::element<2, 1>(J_transport) * getter::element<1, 2>(dFdr) +
      getter::element<4, 1>(J_transport) * getter::element<1, 0>(dFdt) +
      getter::element<5, 1>(J_transport) * getter::element<1, 1>(dFdt) +
      getter::element<6, 1>(J_transport) * getter::element<1, 2>(dFdt);
  getter::element<2, 1>(new_J) =
      getter::element<0, 1>(J_transport) * getter::element<2, 0>(dFdr) +
      getter::element<1, 1>(J_transport) * getter::element<2, 1>(dFdr) +
      getter::element<2, 1>(J_transport) * getter::element<2, 2>(dFdr) +
      getter::element<4, 1>(J_transport) * getter::element<2, 0>(dFdt) +
      getter::element<5, 1>(J_transport) * getter::element<2, 1>(dFdt) +
      getter::element<6, 1>(J_transport) * getter::element<2, 2>(dFdt);

  getter::element<4, 1>(new_J) =
      getter::element<0, 1>(J_transport) * getter::element<0, 0>(dGdr) +
      getter::element<1, 1>(J_transport) * getter::element<0, 1>(dGdr) +
      getter::element<2, 1>(J_transport) * getter::element<0, 2>(dGdr) +
      getter::element<4, 1>(J_transport) * getter::element<0, 0>(dGdt) +
      getter::element<5, 1>(J_transport) * getter::element<0, 1>(dGdt) +
      getter::element<6, 1>(J_transport) * getter::element<0, 2>(dGdt);
  getter::element<5, 1>(new_J) =
      getter::element<0, 1>(J_transport) * getter::element<1, 0>(dGdr) +
      getter::element<1, 1>(J_transport) * getter::element<1, 1>(dGdr) +
      getter::element<2, 1>(J_transport) * getter::element<1, 2>(dGdr) +
      getter::element<4, 1>(J_transport) * getter::element<1, 0>(dGdt) +
      getter::element<5, 1>(J_transport) * getter::element<1, 1>(dGdt) +
      getter::element<6, 1>(J_transport) * getter::element<1, 2>(dGdt);
  getter::element<6, 1>(new_J) =
      getter::element<0, 1>(J_transport) * getter::element<2, 0>(dGdr) +
      getter::element<1, 1>(J_transport) * getter::element<2, 1>(dGdr) +
      getter::element<2, 1>(J_transport) * getter::element<2, 2>(dGdr) +
      getter::element<4, 1>(J_transport) * getter::element<2, 0>(dGdt) +
      getter::element<5, 1>(J_transport) * getter::element<2, 1>(dGdt) +
      getter::element<6, 1>(J_transport) * getter::element<2, 2>(dGdt);

  getter::element<0, 2>(new_J) =
      getter::element<0, 2>(J_transport) * getter::element<0, 0>(dFdr) +
      getter::element<1, 2>(J_transport) * getter::element<0, 1>(dFdr) +
      getter::element<2, 2>(J_transport) * getter::element<0, 2>(dFdr) +
      getter::element<4, 2>(J_transport) * getter::element<0, 0>(dFdt) +
      getter::element<5, 2>(J_transport) * getter::element<0, 1>(dFdt) +
      getter::element<6, 2>(J_transport) * getter::element<0, 2>(dFdt);
  getter::element<1, 2>(new_J) =
      getter::element<0, 2>(J_transport) * getter::element<1, 0>(dFdr) +
      getter::element<1, 2>(J_transport) * getter::element<1, 1>(dFdr) +
      getter::element<2, 2>(J_transport) * getter::element<1, 2>(dFdr) +
      getter::element<4, 2>(J_transport) * getter::element<1, 0>(dFdt) +
      getter::element<5, 2>(J_transport) * getter::element<1, 1>(dFdt) +
      getter::element<6, 2>(J_transport) * getter::element<1, 2>(dFdt);
  getter::element<2, 2>(new_J) =
      getter::element<0, 2>(J_transport) * getter::element<2, 0>(dFdr) +
      getter::element<1, 2>(J_transport) * getter::element<2, 1>(dFdr) +
      getter::element<2, 2>(J_transport) * getter::element<2, 2>(dFdr) +
      getter::element<4, 2>(J_transport) * getter::element<2, 0>(dFdt) +
      getter::element<5, 2>(J_transport) * getter::element<2, 1>(dFdt) +
      getter::element<6, 2>(J_transport) * getter::element<2, 2>(dFdt);

  getter::element<4, 2>(new_J) =
      getter::element<0, 2>(J_transport) * getter::element<0, 0>(dGdr) +
      getter::element<1, 2>(J_transport) * getter::element<0, 1>(dGdr) +
      getter::element<2, 2>(J_transport) * getter::element<0, 2>(dGdr) +
      getter::element<4, 2>(J_transport) * getter::element<0, 0>(dGdt) +
      getter::element<5, 2>(J_transport) * getter::element<0, 1>(dGdt) +
      getter::element<6, 2>(J_transport) * getter::element<0, 2>(dGdt);
  getter::element<5, 2>(new_J) =
      getter::element<0, 2>(J_transport) * getter::element<1, 0>(dGdr) +
      getter::element<1, 2>(J_transport) * getter::element<1, 1>(dGdr) +
      getter::element<2, 2>(J_transport) * getter::element<1, 2>(dGdr) +
      getter::element<4, 2>(J_transport) * getter::element<1, 0>(dGdt) +
      getter::element<5, 2>(J_transport) * getter::element<1, 1>(dGdt) +
      getter::element<6, 2>(J_transport) * getter::element<1, 2>(dGdt);
  getter::element<6, 2>(new_J) =
      getter::element<0, 2>(J_transport) * getter::element<2, 0>(dGdr) +
      getter::element<1, 2>(J_transport) * getter::element<2, 1>(dGdr) +
      getter::element<2, 2>(J_transport) * getter::element<2, 2>(dGdr) +
      getter::element<4, 2>(J_transport) * getter::element<2, 0>(dGdt) +
      getter::element<5, 2>(J_transport) * getter::element<2, 1>(dGdt) +
      getter::element<6, 2>(J_transport) * getter::element<2, 2>(dGdt);

  getter::element<0, 4>(new_J) =
      getter::element<0, 4>(J_transport) * getter::element<0, 0>(dFdr) +
      getter::element<1, 4>(J_transport) * getter::element<0, 1>(dFdr) +
      getter::element<2, 4>(J_transport) * getter::element<0, 2>(dFdr) +
      getter::element<4, 4>(J_transport) * getter::element<0, 0>(dFdt) +
      getter::element<5, 4>(J_transport) * getter::element<0, 1>(dFdt) +
      getter::element<6, 4>(J_transport) * getter::element<0, 2>(dFdt);
  getter::element<1, 4>(new_J) =
      getter::element<0, 4>(J_transport) * getter::element<1, 0>(dFdr) +
      getter::element<1, 4>(J_transport) * getter::element<1, 1>(dFdr) +
      getter::element<2, 4>(J_transport) * getter::element<1, 2>(dFdr) +
      getter::element<4, 4>(J_transport) * getter::element<1, 0>(dFdt) +
      getter::element<5, 4>(J_transport) * getter::element<1, 1>(dFdt) +
      getter::element<6, 4>(J_transport) * getter::element<1, 2>(dFdt);
  getter::element<2, 4>(new_J) =
      getter::element<0, 4>(J_transport) * getter::element<2, 0>(dFdr) +
      getter::element<1, 4>(J_transport) * getter::element<2, 1>(dFdr) +
      getter::element<2, 4>(J_transport) * getter::element<2, 2>(dFdr) +
      getter::element<4, 4>(J_transport) * getter::element<2, 0>(dFdt) +
      getter::element<5, 4>(J_transport) * getter::element<2, 1>(dFdt) +
      getter::element<6, 4>(J_transport) * getter::element<2, 2>(dFdt);

  getter::element<4, 4>(new_J) =
      getter::element<0, 4>(J_transport) * getter::element<0, 0>(dGdr) +
      getter::element<1, 4>(J_transport) * getter::element<0, 1>(dGdr) +
      getter::element<2, 4>(J_transport) * getter::element<0, 2>(dGdr) +
      getter::element<4, 4>(J_transport) * getter::element<0, 0>(dGdt) +
      getter::element<5, 4>(J_transport) * getter::element<0, 1>(dGdt) +
      getter::element<6, 4>(J_transport) * getter::element<0, 2>(dGdt);
  getter::element<5, 4>(new_J) =
      getter::element<0, 4>(J_transport) * getter::element<1, 0>(dGdr) +
      getter::element<1, 4>(J_transport) * getter::element<1, 1>(dGdr) +
      getter::element<2, 4>(J_transport) * getter::element<1, 2>(dGdr) +
      getter::element<4, 4>(J_transport) * getter::element<1, 0>(dGdt) +
      getter::element<5, 4>(J_transport) * getter::element<1, 1>(dGdt) +
      getter::element<6, 4>(J_transport) * getter::element<1, 2>(dGdt);
  getter::element<6, 4>(new_J) =
      getter::element<0, 4>(J_transport) * getter::element<2, 0>(dGdr) +
      getter::element<1, 4>(J_transport) * getter::element<2, 1>(dGdr) +
      getter::element<2, 4>(J_transport) * getter::element<2, 2>(dGdr) +
      getter::element<4, 4>(J_transport) * getter::element<2, 0>(dGdt) +
      getter::element<5, 4>(J_transport) * getter::element<2, 1>(dGdt) +
      getter::element<6, 4>(J_transport) * getter::element<2, 2>(dGdt);

  getter::element<0, 5>(new_J) =
      getter::element<0, 5>(J_transport) * getter::element<0, 0>(dFdr) +
      getter::element<1, 5>(J_transport) * getter::element<0, 1>(dFdr) +
      getter::element<2, 5>(J_transport) * getter::element<0, 2>(dFdr) +
      getter::element<4, 5>(J_transport) * getter::element<0, 0>(dFdt) +
      getter::element<5, 5>(J_transport) * getter::element<0, 1>(dFdt) +
      getter::element<6, 5>(J_transport) * getter::element<0, 2>(dFdt);
  getter::element<1, 5>(new_J) =
      getter::element<0, 5>(J_transport) * getter::element<1, 0>(dFdr) +
      getter::element<1, 5>(J_transport) * getter::element<1, 1>(dFdr) +
      getter::element<2, 5>(J_transport) * getter::element<1, 2>(dFdr) +
      getter::element<4, 5>(J_transport) * getter::element<1, 0>(dFdt) +
      getter::element<5, 5>(J_transport) * getter::element<1, 1>(dFdt) +
      getter::element<6, 5>(J_transport) * getter::element<1, 2>(dFdt);
  getter::element<2, 5>(new_J) =
      getter::element<0, 5>(J_transport) * getter::element<2, 0>(dFdr) +
      getter::element<1, 5>(J_transport) * getter::element<2, 1>(dFdr) +
      getter::element<2, 5>(J_transport) * getter::element<2, 2>(dFdr) +
      getter::element<4, 5>(J_transport) * getter::element<2, 0>(dFdt) +
      getter::element<5, 5>(J_transport) * getter::element<2, 1>(dFdt) +
      getter::element<6, 5>(J_transport) * getter::element<2, 2>(dFdt);

  getter::element<4, 5>(new_J) =
      getter::element<0, 5>(J_transport) * getter::element<0, 0>(dGdr) +
      getter::element<1, 5>(J_transport) * getter::element<0, 1>(dGdr) +
      getter::element<2, 5>(J_transport) * getter::element<0, 2>(dGdr) +
      getter::element<4, 5>(J_transport) * getter::element<0, 0>(dGdt) +
      getter::element<5, 5>(J_transport) * getter::element<0, 1>(dGdt) +
      getter::element<6, 5>(J_transport) * getter::element<0, 2>(dGdt);
  getter::element<5, 5>(new_J) =
      getter::element<0, 5>(J_transport) * getter::element<1, 0>(dGdr) +
      getter::element<1, 5>(J_transport) * getter::element<1, 1>(dGdr) +
      getter::element<2, 5>(J_transport) * getter::element<1, 2>(dGdr) +
      getter::element<4, 5>(J_transport) * getter::element<1, 0>(dGdt) +
      getter::element<5, 5>(J_transport) * getter::element<1, 1>(dGdt) +
      getter::element<6, 5>(J_transport) * getter::element<1, 2>(dGdt);
  getter::element<6, 5>(new_J) =
      getter::element<0, 5>(J_transport) * getter::element<2, 0>(dGdr) +
      getter::element<1, 5>(J_transport) * getter::element<2, 1>(dGdr) +
      getter::element<2, 5>(J_transport) * getter::element<2, 2>(dGdr) +
      getter::element<4, 5>(J_transport) * getter::element<2, 0>(dGdt) +
      getter::element<5, 5>(J_transport) * getter::element<2, 1>(dGdt) +
      getter::element<6, 5>(J_transport) * getter::element<2, 2>(dGdt);

  getter::element<0, 6>(new_J) =
      getter::element<0, 6>(J_transport) * getter::element<0, 0>(dFdr) +
      getter::element<1, 6>(J_transport) * getter::element<0, 1>(dFdr) +
      getter::element<2, 6>(J_transport) * getter::element<0, 2>(dFdr) +
      getter::element<4, 6>(J_transport) * getter::element<0, 0>(dFdt) +
      getter::element<5, 6>(J_transport) * getter::element<0, 1>(dFdt) +
      getter::element<6, 6>(J_transport) * getter::element<0, 2>(dFdt);
  getter::element<1, 6>(new_J) =
      getter::element<0, 6>(J_transport) * getter::element<1, 0>(dFdr) +
      getter::element<1, 6>(J_transport) * getter::element<1, 1>(dFdr) +
      getter::element<2, 6>(J_transport) * getter::element<1, 2>(dFdr) +
      getter::element<4, 6>(J_transport) * getter::element<1, 0>(dFdt) +
      getter::element<5, 6>(J_transport) * getter::element<1, 1>(dFdt) +
      getter::element<6, 6>(J_transport) * getter::element<1, 2>(dFdt);
  getter::element<2, 6>(new_J) =
      getter::element<0, 6>(J_transport) * getter::element<2, 0>(dFdr) +
      getter::element<1, 6>(J_transport) * getter::element<2, 1>(dFdr) +
      getter::element<2, 6>(J_transport) * getter::element<2, 2>(dFdr) +
      getter::element<4, 6>(J_transport) * getter::element<2, 0>(dFdt) +
      getter::element<5, 6>(J_transport) * getter::element<2, 1>(dFdt) +
      getter::element<6, 6>(J_transport) * getter::element<2, 2>(dFdt);

  getter::element<4, 6>(new_J) =
      getter::element<0, 6>(J_transport) * getter::element<0, 0>(dGdr) +
      getter::element<1, 6>(J_transport) * getter::element<0, 1>(dGdr) +
      getter::element<2, 6>(J_transport) * getter::element<0, 2>(dGdr) +
      getter::element<4, 6>(J_transport) * getter::element<0, 0>(dGdt) +
      getter::element<5, 6>(J_transport) * getter::element<0, 1>(dGdt) +
      getter::element<6, 6>(J_transport) * getter::element<0, 2>(dGdt);
  getter::element<5, 6>(new_J) =
      getter::element<0, 6>(J_transport) * getter::element<1, 0>(dGdr) +
      getter::element<1, 6>(J_transport) * getter::element<1, 1>(dGdr) +
      getter::element<2, 6>(J_transport) * getter::element<1, 2>(dGdr) +
      getter::element<4, 6>(J_transport) * getter::element<1, 0>(dGdt) +
      getter::element<5, 6>(J_transport) * getter::element<1, 1>(dGdt) +
      getter::element<6, 6>(J_transport) * getter::element<1, 2>(dGdt);
  getter::element<6, 6>(new_J) =
      getter::element<0, 6>(J_transport) * getter::element<2, 0>(dGdr) +
      getter::element<1, 6>(J_transport) * getter::element<2, 1>(dGdr) +
      getter::element<2, 6>(J_transport) * getter::element<2, 2>(dGdr) +
      getter::element<4, 6>(J_transport) * getter::element<2, 0>(dGdt) +
      getter::element<5, 6>(J_transport) * getter::element<2, 1>(dGdt) +
      getter::element<6, 6>(J_transport) * getter::element<2, 2>(dGdt);

  getter::element<0, 7>(new_J) =
      getter::element<0, 7>(J_transport) * getter::element<0, 0>(dFdr) +
      getter::element<1, 7>(J_transport) * getter::element<0, 1>(dFdr) +
      getter::element<2, 7>(J_transport) * getter::element<0, 2>(dFdr) +
      getter::element<4, 7>(J_transport) * getter::element<0, 0>(dFdt) +
      getter::element<5, 7>(J_transport) * getter::element<0, 1>(dFdt) +
      getter::element<6, 7>(J_transport) * getter::element<0, 2>(dFdt) +
      getter::element<7, 7>(J_transport) * getter::element<0>(dFdqop);
  getter::element<1, 7>(new_J) =
      getter::element<0, 7>(J_transport) * getter::element<1, 0>(dFdr) +
      getter::element<1, 7>(J_transport) * getter::element<1, 1>(dFdr) +
      getter::element<2, 7>(J_transport) * getter::element<1, 2>(dFdr) +
      getter::element<4, 7>(J_transport) * getter::element<1, 0>(dFdt) +
      getter::element<5, 7>(J_transport) * getter::element<1, 1>(dFdt) +
      getter::element<6, 7>(J_transport) * getter::element<1, 2>(dFdt) +
      getter::element<7, 7>(J_transport) * getter::element<1>(dFdqop);
  getter::element<2, 7>(new_J) =
      getter::element<0, 7>(J_transport) * getter::element<2, 0>(dFdr) +
      getter::element<1, 7>(J_transport) * getter::element<2, 1>(dFdr) +
      getter::element<2, 7>(J_transport) * getter::element<2, 2>(dFdr) +
      getter::element<4, 7>(J_transport) * getter::element<2, 0>(dFdt) +
      getter::element<5, 7>(J_transport) * getter::element<2, 1>(dFdt) +
      getter::element<6, 7>(J_transport) * getter::element<2, 2>(dFdt) +
      getter::element<7, 7>(J_transport) * getter::element<2>(dFdqop);

  getter::element<4, 7>(new_J) =
      getter::element<0, 7>(J_transport) * getter::element<0, 0>(dGdr) +
      getter::element<1, 7>(J_transport) * getter::element<0, 1>(dGdr) +
      getter::element<2, 7>(J_transport) * getter::element<0, 2>(dGdr) +
      getter::element<4, 7>(J_transport) * getter::element<0, 0>(dGdt) +
      getter::element<5, 7>(J_transport) * getter::element<0, 1>(dGdt) +
      getter::element<6, 7>(J_transport) * getter::element<0, 2>(dGdt) +
      getter::element<7, 7>(J_transport) * getter::element<0>(dGdqop);
  getter::element<5, 7>(new_J) =
      getter::element<0, 7>(J_transport) * getter::element<1, 0>(dGdr) +
      getter::element<1, 7>(J_transport) * getter::element<1, 1>(dGdr) +
      getter::element<2, 7>(J_transport) * getter::element<1, 2>(dGdr) +
      getter::element<4, 7>(J_transport) * getter::element<1, 0>(dGdt) +
      getter::element<5, 7>(J_transport) * getter::element<1, 1>(dGdt) +
      getter::element<6, 7>(J_transport) * getter::element<1, 2>(dGdt) +
      getter::element<7, 7>(J_transport) * getter::element<1>(dGdqop);
  getter::element<6, 7>(new_J) =
      getter::element<0, 7>(J_transport) * getter::element<2, 0>(dGdr) +
      getter::element<1, 7>(J_transport) * getter::element<2, 1>(dGdr) +
      getter::element<2, 7>(J_transport) * getter::element<2, 2>(dGdr) +
      getter::element<4, 7>(J_transport) * getter::element<2, 0>(dGdt) +
      getter::element<5, 7>(J_transport) * getter::element<2, 1>(dGdt) +
      getter::element<6, 7>(J_transport) * getter::element<2, 2>(dGdt) +
      getter::element<7, 7>(J_transport) * getter::element<2>(dGdqop);
  getter::element<7, 7>(new_J) = dqopqop * getter::element<7, 7>(J_transport);
}
}  // namespace detray::detail
