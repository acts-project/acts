// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#ifndef ACTS_SEEDING_PERIODIC_HPP
#define ACTS_SEEDING_PERIODIC_HPP

#include <cmath>

namespace Acts {
namespace Seeding {
  namespace detail {

    /** @brief Wrap a periodic value back into the nominal range. */
    template <typename T>
    inline T
    wrap_periodic(T value, T start, T range)
    {
      using std::floor;
      return value - range * floor((value - start) / range);
    }

    /** @brief Calculate the equivalent angle in the [0, 2*pi) range */
    template <typename T>
    inline T
    radian_pos(T x)
    {
      return wrap_periodic<T>(x, T(0), T(2 * M_PI));
    }

    /** @brief Calculate the equivalent angle in the [-pi, pi] range */
    template <typename T>
    inline T
    radian_sym(T x)
    {
      return wrap_periodic<T>(x, T(-M_PI), T(2 * M_PI));
    }

  }  // namespace detail
}  // namespace Seeding
}  // namespace Acts

#endif  // ACTS_SEEDING_PERIODIC_HPP
