// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once
// Acts include(s)
#include "Acts/Utilities/Definitions.hpp"

namespace Acts {
/// @cond detail
namespace detail {
  /// @brief  initialize projection matrices
  ///
  /// @tparam columns number of columns (= dimension of full parameter space)
  /// @tparam rows template parameter pack containing the indices of the
  ///         parameters to be projected
  ///
  /// This struct provides an initialization method for a projection matrix M
  /// such that only the entries with the given indices are selected from a full
  /// parameter vector. That means, M is a mapping M: (Nx1) --> (Sx1) if N is
  /// the total number of parameters and S is the number of given indices.
  ///
  /// @return make_projection_matrix<columns,rows...>::init() returns a matrix
  ///         with dimensions (`sizeof...(rows)` x columns)
  template <unsigned int columns, unsigned int... rows>
  struct make_projection_matrix;

  /// @cond
  // build projection matrix by iteratively stacking row vectors
  template <unsigned int columns, unsigned int i, unsigned int... N>
  struct make_projection_matrix<columns, i, N...>
  {
    static ActsMatrixD<sizeof...(N) + 1, columns>
    init()
    {
      ActsRowVectorD<columns> v;
      v.setZero();
      v(i) = 1;

      ActsMatrixD<sizeof...(N) + 1, columns> m;
      m.row(0) << v;
      m.block(1, 0, sizeof...(N), columns)
          << make_projection_matrix<columns, N...>::init();

      return m;
    }
  };

  // projection matrix for a single local parameter is a simple row vector
  template <unsigned int columns, unsigned int i>
  struct make_projection_matrix<columns, i>
  {
    static ActsRowVectorD<columns>
    init()
    {
      ActsRowVectorD<columns> v;
      v.setZero();
      v(i) = 1;
      return v;
    }
  };
  /// @endcond
}  // namespace detail
/// @endcond
}  // namespace Acts