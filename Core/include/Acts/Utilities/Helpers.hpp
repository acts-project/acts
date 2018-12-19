// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// AlgebraHelper.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

// libc/STL include(s)
#include <cmath>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>
#include "strings.h"

// Acts include(s)
#include "Definitions.hpp"

#include <boost/tti/has_member_function.hpp>

#ifndef ACTS_BIT_CODING
#define ACTS_BIT_CODING 1
#if (__GNUC__ >= 4)
#define ACTS_BIT_SHIFT(mask) __builtin_ctzl(mask)
#else
#define ACTS_BIT_SHIFT(mask) (ffsl(mask) - 1)
#endif
#define ACTS_BIT_ENCODE(value, mask) (value << ACTS_BIT_SHIFT(mask))
#define ACTS_BIT_DECODE(code, mask) ((code & mask) >> ACTS_BIT_SHIFT(mask))
#endif

/** Geometry primitives helper functions
 */

namespace Acts {

/** EventPrimitvesToStringConverter

    inline methods for conversion of EventPrimitives (Matrix)
    to std::string.

    This is to enhance formatted screen ouput and for ASCII based
    testing.

    The offset can be used to offset the lines (starting from line 2) wrt to the
    zero position for formatting reasons.


 */

namespace VectorHelpers {

  namespace detail {
    // helper to figure out if a type has a member called phi
    BOOST_TTI_HAS_MEMBER_FUNCTION(phi)
    template <typename T>
    using has_phi_method
        = has_member_function_phi<T,
                                  double,
                                  boost::mpl::vector<>,
                                  boost::function_types::const_qualified>;
  }

  // default call on Eigen types, calculate radius
  template <typename Derived>
  double
  phi(const Eigen::MatrixBase<Derived>& v)
  {
    if (v.rows() < 2) {
      return 0.;
    }
    return std::atan2(v[1], v[0]);
  }

  // if called-upon type has phi method, call that
  template <typename T,
            std::enable_if_t<detail::has_phi_method<T>::value, int> = 0>
  double
  phi(const T& v)
  {
    return v.phi();
  }

  template <typename Derived>
  double
  perp(const Eigen::MatrixBase<Derived>& v)
  {
    if (v.rows() < 2) {
      return 0.;
    }
    return std::sqrt(v[0] * v[0] + v[1] * v[1]);
  }

  template <typename Derived>
  double
  theta(const Eigen::MatrixBase<Derived>& v)
  {
    if (v.rows() < 3) {
      return 0.;
    }
    return std::atan2(std::sqrt(v[0] * v[0] + v[1] * v[1]), v[2]);
  }

  template <typename Derived>
  double
  eta(const Eigen::MatrixBase<Derived>& v)
  {
    return std::atanh(v[2] / v.norm());
  }

  /// @brief Calculates column-wise cross products of a matrix and a vector and
  /// stores the result column-wise in a matrix.
  ///
  /// @param [in] m Matrix that will be used for cross products
  /// @param [in] v Vector for cross products
  /// @return Constructed matrix
  inline ActsMatrixD<3, 3>
  cross(const ActsMatrixD<3, 3>& m, const Vector3D& v)
  {
    ActsMatrixD<3, 3> r;
    r.col(0) = m.col(0).cross(v);
    r.col(1) = m.col(1).cross(v);
    r.col(2) = m.col(2).cross(v);

    return r;
  }
}

namespace detail {

  inline double
  roundWithPrecision(double val, int precision)
  {
    if (val < 0 && std::abs(val) * std::pow(10, precision) < 1.) {
      return -val;
    }
    return val;
  }
}

inline std::string
toString(const ActsMatrixXd& matrix,
         int                 precision = 4,
         const std::string&  offset    = "")
{
  std::ostringstream sout;

  sout << std::setiosflags(std::ios::fixed) << std::setprecision(precision);
  if (matrix.cols() == 1) {
    sout << "(";
    for (int i = 0; i < matrix.rows(); ++i) {
      double val = detail::roundWithPrecision(matrix(i, 0), precision);
      sout << val;
      if (i != matrix.rows() - 1) {
        sout << ", ";
      }
    }
    sout << ")";
  } else {
    for (int i = 0; i < matrix.rows(); ++i) {
      for (int j = 0; j < matrix.cols(); ++j) {
        if (j == 0) {
          sout << "(";
        }
        double val = detail::roundWithPrecision(matrix(i, j), precision);
        sout << val;
        if (j == matrix.cols() - 1) {
          sout << ")";
        } else {
          sout << ", ";
        }
      }
      if (i
          != matrix.rows()
              - 1) {  // make the end line and the offset in the next line
        sout << std::endl;
        sout << offset;
      }
    }
  }
  return sout.str();
}

inline std::string
toString(const Acts::Translation3D& translation, int precision = 4)
{
  Acts::Vector3D trans;
  trans[0] = translation.x();
  trans[1] = translation.y();
  trans[2] = translation.z();
  return toString(trans, precision);
}

inline std::string
toString(const Acts::Transform3D& transform,
         int                      precision = 4,
         const std::string&       offset    = "")
{
  std::ostringstream sout;
  sout << "Translation : " << toString(transform.translation(), precision)
       << std::endl;
  std::string rotationOffset = offset + "              ";
  sout << offset << "Rotation    : "
       << toString(transform.rotation(), precision + 2, rotationOffset);
  return sout.str();
}

template <typename T>
std::vector<T*>
unpack_shared_vector(const std::vector<std::shared_ptr<T>>& items)
{
  std::vector<T*> rawPtrs;
  rawPtrs.reserve(items.size());
  for (const std::shared_ptr<T>& item : items) {
    rawPtrs.push_back(item.get());
  }
  return rawPtrs;
}

template <typename T>
std::vector<const T*>
unpack_shared_vector(const std::vector<std::shared_ptr<const T>>& items)
{
  std::vector<const T*> rawPtrs;
  rawPtrs.reserve(items.size());
  for (const std::shared_ptr<const T>& item : items) {
    rawPtrs.push_back(item.get());
  }
  return rawPtrs;
}

}  // end of Acts namespace
