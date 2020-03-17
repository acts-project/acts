// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// AlgebraHelper.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

#include <bitset>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <memory>
#include <string>
#include <vector>

#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/TypeTraits.hpp"
#include "Acts/Utilities/detail/DefaultParameterDefinitions.hpp"

#define ACTS_CHECK_BIT(value, mask) ((value & mask) == mask)

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
template <class T>
using phi_method_t = decltype(std::declval<const T>().phi());

template <class T>
using has_phi_method = concept ::is_detected<phi_method_t, T>;

}  // namespace detail

/// Calculate phi (transverse plane angle) from compatible Eigen types
/// @tparam Derived Eigen derived concrete type
/// @param v Any vector like Eigen type, static or dynamic
/// @note Will static assert that the number of rows of @p v is at least 2, or
/// in case of dynamic size, will abort execution if that is not the case.
/// @return The value of the angle in the transverse plane.
template <typename Derived>
double phi(const Eigen::MatrixBase<Derived>& v) noexcept {
  constexpr int rows = Eigen::MatrixBase<Derived>::RowsAtCompileTime;
  if constexpr (rows != -1) {
    // static size, do compile time check
    static_assert(rows >= 2,
                  "Phi function not valid for vectors not at least 2D");
  } else {
    // dynamic size
    if (v.rows() < 2) {
      std::cerr << "Phi function not valid for vectors not at least 2D"
                << std::endl;
      std::abort();
    }
  }

  return std::atan2(v[1], v[0]);
}

/// Calculate phi (transverse plane angle) from anything implementing a method
/// like `phi()` returing anything convertible to `double`.
/// @tparam T anything that has a phi method
/// @param v Any type that implements a phi method
/// @return The phi value
template <typename T,
          std::enable_if_t<detail::has_phi_method<T>::value, int> = 0>
double phi(const T& v) noexcept {
  return v.phi();
}

/// Calculate radius in the transverse (xy) plane of a vector
/// @tparam Derived Eigen derived concrete type
/// @param v Any vector like Eigen type, static or dynamic
/// @note Will static assert that the number of rows of @p v is at least 2, or
/// in case of dynamic size, will abort execution if that is not the case.
/// @return The transverse radius value.
template <typename Derived>
double perp(const Eigen::MatrixBase<Derived>& v) noexcept {
  constexpr int rows = Eigen::MatrixBase<Derived>::RowsAtCompileTime;
  if constexpr (rows != -1) {
    // static size, do compile time check
    static_assert(rows >= 2,
                  "Perp function not valid for vectors not at least 2D");
  } else {
    // dynamic size
    if (v.rows() < 2) {
      std::cerr << "Perp function not valid for vectors not at least 2D"
                << std::endl;
      std::abort();
    }
  }
  return std::sqrt(v[0] * v[0] + v[1] * v[1]);
}

/// Calculate the theta angle (longitudinal w.r.t. z axis) of a vector
/// @tparam Derived Eigen derived concrete type
/// @param v Any vector like Eigen type, static or dynamic
/// @note Will static assert that the number of rows of @p v is at least 3, or
/// in case of dynamic size, will abort execution if that is not the case.
/// @return The theta value
template <typename Derived>
double theta(const Eigen::MatrixBase<Derived>& v) noexcept {
  constexpr int rows = Eigen::MatrixBase<Derived>::RowsAtCompileTime;
  if constexpr (rows != -1) {
    // static size, do compile time check
    static_assert(rows >= 3, "Theta function not valid for non-3D vectors.");
  } else {
    // dynamic size
    if (v.rows() < 3) {
      std::cerr << "Theta function not valid for non-3D vectors." << std::endl;
      std::abort();
    }
  }

  return std::atan2(std::sqrt(v[0] * v[0] + v[1] * v[1]), v[2]);
}

/// Calculate the pseudorapidity for a vector.
/// @tparam Derived Eigen derived concrete type
/// @param v Any vector like Eigen type, static or dynamic
/// @note Will static assert that the number of rows of @p v is at least 3, or
/// in case of dynamic size, will abort execution if that is not the case.
/// @return The pseudorapidity value
template <typename Derived>
double eta(const Eigen::MatrixBase<Derived>& v) noexcept {
  constexpr int rows = Eigen::MatrixBase<Derived>::RowsAtCompileTime;
  if constexpr (rows != -1) {
    // static size, do compile time check
    static_assert(rows >= 3, "Eta function not valid for non-3D vectors.");
  } else {
    // dynamic size
    if (v.rows() < 3) {
      std::cerr << "Eta function not valid for non-3D vectors." << std::endl;
      std::abort();
    }
  }

  return std::atanh(v[2] / v.norm());
}

/// Helper method to cast out the binning value from a 3D Vector
///
/// For this method a 3D vector is required to guarantee all potential
/// binning values
///
static double cast(const Vector3D& position, BinningValue bval) {
  if (bval < 3)
    return position[bval];
  switch (bval) {
    case binR:
      return perp(position);
      break;
    case binPhi:
      return phi(position);
      break;
    case binH:
      return theta(position);
      break;
    case binEta:
      return eta(position);
      break;
    case binMag:
      return position.norm();
      break;
  }
  return 0.;
}

/// @brief Calculates column-wise cross products of a matrix and a vector and
/// stores the result column-wise in a matrix.
///
/// @param [in] m Matrix that will be used for cross products
/// @param [in] v Vector for cross products
/// @return Constructed matrix
inline ActsMatrixD<3, 3> cross(const ActsMatrixD<3, 3>& m, const Vector3D& v) {
  ActsMatrixD<3, 3> r;
  r.col(0) = m.col(0).cross(v);
  r.col(1) = m.col(1).cross(v);
  r.col(2) = m.col(2).cross(v);

  return r;
}

/// @brief Access to the time component of input parameter
///
/// @param spacePointVec The SpacePointVector
/// @return Reference to the time component
inline ParValue_t& time(SpacePointVector& spacePointVec) {
  return spacePointVec[3];
}

/// @brief Const overload access to the
/// time component of input parameter
///
/// @param spacePointVec The SpacePointVector
/// @return Reference to the time component
inline const ParValue_t& time(const SpacePointVector& spacePointVec) {
  return spacePointVec[3];
}

/// @brief Access to the time component of input parameter
///
/// @param boundVec The BoundVector
/// @return Reference to the time component
inline ParValue_t& time(BoundVector& boundVec) {
  return boundVec[eT];
}

/// @brief Const overload access to the
/// time component of input parameter
///
/// @param boundVec The BoundVector
/// @return Reference to the time component
inline const ParValue_t& time(const BoundVector& boundVec) {
  return boundVec[eT];
}

/// @brief Access to the time component of input parameter
///
/// @param freeVec The FreeVector
/// @return Reference to the time component
inline ParValue_t& time(FreeVector& freeVec) {
  return freeVec[7];
}

/// @brief Const overload access to the
/// time component of input parameter
///
/// @param freeVec The FreeVector
/// @return Reference to the time component
inline const ParValue_t& time(const FreeVector& freeVec) {
  return freeVec[7];
}

/// @brief Access to the position components of input parameter
///
/// @param spacePointVec The SpacePointVector
/// @return Reference to the position components
inline auto position(SpacePointVector& spacePointVec) {
  return spacePointVec.head<3>();
}

/// @brief Const overload access to the
/// position components of input parameter
///
/// @param spacePointVec The SpacePointVector
/// @return Reference to the position components
inline auto position(const SpacePointVector& spacePointVec) {
  return spacePointVec.head<3>();
}

/// @brief Access to the
/// position components of input parameter
///
/// @param freeVec The SpacePointVector
/// @return Reference to the position components
inline auto position(FreeVector& freeVec) {
  return freeVec.head<3>();
}

/// @brief Const overload access to the
/// position components of input parameter
///
/// @param freeVec The SpacePointVector
/// @return Reference to the position components
inline auto position(const FreeVector& freeVec) {
  return freeVec.head<3>();
}

}  // namespace VectorHelpers

namespace detail {

inline double roundWithPrecision(double val, int precision) {
  if (val < 0 && std::abs(val) * std::pow(10, precision) < 1.) {
    return -val;
  }
  return val;
}
}  // namespace detail

/// Print out a matrix in a structured way.
/// @param matrix The matrix to print
/// @param precision Numeric output precision
/// @param offset Offset in front of matrix lines
/// @return The printed string
inline std::string toString(const ActsMatrixXd& matrix, int precision = 4,
                            const std::string& offset = "") {
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
      if (i != matrix.rows() -
                   1) {  // make the end line and the offset in the next line
        sout << std::endl;
        sout << offset;
      }
    }
  }
  return sout.str();
}

/// Print out a translation in a structured way.
/// @param matrix The translation to print
/// @param precision Numeric output precision
/// @return The printed string
inline std::string toString(const Acts::Translation3D& translation,
                            int precision = 4) {
  Acts::Vector3D trans;
  trans[0] = translation.x();
  trans[1] = translation.y();
  trans[2] = translation.z();
  return toString(trans, precision);
}

/// Print out a transform in a structured way.
/// @param matrix The transform to print
/// @param precision Numeric output precision
/// @param offset Offset in front of matrix lines
/// @return The printed string
inline std::string toString(const Acts::Transform3D& transform,
                            int precision = 4, const std::string& offset = "") {
  std::ostringstream sout;
  sout << "Translation : " << toString(transform.translation(), precision)
       << std::endl;
  std::string rotationOffset = offset + "              ";
  sout << offset << "Rotation    : "
       << toString(transform.rotation(), precision + 2, rotationOffset);
  return sout.str();
}

/// Helper function to unpack a vector of @c shared_ptr into a vector of raw
/// pointers
/// @tparam T the stored type
/// @param items The vector of @c shared_ptr
/// @return The unpacked vector
template <typename T>
std::vector<T*> unpack_shared_vector(
    const std::vector<std::shared_ptr<T>>& items) {
  std::vector<T*> rawPtrs;
  rawPtrs.reserve(items.size());
  for (const std::shared_ptr<T>& item : items) {
    rawPtrs.push_back(item.get());
  }
  return rawPtrs;
}

/// Helper function to unpack a vector of @c shared_ptr into a vector of raw
/// pointers (const version)
/// @tparam T the stored type
/// @param items The vector of @c shared_ptr
/// @return The unpacked vector
template <typename T>
std::vector<const T*> unpack_shared_vector(
    const std::vector<std::shared_ptr<const T>>& items) {
  std::vector<const T*> rawPtrs;
  rawPtrs.reserve(items.size());
  for (const std::shared_ptr<const T>& item : items) {
    rawPtrs.push_back(item.get());
  }
  return rawPtrs;
}

/// @brief Dispatch a call based on a runtime value on a function taking the
/// value at compile time.
///
/// This function allows to write a templated functor, which accepts a @c size_t
/// like paramater at compile time. It is then possible to make a call to the
/// corresponding instance of the functor based on a runtime value. To achieve
/// this, the function essentially created a if cascade between @c N and @c
/// NMAX, attempting to find the right instance. Because the cascade is visible
/// to the compiler entirely, it should be able to optimize.
///
/// @tparam Callable Type which takes a size_t as a compile time param
/// @tparam N Value from which to start the dispatch chain, i.e. 0 in most cases
/// @tparam NMAX Maximum value up to which to attempt a dispatch
/// @param v The runtime value to dispatch on
/// @param args Additional arguments passed to @c Callable::invoke().
/// @note @c Callable is expected to have a static member function @c invoke
/// that is callable with @c Args
template <template <size_t> class Callable, size_t N, size_t NMAX,
          typename... Args>
decltype(Callable<N>::invoke(std::declval<Args>()...)) template_switch(
    size_t v, Args&&... args) {
  if (v == N) {
    return Callable<N>::invoke(std::forward<Args>(args)...);
  }
  if constexpr (N < NMAX) {
    return template_switch<Callable, N + 1, NMAX>(v,
                                                  std::forward<Args>(args)...);
  }
  std::cerr << "template_switch<Fn, " << N << ", " << NMAX << ">(v=" << v
            << ") is not valid (v > NMAX)" << std::endl;
  std::abort();
}

/// Convert a bitset to a matrix of integers, with each element set to the bit
/// value.
/// @note How the bits are assigned to matrix elements depends on the storage
/// type of the matrix being converted (row-major or col-major)
/// @tparam MatrixType Matrix type that is produced
/// @param bs The bitset to convert
/// @return A matrix with the integer values of the bits from @p bs
template <typename MatrixType>
MatrixType bitsetToMatrix(const std::bitset<MatrixType::RowsAtCompileTime *
                                            MatrixType::ColsAtCompileTime>
                              bs) {
  constexpr int rows = MatrixType::RowsAtCompileTime;
  constexpr int cols = MatrixType::ColsAtCompileTime;

  static_assert(rows != -1 && cols != -1,
                "bitsetToMatrix does not support dynamic matrices");

  MatrixType m;
  auto* p = m.data();
  for (size_t i = 0; i < rows * cols; i++) {
    p[i] = bs[rows * cols - 1 - i];
  }
  return m;
}

/// Convert an integer matrix to a bitset.
/// @note How the bits are ordered depends on the storage type of the matrix
/// being converted (row-major or col-major)
/// @tparam Derived Eigen base concrete type
/// @param m Matrix that is converted
/// @return The converted bitset.
template <typename Derived>
auto matrixToBitset(const Eigen::PlainObjectBase<Derived>& m) {
  using MatrixType = Eigen::PlainObjectBase<Derived>;
  constexpr size_t rows = MatrixType::RowsAtCompileTime;
  constexpr size_t cols = MatrixType::ColsAtCompileTime;

  std::bitset<rows * cols> res;

  auto* p = m.data();
  for (size_t i = 0; i < rows * cols; i++) {
    res[rows * cols - 1 - i] = p[i];
  }

  return res;
}

}  // namespace Acts
