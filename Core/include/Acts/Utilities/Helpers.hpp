// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Common.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Utilities/BinningType.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

#include <array>
#include <bitset>
#include <cassert>
#include <cmath>
#include <cstdlib>
#include <iomanip>
#include <iostream>
#include <limits>
#include <memory>
#include <string>
#include <tuple>
#include <type_traits>
#include <vector>

#define ACTS_CHECK_BIT(value, mask) ((value & mask) == mask)

namespace Acts {

namespace VectorHelpers {
namespace detail {
template <class T>
using phi_method_t = decltype(std::declval<const T>().phi());

template <class T>
using has_phi_method = Concepts ::is_detected<phi_method_t, T>;

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
    static_assert(rows == 3, "Theta function not valid for non-3D vectors.");
  } else {
    // dynamic size
    if (v.rows() != 3) {
      std::cerr << "Theta function not valid for non-3D vectors." << std::endl;
      std::abort();
    }
  }

  return std::atan2(std::sqrt(v[0] * v[0] + v[1] * v[1]), v[2]);
}

/// @brief Fast evaluation of trigonomic functions.
///
/// @param direction for this evaluatoin
///
/// @return cos(phi), sin(phi), cos(theta), sin(theta), 1/sin(theta)
static inline const std::array<ActsScalar, 5> evaluateTrigonomics(
    const Vector3& direction) {
  const ActsScalar x = direction(0);  // == cos(phi) * sin(theta)
  const ActsScalar y = direction(1);  // == sin(phi) * sin(theta)
  const ActsScalar z = direction(2);  // == cos(theta)
  // can be turned into cosine/sine
  const ActsScalar cosTheta = z;
  const ActsScalar sinTheta = std::sqrt(x * x + y * y);
  const ActsScalar invSinTheta = 1. / sinTheta;
  const ActsScalar cosPhi = x * invSinTheta;
  const ActsScalar sinPhi = y * invSinTheta;

  return {cosPhi, sinPhi, cosTheta, sinTheta, invSinTheta};
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
    static_assert(rows == 3, "Eta function not valid for non-3D vectors.");
  } else {
    // dynamic size
    if (v.rows() != 3) {
      std::cerr << "Eta function not valid for non-3D vectors." << std::endl;
      std::abort();
    }
  }

  return std::atanh(v[2] / v.norm());
}

/// Helper method to extract the binning value from a 3D vector.
///
/// For this method a 3D vector is required to guarantee all potential
/// binning values.
inline double cast(const Vector3& position, BinningValue bval) {
  switch (bval) {
    case binX:
      return position[0];
    case binY:
      return position[1];
    case binZ:
      return position[2];
    case binR:
      return perp(position);
    case binPhi:
      return phi(position);
    case binRPhi:
      return perp(position) * phi(position);
    case binH:
      return theta(position);
    case binEta:
      return eta(position);
    case binMag:
      return position.norm();
    default:
      assert(false and "Invalid BinningValue enum value");
      return std::numeric_limits<double>::quiet_NaN();
  }
}

/// @brief Calculates column-wise cross products of a matrix and a vector and
/// stores the result column-wise in a matrix.
///
/// @param [in] m Matrix that will be used for cross products
/// @param [in] v Vector for cross products
/// @return Constructed matrix
inline ActsMatrix<3, 3> cross(const ActsMatrix<3, 3>& m, const Vector3& v) {
  ActsMatrix<3, 3> r;
  r.col(0) = m.col(0).cross(v);
  r.col(1) = m.col(1).cross(v);
  r.col(2) = m.col(2).cross(v);

  return r;
}

/// Access the three-position components in a four-position vector.
inline auto position(const Vector4& pos4) {
  return pos4.segment<3>(ePos0);
}

/// Access the three-position components in a free parameters vector.
inline auto position(const FreeVector& params) {
  return params.segment<3>(eFreePos0);
}

/// Construct a four-vector from a three-vector and scalar fourth component.
template <typename vector3_t>
inline auto makeVector4(const Eigen::MatrixBase<vector3_t>& vec3,
                        typename vector3_t::Scalar w)
    -> Eigen::Matrix<typename vector3_t::Scalar, 4, 1> {
  EIGEN_STATIC_ASSERT_VECTOR_SPECIFIC_SIZE(vector3_t, 3);

  Eigen::Matrix<typename vector3_t::Scalar, 4, 1> vec4;
  vec4[ePos0] = vec3[ePos0];
  vec4[ePos1] = vec3[ePos1];
  vec4[ePos2] = vec3[ePos2];
  vec4[eTime] = w;
  return vec4;
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
///
/// @tparam derived_t Type of the matrix
/// @param matrix The matrix to print
/// @param precision Numeric output precision
/// @param offset Offset in front of matrix lines
/// @return The printed string
template <typename derived_t>
inline std::string toString(const Eigen::MatrixBase<derived_t>& matrix,
                            int precision = 4, const std::string& offset = "") {
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
/// @param translation The translation to print
/// @param precision Numeric output precision
/// @return The printed string
inline std::string toString(const Acts::Translation3& translation,
                            int precision = 4) {
  Acts::Vector3 trans;
  trans[0] = translation.x();
  trans[1] = translation.y();
  trans[2] = translation.z();
  return toString(trans, precision);
}

/// Print out a transform in a structured way.
/// @param transform The transform to print
/// @param precision Numeric output precision
/// @param offset Offset in front of matrix lines
/// @return The printed string
inline std::string toString(const Acts::Transform3& transform,
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

/// Helper function to unpack a vector of @c shared_ptr into a vector of raw
/// pointers
/// @tparam T the stored type
/// @param items The vector of @c shared_ptr
/// @return The unpacked vector
template <typename T>
std::vector<const T*> unpack_shared_const_vector(
    const std::vector<std::shared_ptr<T>>& items) {
  std::vector<const T*> rawPtrs;
  rawPtrs.reserve(items.size());
  for (const std::shared_ptr<T>& item : items) {
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
auto template_switch(size_t v, Args&&... args) {
  if (v == N) {
    return Callable<N>::invoke(std::forward<Args>(args)...);
  }
  if (v == 0) {
    std::cerr << "template_switch<Fn, " << N << ", " << NMAX << ">(v=" << v
              << ") is not valid (v == 0 and N != 0)" << std::endl;
    std::abort();
  }
  if constexpr (N < NMAX) {
    return template_switch<Callable, N + 1, NMAX>(v,
                                                  std::forward<Args>(args)...);
  }
  std::cerr << "template_switch<Fn, " << N << ", " << NMAX << ">(v=" << v
            << ") is not valid (v > NMAX)" << std::endl;
  std::abort();
}

/// Alternative version of @c template_switch which accepts a generic
/// lambda and communicates the dimension via an integral constant type
/// @tparam N Value from which to start the dispatch chain, i.e. 0 in most cases
/// @tparam NMAX Maximum value up to which to attempt a dispatch
/// @param v The runtime value to dispatch on
/// @param func The lambda to invoke
/// @param args Additional arguments passed to @p func
template <size_t N, size_t NMAX, typename Lambda, typename... Args>
auto template_switch_lambda(size_t v, Lambda&& func, Args&&... args) {
  if (v == N) {
    return func(std::integral_constant<size_t, N>{},
                std::forward<Args>(args)...);
  }
  if (v == 0) {
    std::cerr << "template_switch<Fn, " << N << ", " << NMAX << ">(v=" << v
              << ") is not valid (v == 0 and N != 0)" << std::endl;
    std::abort();
  }
  if constexpr (N < NMAX) {
    return template_switch_lambda<N + 1, NMAX>(v, func,
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
    res[rows * cols - 1 - i] = static_cast<bool>(p[i]);
  }

  return res;
}

/// @brief Perform a blocked matrix multiplication, avoiding Eigen GEMM
/// methods.
///
/// @tparam A The type of the first matrix, which should be MxN
/// @tparam B The type of the second matrix, which should be NxP
///
/// @param[in] a An MxN matrix of type A
/// @param[in] b An NxP matrix of type P
///
/// @returns The product ab
template <typename A, typename B>
inline ActsMatrix<A::RowsAtCompileTime, B::ColsAtCompileTime> blockedMult(
    const A& a, const B& b) {
  // Extract the sizes of the matrix types that we receive as template
  // parameters.
  constexpr int M = A::RowsAtCompileTime;
  constexpr int N = A::ColsAtCompileTime;
  constexpr int P = B::ColsAtCompileTime;

  // Ensure that the second dimension of our first matrix equals the first
  // dimension of the second matrix, otherwise we cannot multiply.
  static_assert(N == B::RowsAtCompileTime);

  if constexpr (M <= 4 && N <= 4 && P <= 4) {
    // In cases where the matrices being multiplied are small, we can rely on
    // Eigen do to a good job, and we don't really need to do any blocking.
    return a * b;
  } else {
    // Here, we want to calculate the expression: C = AB, Eigen, natively,
    // doesn't do a great job at this if the matrices A and B are large
    // (roughly M >= 8, N >= 8, or P >= 8), and applies a slow GEMM operation.
    // We apply a blocked matrix multiplication operation to decompose the
    // multiplication into smaller operations, relying on the fact that:
    //
    // ┌         ┐   ┌         ┐ ┌         ┐
    // │ C₁₁ C₁₂ │ = │ A₁₁ A₁₂ │ │ B₁₁ B₁₂ │
    // │ C₂₁ C₂₂ │ = │ A₂₁ A₂₂ │ │ B₂₁ B₂₂ │
    // └         ┘   └         ┘ └         ┘
    //
    // where:
    //
    // C₁₁ = A₁₁ * B₁₁ + A₁₂ * B₂₁
    // C₁₂ = A₁₁ * B₁₂ + A₁₂ * B₂₂
    // C₂₁ = A₂₁ * B₁₁ + A₂₂ * B₂₁
    // C₂₂ = A₂₁ * B₁₂ + A₂₂ * B₂₂
    //
    // The sizes of these submatrices are roughly half (in each dimension) that
    // of the parent matrix. If the size of the parent matrix is even, we can
    // divide it exactly, If the size of the parent matrix is odd, then some
    // of the submatrices will be one larger than the others. In general, for
    // any matrix Q, the sizes of the submatrices are (where / denotes integer
    // division):
    //
    // Q₁₁ : M / 2 × P / 2
    // Q₁₂ : M / 2 × (P + 1) / 2
    // Q₂₁ : (M + 1) / 2 × P / 2
    // Q₂₂ : (M + 1) / 2 × (P + 1) / 2
    //
    // See https://csapp.cs.cmu.edu/public/waside/waside-blocking.pdf for a
    // more in-depth explanation of blocked matrix multiplication.
    constexpr int M1 = M / 2;
    constexpr int M2 = (M + 1) / 2;
    constexpr int N1 = N / 2;
    constexpr int N2 = (N + 1) / 2;
    constexpr int P1 = P / 2;
    constexpr int P2 = (P + 1) / 2;

    // Construct the end result in this matrix, which destroys a few of Eigen's
    // built-in optimization techniques, but sadly this is necessary.
    ActsMatrix<M, P> r;

    // C₁₁ = A₁₁ * B₁₁ + A₁₂ * B₂₁
    r.template topLeftCorner<M1, P1>().noalias() =
        a.template topLeftCorner<M1, N1>() *
            b.template topLeftCorner<N1, P1>() +
        a.template topRightCorner<M1, N2>() *
            b.template bottomLeftCorner<N2, P1>();

    // C₁₂ = A₁₁ * B₁₂ + A₁₂ * B₂₂
    r.template topRightCorner<M1, P2>().noalias() =
        a.template topLeftCorner<M1, N1>() *
            b.template topRightCorner<N1, P2>() +
        a.template topRightCorner<M1, N2>() *
            b.template bottomRightCorner<N2, P2>();

    // C₂₁ = A₂₁ * B₁₁ + A₂₂ * B₂₁
    r.template bottomLeftCorner<M2, P1>().noalias() =
        a.template bottomLeftCorner<M2, N1>() *
            b.template topLeftCorner<N1, P1>() +
        a.template bottomRightCorner<M2, N2>() *
            b.template bottomLeftCorner<N2, P1>();

    // C₂₂ = A₂₁ * B₁₂ + A₂₂ * B₂₂
    r.template bottomRightCorner<M2, P2>().noalias() =
        a.template bottomLeftCorner<M2, N1>() *
            b.template topRightCorner<N1, P2>() +
        a.template bottomRightCorner<M2, N2>() *
            b.template bottomRightCorner<N2, P2>();

    return r;
  }
}

/// Clamp a numeric value to another type, respecting range of the target type
/// @tparam T the target type
/// @tparam U the source type
/// @param value the value to clamp
/// @return the clamped value
template <typename T, typename U>
T clampValue(U value) {
  return std::clamp(value, static_cast<U>(std::numeric_limits<T>::lowest()),
                    static_cast<U>(std::numeric_limits<T>::max()));
}

/// Return min/max from a (optionally) sorted series, obsolete with C++20
/// (ranges)
///
/// @tparam T a numeric series
///
/// @param tseries is the number series
///
/// @return [ min, max ] in an array of length 2
template <typename T>
std::array<typename T::value_type, 2u> min_max(const T& tseries) {
  return {*std::min_element(tseries.begin(), tseries.end()),
          *std::max_element(tseries.begin(), tseries.end())};
}

/// Return range and medium of a sorted numeric series
///
/// @tparam T a numeric series
///
/// @param tseries is the number series
///
/// @return [ range, medium ] in an tuple
template <typename T>
std::tuple<typename T::value_type, ActsScalar> range_medium(const T& tseries) {
  auto [min, max] = min_max(tseries);
  typename T::value_type range = (max - min);
  ActsScalar medium = static_cast<ActsScalar>((max + min) * 0.5);
  return std::tie(range, medium);
}

}  // namespace Acts
