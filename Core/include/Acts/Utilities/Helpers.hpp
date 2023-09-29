// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <algorithm>
#include <iostream>
#include <limits>
#include <memory>
#include <optional>
#include <string>
#include <vector>

#define ACTS_CHECK_BIT(value, mask) ((value & mask) == mask)

namespace Acts {

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

/// This can be abandoned with C++20 to use the std::to_array method
///
/// @note only the first kDIM elements will obviously be filled, if the
/// vector tends to be longer, it is truncated
///
/// @param vecvals the vector of bound values to be converted
/// @return an array with the filled values
template <std::size_t kDIM, typename value_type>
std::array<value_type, kDIM> to_array(const std::vector<value_type>& vecvals) {
  std::array<value_type, kDIM> rarray = {};
  for (const auto [iv, v] : enumerate(vecvals)) {
    if (iv < kDIM) {
      rarray[iv] = v;
    }
  }
  return rarray;
}

/// @brief Dispatch a call based on a runtime value on a function taking the
/// value at compile time.
///
/// This function allows to write a templated functor, which accepts a @c size_t
/// like parameter at compile time. It is then possible to make a call to the
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
