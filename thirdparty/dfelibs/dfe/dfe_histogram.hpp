// SPDX-License-Identifier: MIT
// Copyright 2018 Moritz Kiehn
//
// Permission is hereby granted, free of charge, to any person obtaining a copy
// of this software and associated documentation files (the "Software"), to deal
// in the Software without restriction, including without limitation the rights
// to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
// copies of the Software, and to permit persons to whom the Software is
// furnished to do so, subject to the following conditions:
//
// The above copyright notice and this permission notice shall be included in
// all copies or substantial portions of the Software.
//
// THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
// IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
// FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
// AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
// LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
// OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
// SOFTWARE.

/// \file
/// \brief   N-dimensional histograms with composable axes
/// \author  Moritz Kiehn <msmk@cern.ch>
/// \date    2018-05-22

#pragma once

#include <algorithm>
#include <array>
#include <functional>
#include <initializer_list>
#include <numeric>
#include <stdexcept>
#include <type_traits>
#include <utility>
#include <vector>

namespace dfe {
namespace histogram_impl {
namespace {

/// A simple n-dimensional array.
///
/// Data must be accessed through n-dimensional indices. The internal
/// storage format is considered an implementation detail. The size along
/// each dimension is set at run-time.
template<typename T, std::size_t NDimensions>
class NArray {
public:
  using Index = std::array<std::size_t, NDimensions>;

  /// Construct default-initialized NArray with given size along each dimension.
  NArray(Index size, const T& value = T());
  NArray(const NArray&) = default;
  NArray(NArray&&) = default;
  NArray& operator=(const NArray&) = default;
  NArray& operator=(NArray&&) = default;
  ~NArray() = default;

  /// Size along all dimensions.
  constexpr const Index& size() const { return m_size; }
  /// Read-only access element without boundary check.
  constexpr const T& operator[](Index idx) const { return m_data[linear(idx)]; }
  /// Access element without boundary check.
  constexpr T& operator[](Index idx) { return m_data[linear(idx)]; }
  /// Read-only access element with boundary check.
  const T& at(Index idx) const;
  /// Access element with boundary check.
  T& at(Index idx);

private:
  constexpr std::size_t linear(Index idx) const;
  constexpr bool within_bounds(Index idx) const;

  Index m_size;
  std::vector<T> m_data;
};

} // namespace
} // namespace histogram_impl

/// Uniform binning without under/overflow bins.
template<typename T>
class UniformAxis {
public:
  using Value = T;

  /// \param lower Lower inclusive boundary
  /// \param upper Upper exclusive boundary
  /// \param nbins Number of data bins within those boundaries
  UniformAxis(T lower, T upper, std::size_t nbins);

  /// Total number of bins along this axis including under/overflow bins.
  constexpr std::size_t nbins() const { return m_nbins; }
  /// Compute bin number for a test value.
  std::size_t index(T value) const;

private:
  std::size_t m_nbins;
  T m_lower;
  T m_upper;
};

/// Uniform binning with under/overflow bins.
///
/// The first and last bin index correspond to the under/overflow bin.
template<typename T>
class OverflowAxis {
public:
  using Value = T;

  /// \param lower Lower inclusive boundary
  /// \param upper Upper exclusive boundary
  /// \param nbins Number of data bins within those boundaries
  OverflowAxis(T lower, T upper, std::size_t nbins);

  /// Total number of bins along this axis including under/overflow bins.
  constexpr std::size_t nbins() const { return 2 + m_ndatabins; }
  /// Compute bin number for a test value.
  constexpr std::size_t index(T value) const;

private:
  std::size_t m_ndatabins;
  T m_lower;
  T m_upper;
};

/// Variable binninng defined by arbitrary bin edges.
template<typename T>
class VariableAxis {
public:
  using Value = T;

  /// \param edges Bin edges, lower ones inclusive, upper ones exclusive.
  explicit VariableAxis(std::vector<T>&& edges);
  /// \param edges Bin edges, lower ones inclusive, upper ones exclusive.
  VariableAxis(std::initializer_list<T> edges);

  /// Total number of bins along this axis including under/overflow bins.
  constexpr std::size_t nbins() const { return m_edges.size() - 1; }
  /// Compute bin number for a test value.
  std::size_t index(T value) const;

private:
  std::vector<T> m_edges;
};

/// A generic histogram with configurable axes.
///
/// \tparam T    The type of the data stored per bin
/// \tparam Axes Types must provide `::Value`, `.nbins()` and `.index(...)`
template<typename T, typename... Axes>
class Histogram {
public:
  using Data = dfe::histogram_impl::NArray<T, sizeof...(Axes)>;
  using Index = typename Data::Index;

  Histogram(Axes&&... axes);

  /// Get the number of bins along all axes.
  constexpr const Index& size() const { return m_data.size(); }
  /// Get the current entry value in the given bin.
  const T& value(Index idx) const { return m_data.at(idx); }
  /// Fill an entry into the histogram.
  ///
  /// \param values
  /// \param weight Associated weight, the default of 1 just counts entries.
  void fill(typename Axes::Value... values, T weight = static_cast<T>(1)) {
    // TODO 2018-11-28 how to typedef parameter pack Axes::Value...?
    m_data[index(std::index_sequence_for<Axes...>(), values...)] += weight;
  }

private:
  template<std::size_t... Is>
  constexpr Index index(
    std::index_sequence<Is...>, typename Axes::Value... values) const {
    return Index{std::get<Is>(m_axes).index(values)...};
  }

  Data m_data;
  std::tuple<Axes...> m_axes;
};

// predefined histogram types

using Histogram1 = Histogram<double, OverflowAxis<double>>;
using Histogram2 =
  Histogram<double, OverflowAxis<double>, OverflowAxis<double>>;

// implementation NArray

template<typename T, std::size_t NDimensions>
inline histogram_impl::NArray<T, NDimensions>::NArray(
  Index size, const T& value)
  : m_size(size)
  , m_data(
      std::accumulate(
        size.begin(), size.end(), static_cast<std::size_t>(1),
        std::multiplies<std::size_t>()),
      value) {}

// construct linear column-major index from n-dimensional index
template<typename T, std::size_t NDimensions>
constexpr std::size_t
histogram_impl::NArray<T, NDimensions>::linear(Index idx) const {
  std::size_t result = 0;
  std::size_t step = 1;
  for (std::size_t i = 0; i < NDimensions; ++i) {
    result += step * idx[i];
    step *= m_size[i];
  }
  return result;
}

template<typename T, std::size_t NDimensions>
constexpr bool
histogram_impl::NArray<T, NDimensions>::within_bounds(Index idx) const {
  for (std::size_t i = 0; i < NDimensions; ++i) {
    if (m_size[i] <= idx[i]) {
      return false;
    }
  }
  return true;
}

template<typename T, std::size_t NDimensions>
inline const T&
histogram_impl::NArray<T, NDimensions>::at(Index idx) const {
  if (!within_bounds(idx)) {
    throw std::out_of_range("NArray index is out of valid range");
  }
  return m_data[linear(idx)];
}

template<typename T, std::size_t NDimensions>
inline T&
histogram_impl::NArray<T, NDimensions>::at(Index idx) {
  if (!within_bounds(idx)) {
    throw std::out_of_range("NArray index is out of valid range");
  }
  return m_data[linear(idx)];
}

// implementation UniformAxis

template<typename T>
inline UniformAxis<T>::UniformAxis(T lower, T upper, std::size_t nbins)
  : m_nbins(nbins), m_lower(lower), m_upper(upper) {}

template<typename T>
inline std::size_t
UniformAxis<T>::index(T value) const {
  if (value < this->m_lower) {
    throw std::out_of_range("Value is smaller than lower axis limit");
  }
  if (m_upper <= value) {
    throw std::out_of_range("Value is equal or larger than upper axis limit");
  }
  // cast truncates to integer part; should work since index is always > 0.
  return static_cast<std::size_t>(
    m_nbins * (value - m_lower) / (m_upper - m_lower));
}

// implementation OverflowAxis

template<typename T>
inline OverflowAxis<T>::OverflowAxis(
  Value lower, Value upper, std::size_t nbins)
  : m_ndatabins(nbins), m_lower(lower), m_upper(upper) {}

template<typename T>
constexpr std::size_t
OverflowAxis<T>::index(T value) const {
  if (value < m_lower) {
    return 0;
  }
  if (m_upper <= value) {
    return m_ndatabins + 1;
  }
  // cast truncates to integer part; should work since index is always > 0.
  return 1
         + static_cast<std::size_t>(
           m_ndatabins * (value - m_lower) / (m_upper - m_lower));
}

// implementation VariableAxis

template<typename T>
inline VariableAxis<T>::VariableAxis(std::vector<Value>&& edges)
  : m_edges(std::move(edges)) {
  if (m_edges.size() < 2) {
    throw std::invalid_argument("Less than two bin edges");
  }
  // edges must be sorted and unique
  if (!std::is_sorted(
        m_edges.begin(), m_edges.end(), std::less_equal<Value>())) {
    throw std::invalid_argument("Bin edges are not sorted or have duplicates");
  }
}

template<typename T>
inline VariableAxis<T>::VariableAxis(std::initializer_list<Value> edges)
  : VariableAxis(std::vector<Value>(edges)) {}

template<typename T>
inline std::size_t
VariableAxis<T>::index(T value) const {
  // find upper edge of the corresponding bin
  auto it = std::upper_bound(m_edges.begin(), m_edges.end(), value);
  if (it == m_edges.begin()) {
    throw std::out_of_range("Value is smaller than lower axis limit");
  }
  if (it == m_edges.end()) {
    throw std::out_of_range("Value is equal or larger than upper axis limit");
  }
  return std::distance(m_edges.begin(), it) - 1;
}

// implementation Histogram

template<typename T, typename... Axes>
inline Histogram<T, Axes...>::Histogram(Axes&&... axes)
  // access nbins *before* moving the axes, otherwise the axes are invalid.
  : m_data(Index{axes.nbins()...}, static_cast<T>(0))
  , m_axes(std::move(axes)...) {}

} // namespace dfe
