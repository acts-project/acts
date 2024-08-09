// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <array>
#include <memory>
#include <stdexcept>
#include <vector>

namespace Acts {

template <class T, std::size_t D>
class DynamicArray {
 public:
  /// Abbrivation of the data-type being stored
  using value_type = T;
  ///  Abbrivation of the underlying vector type
  using vec_type = std::vector<value_type>;

  /**
   * @brief Construct a new Dynamical Array object without allocating any memory.
   *        Any trial to access the data later will lead to an exception later
   */
  DynamicArray() = default;

  /**
   * @brief Construct a new Dynamical Array object with the specification of the array ranges in all dimensions         *
   * @tparam args  list of integers specyfing the length in each direction. E.g.
   *               DynamicArray<double,4> array{4,7,9};
   *               constructs an array with maximum 4, 7 and 9 elements along
   * the first, second and third direction
   */
  template <typename... args>
  DynamicArray(args... jks);
  /// Standard assignment operator
  DynamicArray& operator=(const DynamicArray& other) = default;
  /// Standard copy constructor
  DynamicArray(const DynamicArray& other) = default;
  /// Standard move assignment operator
  DynamicArray& operator=(DynamicArray&& other) = default;
  /// Standard move constructor
  DynamicArray(DynamicArray&& other) = default;

  /// Standard begin iterator provided by the underlying vector
  typename vec_type::iterator begin();
  typename vec_type::const_iterator begin() const;
  /// Standard end iterator
  typename vec_type::iterator end();
  typename vec_type::const_iterator end() const;

  /// Reverse begin iterator
  typename vec_type::reverse_iterator rbegin();
  typename vec_type::const_reverse_iterator rbegin() const;
  /// Reverse end iterator
  typename vec_type::reverse_iterator rend();
  typename vec_type::const_reverse_iterator rend() const;

  /// Returns the total number of elements that can be stored in the container
  std::size_t size() const;
  /// Returns the array containing the defined stripes in each direction
  const std::array<std::size_t, D>& getStripes() const;
  /// Returns the array containing the offsets in memory between two neighboring
  /// elements in the same dimension
  const std::array<std::size_t, D>& getMemoryOffsets() const;

  /**
   * @brief Getter methods to access the data in (non-) const mode using the D- local dimensional indices. The number of indices given
   * to the method must match the dimension of the array
   */
  template <typename... args>
  value_type& get(args... jks);
  template <typename... args>
  const value_type& get(args... jks) const;
  /**
   * @brief Translation of the D- local dimensional indices to the global index finding the corresponding element in memory
   *
   */
  inline std::size_t index(const std::array<std::size_t, D>& indices) const;
  template <typename... args>
  std::size_t index(args... jks) const;
  /**
   *  @brief Translation of the global index into the axis-index
   */
  std::array<std::size_t, D> axisIndices(std::size_t globIdx) const;
  /**
   * @brief Getter methods using the global index if the data is represented as a 1D array.
   *
   */
  value_type& get_val(std::size_t idx);
  const value_type& get_val(std::size_t idx) const;

  /// Sets all elements contained in the array to the passed value
  void assign(const T& value);

  /**
   *  @brief method to change the memory layout of the array after the array is constructed
   *          e.g. DynamicArray<int,3> array;
   *               array.changeStripes(2,4,5);
   *          will assign memory to store 2x4x5 elements. If one of the values
   * is zero, then the memory will be cleared
   *
   */
  template <typename... args>
  void changeStripes(args... jks);
  void changeStripes(const std::array<std::size_t, D>& stripes);

  /// Standard logical operators. The operator() returns true if the array is
  /// not empty.
  operator bool() const;
  bool operator!() const;
  bool empty() const;

  /**
   * @brief The (Const)ArrayNavigator is a lightweight helper class to navigate through the array via the multidimensional [][]
   * operator. For every D dimensional array, there exist D-1 ArrayNavigators
   * that create other ArrayNavigators having one dimension less until N = 1 is
   * reached. This navigator grants the direct access to the data. Each
   * Navigator stores the index in the corresponding bracket and is directly
   * connected to its parent which is the DynamicArray for D-1.
   */
  template <std::size_t N>
  class ArrayNavigator {
   public:
    using parent_type = typename std::conditional<N == D - 1, DynamicArray,
                                                  ArrayNavigator<N + 1>>::type;
    using value_type =
        typename std::conditional<N == 1, T&, ArrayNavigator<N - 1>>::type;

    ArrayNavigator(parent_type* parent, std::size_t idx);
    /// Access to the data
    value_type operator[](std::size_t i);
    /// Global index calculator upstream
    template <typename... args>
    std::size_t index(args... indices) const;
    /// Getter method to stream the data
    T& get_val(std::size_t idx);

   private:
    parent_type* m_parent{nullptr};
    std::size_t m_idx{0};
  };

  template <std::size_t N>
  class ConstArrayNavigator {
   public:
    using parent_type =
        typename std::conditional<N == D - 1, DynamicArray,
                                  ConstArrayNavigator<N + 1>>::type;
    using value_type =
        typename std::conditional<N == 1, const T&,
                                  ConstArrayNavigator<N - 1>>::type;

    ConstArrayNavigator(const parent_type* parent, std::size_t idx);
    /// Access to the data
    value_type operator[](std::size_t i) const;
    /// Global index calculator upstream
    template <typename... args>
    std::size_t index(args... indices) const;
    /// Getter method to stream the data
    const T& get_val(std::size_t idx) const;

   private:
    const parent_type* m_parent{nullptr};
    std::size_t m_idx{0};
  };

  /**
   * @brief Data access to the underlying data. It returns the i-th ArrayNavigator. If the index is out of bounds, i.e. i >=
   * stripSize() or the array does not store any data, an exception is thrown.
   * To access the specific element one needs to specify all indices
   * MultiDimensionalArray<double,5> arr{4,6,8,9,10};         *
   * arr[0][1][2][3][4] = 24.;
   * @param i: Index of the first requested row
   * @return array_type& : D-1 dimensional ArrayNavigator if D > 1 other wise the underlying data type T
   */

  /// Depending on the dimension of the array the [] operator needs to return
  /// either the lightweight ArrayNavigators or the value directly. The latter
  /// is only true if the dimension is one which make this class essentially to
  /// a std::vector<T>
  using array_type =
      typename std::conditional<D != 1, ArrayNavigator<D - 1>, T&>::type;
  using const_array_type =
      typename std::conditional<D != 1, ConstArrayNavigator<D - 1>,
                                const T&>::type;

  inline array_type operator[](std::size_t i);
  inline const_array_type operator[](std::size_t i) const;

 private:
  template <std::size_t pos, typename... args>
  std::size_t index(std::size_t arg_val, args... jks) const;
  template <std::size_t pos, typename... args>
  void assign_stripes(std::size_t val, args... jks);

  template <std::size_t pos>
  std::size_t index(std::size_t arg_val) const;
  template <std::size_t pos>
  void assign_stripes(std::size_t val);
  void allocate();

  vec_type m_data_ptr{};
  std::size_t m_size{0};
  std::array<std::size_t, D> m_stripe_lengths{};
  std::array<std::size_t, D> m_stripe_offsets{};
};

}  // namespace Acts

#include "Acts/Utilities/DynamicArray.ipp"
