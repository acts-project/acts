// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <any>
#include <cassert>
#include <memory>
#include <vector>

namespace Acts::detail {

/// Interface of the dynamic vector column used in the track container backend.
/// The interface provides basic methods to retrieve the data in a type erase
/// form. Further methods for memory pre-allocation and memory release are
/// defined
class DynamicColumnBase {
 public:
  virtual ~DynamicColumnBase() = default;

  /// Retrieve the stored data for the i-th proxy object stored in the container
  /// @param i: The index of the object within the container
  /// @return The stored data in a type-erased form for mutable access
  virtual std::any get(std::size_t i) = 0;

  /// Retrieve the stored data for the i-th proxy object stored in the container
  /// @param i: The index of the object within the container
  /// @return The stored data in a type-erased form for read access
  virtual std::any get(std::size_t i) const = 0;

  /// Allocate new memory for the proxy object held by the container
  virtual void add() = 0;

  /// Release all stored memory
  virtual void clear() = 0;

  /// Pre-reserve the memory for later allocation
  /// @param size: The number of memory blocks to allocate
  virtual void reserve(std::size_t size) = 0;

  /// Allocate memory for a given number of proxy objects
  /// @param size: The number of elements to be allocated
  virtual void resize(std::size_t size) = 0;

  /// Remove the memory associated with the i-th proxy object from the column
  /// @param i: The index of the object within the container
  virtual void erase(std::size_t i) = 0;

  /// The number of the currently allocated elements
  /// @return The size of the underlying storage
  virtual std::size_t size() const = 0;

  /// Copy the memory stored at the srcIdx of the passed DynamicColumn
  /// to the dstIdx of this instance
  /// @param dstIdx: Index of the memory slot of this instance
  ///                the copied information will be stored
  /// @param src: The source from where the memory can be copy
  /// @param srcIdx: Index at which element the memory to copy resides in src
  virtual void copyFrom(std::size_t dstIdx, const DynamicColumnBase& src,
                        std::size_t srcIdx) = 0;

  /// Copy some external data and store it in this instance
  /// @param dstIdx: Index if the memory slote where the passed
  ///                information will be stored
  /// @param srcPtr: Type-erased reference to the information to store
  virtual void copyFrom(std::size_t dstIdx, const std::any& srcPtr) = 0;

  /// Create a clone of this DynamicColumnBase instance
  /// @param empty: If toggled to true the content will not be
  ///               copied to the clone
  /// @returns A unique_ptr to a fresh DynamicColumnBase instance
  virtual std::unique_ptr<DynamicColumnBase> clone(
      bool empty = false) const = 0;
};

/// Actual implementation of the Column memory using a std::vector as backend.
/// The template parameter can be anything which is default constrictible and
/// copy assignable.
/// @tparam T: Data type of the data to be stored by the column
template <typename T>
class DynamicColumn : public DynamicColumnBase {
 public:
  /// @copydoc DynamicColumnBase::get
  std::any get(std::size_t i) override {
    assert(i < m_vector.size() && "DynamicColumn out of bounds");
    return &m_vector[i];
  }

  /// @copydoc DynamicColumnBase::get
  std::any get(std::size_t i) const override {
    assert(i < m_vector.size() && "DynamicColumn out of bounds");
    return &m_vector[i];
  }

  /// @copydoc DynamicColumnBase::add
  void add() override { m_vector.emplace_back(); }

  /// @copydoc DynamicColumnBase::clear
  void clear() override { m_vector.clear(); }

  /// @copydoc DynamicColumnBase::resize
  void resize(std::size_t size) override { m_vector.resize(size); }

  /// @copydoc DynamicColumnBase::reserve
  void reserve(std::size_t size) override { m_vector.reserve(size); }

  /// @copydoc DynamicColumnBase::erase
  void erase(std::size_t i) override { m_vector.erase(m_vector.begin() + i); }

  /// @copydoc DynamicColumnBase::size
  std::size_t size() const override { return m_vector.size(); }

  /// @copydoc DynamicColumnBase::clone
  std::unique_ptr<DynamicColumnBase> clone(bool empty) const override {
    if (empty) {
      return std::make_unique<DynamicColumn<T>>();
    }
    return std::make_unique<DynamicColumn<T>>(*this);
  }

  /// @copydoc DynamicColumnBase::copyFrom
  void copyFrom(std::size_t dstIdx, const DynamicColumnBase& src,
                std::size_t srcIdx) override {
    const auto* other = dynamic_cast<const DynamicColumn<T>*>(&src);
    assert(other != nullptr &&
           "Source column is not of same type as destination");
    m_vector.at(dstIdx) = other->m_vector.at(srcIdx);
  }

  /// @copydoc DynamicColumnBase::copyFrom
  void copyFrom(std::size_t dstIdx, const std::any& srcPtr) override {
    const auto* other = std::any_cast<const T*>(srcPtr);
    assert(other != nullptr &&
           "Source column is not of same type as destination");
    m_vector.at(dstIdx) = *other;
  }

 private:
  std::vector<T> m_vector;
};

/// Template specification to circumvent the flaws of the bool implementation
/// of a std::vector.
template <>
class DynamicColumn<bool> : public DynamicColumnBase {
 public:
  /// @copydoc DynamicColumnBase::get
  std::any get(std::size_t i) override {
    assert(i < m_vector.size() && "DynamicColumn out of bounds");
    return &m_vector[i].value;
  }

  /// @copydoc DynamicColumnBase::get
  std::any get(std::size_t i) const override {
    assert(i < m_vector.size() && "DynamicColumn out of bounds");
    return &m_vector[i].value;
  }

  /// @copydoc DynamicColumnBase::add
  void add() override { m_vector.emplace_back(); }

  /// @copydoc DynamicColumnBase::reserve
  void reserve(std::size_t size) override { m_vector.reserve(size); }

  /// @copydoc DynamicColumnBase::resize
  void resize(std::size_t size) override { m_vector.resize(size); }

  /// @copydoc DynamicColumnBase::clear
  void clear() override { m_vector.clear(); }

  /// @copydoc DynamicColumnBase::erase
  void erase(std::size_t i) override { m_vector.erase(m_vector.begin() + i); }

  /// @copydoc DynamicColumnBase::size
  std::size_t size() const override { return m_vector.size(); }

  /// @copydoc DynamicColumnBase::clone
  std::unique_ptr<DynamicColumnBase> clone(bool empty) const override {
    if (empty) {
      return std::make_unique<DynamicColumn<bool>>();
    }
    return std::make_unique<DynamicColumn<bool>>(*this);
  }

  /// @copydoc DynamicColumnBase::copyFrom
  void copyFrom(std::size_t dstIdx, const DynamicColumnBase& src,
                std::size_t srcIdx) override {
    const auto* other = dynamic_cast<const DynamicColumn<bool>*>(&src);
    assert(other != nullptr &&
           "Source column is not of same type as destination");
    m_vector.at(dstIdx) = other->m_vector.at(srcIdx);
  }

  /// @copydoc DynamicColumnBase::copyFrom
  void copyFrom(std::size_t dstIdx, const std::any& srcPtr) override {
    const auto* other = std::any_cast<const bool*>(srcPtr);
    assert(other != nullptr &&
           "Source column is not of same type as destination");
    m_vector.at(dstIdx).value = *other;
  }

 private:
  /// Auxiliary struct to wrap a boolean for use in a vector
  struct Wrapper {
    bool value{false};
  };
  std::vector<Wrapper> m_vector;
};

}  // namespace Acts::detail
