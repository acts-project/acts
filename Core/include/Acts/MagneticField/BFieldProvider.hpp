// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/MagneticField/MagneticFieldContext.hpp"

#include <any>
#include <array>
#include <iostream>
#include <memory>

namespace Acts {

namespace detail {

/**
 * Two alternative opaque cache types are given, one using std::any, and the
 * other using a fixed maximum sized allocation on the stack. AnyCache did not
 * seem to perform worse even with heap allocation and checked casts. Both
 * implemenations are kept so we can evaluate in the future. It is selected
 * inside `BFieldBase`.
 */

/// Small opaque cache type which uses small buffer optimization
class SBOCache {
 public:
  template <typename T, typename... Args>
  static SBOCache make(Args&&... args) {
    SBOCache cache{};

    static_assert(std::is_same_v<T, std::decay_t<T>>,
                  "Please pass the raw type, no const or ref");
    static_assert(sizeof(T) <= sizeof(cache.m_data),
                  "Passed type is too large");

    /*T* ptr =*/new (cache.m_data.data()) T(std::forward<Args>(args)...);
    static Handler<T> static_handler{};
    cache.m_handler = &static_handler;

    return cache;
  }

  template <typename T>
  T& get() {
    static_assert(std::is_same_v<T, std::decay_t<T>>,
                  "Please pass the raw type, no const or ref");
    static_assert(sizeof(T) <= sizeof(m_data), "Passed type is too large");
    return *reinterpret_cast<T*>(m_data.data());
  }

  ~SBOCache() {
    assert(m_handler && "Handler cannot be dead");
    m_handler->destroy(m_data.data());
  }

  SBOCache(SBOCache&& other) {
    m_handler = other.m_handler;
    assert(m_handler && "Handler is null");
    m_handler->move(other.m_data.data(), m_data.data());
  }

  SBOCache& operator=(SBOCache&& other) {
    m_handler = other.m_handler;
    assert(m_handler && "Handler is null");
    m_handler->move(other.m_data.data(), m_data.data());
    return *this;
  }

  SBOCache(const SBOCache& other) {
    m_handler = other.m_handler;
    assert(m_handler && "Handler is null");
    m_handler->copy(other.m_data.data(), m_data.data());
  };

  SBOCache& operator=(const SBOCache& other) {
    m_handler = other.m_handler;
    assert(m_handler && "Handler is null");
    m_handler->copy(other.m_data.data(), m_data.data());
    return *this;
  };

 private:
  SBOCache(){};

  struct HandlerBase {
    virtual void destroy(void* ptr) = 0;
    virtual void move(void* from, void* to) = 0;
    virtual void copy(const void* from, void* to) = 0;
    virtual ~HandlerBase() = default;
  };

  template <typename T>
  struct Handler final : public HandlerBase {
    void destroy(void* ptr) override {
      assert(ptr != nullptr && "Address to destroy is nullptr");
      T* obj = static_cast<T*>(ptr);
      obj->~T();
    }

    void move(void* from, void* to) override {
      assert(from != nullptr && "Source is null");
      assert(to != nullptr && "Target is null");

      T* _from = static_cast<T*>(from);
      T* _to = static_cast<T*>(to);

      (*_to) = std::move(*_from);
    }

    void copy(const void* from, void* to) override {
      assert(from != nullptr && "Source is null");
      assert(to != nullptr && "Target is null");

      const T* _from = static_cast<const T*>(from);
      T* _to = static_cast<T*>(to);

      (*_to) = (*_from);
    }
  };

  alignas(std::max_align_t) std::array<char, 512> m_data;
  HandlerBase* m_handler{nullptr};
};

/// Opaque cache type using std::any under the hood
/// @note Currently unused
class AnyCache {
 public:
  template <typename T, typename... Args>
  static AnyCache make(Args&&... args) {
    AnyCache cache{};

    static_assert(std::is_same_v<T, std::decay_t<T>>,
                  "Please pass the raw type, no const or ref");

    cache.m_any.emplace<T>(std::forward<Args>(args)...);

    return cache;
  }

  template <typename T>
  T& get() {
    return std::any_cast<T&>(m_any);
  }

 private:
  AnyCache() = default;

  std::any m_any;
};
}  // namespace detail

/// Base class for all magnetic field providers
class BFieldProvider {
 public:
  using Cache = detail::SBOCache;

  /// @brief Make an opaque cache for the magnetic field
  ///
  /// @param mctx The magnetic field context to generate cache for
  /// @return Cache The opaque cache object
  virtual Cache makeCache(const MagneticFieldContext& mctx) const = 0;

  /// @brief retrieve magnetic field value
  ///
  /// @param [in] position global 3D position
  ///
  /// @return magnetic field vector at given position
  virtual Vector3 getField(const Vector3& position, Cache& cache) const = 0;

  /// @brief retrieve magnetic field value
  ///
  /// @param [in] position global 3D position
  /// @param [in,out] cache Cache object. Contains field cell used for
  /// interpolation
  ///
  /// @return magnetic field vector at given position
  virtual Vector3 getField(const Vector3& position) const = 0;

  /// @brief retrieve magnetic field value & its gradient
  ///
  /// @param [in]  position   global 3D position
  /// @param [out] derivative gradient of magnetic field vector as (3x3) matrix
  /// @return magnetic field vector
  virtual Vector3 getFieldGradient(const Vector3& position,
                                   ActsMatrix<3, 3>& derivative) const = 0;

  /// @brief retrieve magnetic field value & its gradient
  ///
  /// @param [in]  position   global 3D position
  /// @param [out] derivative gradient of magnetic field vector as (3x3) matrix
  /// @param [in,out] cache Cache object. Contains field cell used for
  /// @return magnetic field vector
  virtual Vector3 getFieldGradient(const Vector3& position,
                                   ActsMatrix<3, 3>& derivative,
                                   Cache& cache) const = 0;
};  // namespace Acts

}  // namespace Acts