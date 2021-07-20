// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <type_traits>

namespace Acts {

namespace detail {
namespace PolymorphicValue {

/// @brief Control block base class which defines the needed interface
///
/// @tparam T The returnable type
template <typename T>
struct ControlBlockBase {
  virtual ~ControlBlockBase() = default;

  virtual std::unique_ptr<ControlBlockBase<T>> clone() const = 0;

  virtual T* pointer() = 0;
  virtual const T* pointer() const = 0;

  virtual T* releaseValue() = 0;
};

/// @brief Control block which directly stores a @c std::unique_ptr to the
/// actual value
///
/// @tparam T The returnable type
/// @tparam U The concrete type whose copy behavior is preserved
template <typename T, typename U>
struct ControlBlock : public ControlBlockBase<T> {
  static_assert(std::is_convertible_v<U*, T*>, "Types must be convertible");

  explicit ControlBlock(std::unique_ptr<U> u) : m_value{std::move(u)} {}

  std::unique_ptr<ControlBlockBase<T>> clone() const override {
    return std::make_unique<ControlBlock>(std::make_unique<U>(*m_value));
  }

  T* pointer() override { return m_value.get(); }
  const T* pointer() const override { return m_value.get(); }

  T* releaseValue() override { return m_value.release(); }

  std::unique_ptr<U> m_value;
};

/// @brief Control block which wraps another control block and delegates to it
/// value
///
/// @tparam T The returnable type
/// @tparam U The concrete type whose copy behavior is preserved
template <typename T, typename U>
struct DelegatingControlBlock : public ControlBlockBase<T> {
  DelegatingControlBlock(std::unique_ptr<ControlBlockBase<U>> delegate)
      : m_delegate{std::move(delegate)} {}

  std::unique_ptr<ControlBlockBase<T>> clone() const override {
    return std::make_unique<DelegatingControlBlock>(m_delegate->clone());
  }

  T* pointer() override { return m_delegate->pointer(); }
  const T* pointer() const override { return m_delegate->pointer(); }

  T* releaseValue() override { return m_delegate->releaseValue(); }

  std::unique_ptr<ControlBlockBase<U>> m_delegate;
};

/// @brief Control block which throws a runtime exception when the contained
/// value is to be cloned
///
/// @tparam T The returnable type
/// @tparam U The concrete type
template <typename T, typename U>
struct NonCopyableControlBlock : public ControlBlockBase<T> {
  static_assert(std::is_convertible_v<U*, T*>, "Types must be convertible");

  explicit NonCopyableControlBlock(std::unique_ptr<U> u)
      : m_value{std::move(u)} {}

  std::unique_ptr<ControlBlockBase<T>> clone() const override {
    throw std::runtime_error{"Stored PolymorphicValue is not copyable"};
  }

  T* pointer() override { return m_value.get(); }
  const T* pointer() const override { return m_value.get(); }

  T* releaseValue() override { return m_value.release(); }

  std::unique_ptr<U> m_value;
};

}  // namespace PolymorphicValue
}  // namespace detail

/// @brief Type trait to detect whether a type is a polymorphic value
///
/// @tparam T The type to check
template <class T>
class PolymorphicValue;

/// @cond
/// @brief Base specialization: false
///
/// @tparam T The type to check
template <class T>
struct IsPolymorphicValue : std::false_type {};
/// @endcond

///@brief The true specialization
///
///@tparam T The type to test
template <class T>
struct IsPolymorphicValue<PolymorphicValue<T>> : std::true_type {};

/// @brief Polymorphic value container
///
/// This container is similar to @c std::unique_ptr in that it manages the
/// lifetime of a contained object. In contrast to it, this type retains
/// knowledge of how to copy the stored type, and can be freely be copied in
/// turn.
///
/// @note Any type deriving from @c T can be stored in the polymorphic value
/// that is created this way
///
///@tparam T The type to store.
template <typename T>
class PolymorphicValue {
  // Friend declarations that are needed
  template <typename U>
  friend class PolymorphicValue;

  /// @cond
  template <typename U, typename... Args>
  friend PolymorphicValue<U> makePolymorphicValue(Args&&... args);

  template <typename U, typename V, typename... Args>
  friend PolymorphicValue<U> makePolymorphicValue(Args&&... args);
  /// @endcond

  /// @brief Private constructor to construct a new Polymorphic Value object in
  /// place
  ///
  /// @note This constructor is private. Use @c Acts::makePolymorphicValue from
  /// external
  /// @tparam U The derived type to construct
  /// @tparam Args Template pack to construct the value with
  /// @param args Parameter pack to construct with
  template <typename U, typename _U = U, typename _T = T,
            typename = std::enable_if_t<std::is_convertible_v<_U*, _T*> &&
                                        !IsPolymorphicValue<_U>::value>,
            typename... Args>
  explicit PolymorphicValue(std::in_place_type_t<U>, Args&&... args)
      : m_controlBlock{std::make_unique<
            detail::PolymorphicValue::ControlBlock<T, U>>(
            std::make_unique<U>(std::forward<Args>(args)...))},
        m_pointer{m_controlBlock->pointer()} {}

 public:
  /// @brief Default constructor
  PolymorphicValue() {}

  /// @brief Construct a new Polymorphic Value object from an rvalue reference
  ///
  /// @tparam U Concrete type of the rvalue reference
  /// @param u The rvalue reference argument
  template <typename U, typename _U = U, typename _T = T,
            typename = std::enable_if_t<std::is_convertible_v<_U*, _T*> &&
                                        !IsPolymorphicValue<_U>::value>>
  explicit PolymorphicValue(U&& u)
      : m_controlBlock{std::make_unique<
            detail::PolymorphicValue::ControlBlock<T, U>>(
            std::make_unique<U>(u))},
        m_pointer{m_controlBlock->pointer()} {}

  /// @brief Construct a new Polymorphic Value object from a raw pointer
  ///
  /// @note @c PolymorphicValue assumes ownership of the raw pointer
  /// @note If @c U is **not** copyable, this @c PolymorphicValue will also not
  /// be copyable. This does not lead to a compile-time error, but will throw a
  /// runtime exception when a copy is requested. This is needed so that @c
  /// PolymorphicValue can assume ownership of a raw pointer, but cannot copy it
  /// through e.g. a pointer to a virtual base class, even though derived
  /// classes are copyable.
  ///
  /// @tparam U Concrete type of the pointed at variable
  /// @param u The pointer argument
  template <typename U, typename _U = U, typename _T = T,
            typename = std::enable_if_t<std::is_convertible_v<_U*, _T*> &&
                                        !IsPolymorphicValue<_U>::value>>
  explicit PolymorphicValue(U* u) {
    if constexpr (std::is_copy_constructible_v<U>) {
      m_controlBlock =
          std::make_unique<detail::PolymorphicValue::ControlBlock<T, U>>(
              std::unique_ptr<U>(u));
    } else {
      m_controlBlock = std::make_unique<
          detail::PolymorphicValue::NonCopyableControlBlock<T, U>>(
          std::unique_ptr<U>(u));
    }
    m_pointer = m_controlBlock->pointer();
  }

  /// @brief Construct a new Polymorphic Value object from another @c
  /// PolymorphicValue by copying
  ///
  /// @note This only works if the base type of the other polymorphic value can
  /// be converted to the base type of this polymorphic value.
  ///
  /// @tparam U The base type of the other polymorphic value
  /// @param other The other polymorphic value
  template <typename U, typename _U = U, typename _T = T,
            typename = std::enable_if_t<!std::is_same_v<_U, _T> &&
                                        std::is_convertible_v<_U*, _T*>>>
  PolymorphicValue(const PolymorphicValue<U>& other) {
    if (!other.m_controlBlock) {
      reset();
      return;
    }
    auto cbTmp = std::move(m_controlBlock);
    m_controlBlock = std::make_unique<
        detail::PolymorphicValue::DelegatingControlBlock<T, U>>(
        other.m_controlBlock->clone());
    m_pointer = m_controlBlock->pointer();
  }

  /// @brief Construct a new Polymorphic Value object from an rvalue reference
  /// to another @c PolymorphicValue
  ///
  /// @note This only works if the base type of the other polymorphic value can
  /// be converted to the base type of this polymorphic value.
  ///
  /// @tparam U The base type of the other polymorphic value
  /// @param other The other polymorphic value
  template <typename U, typename _U = U, typename _T = T,
            typename = std::enable_if_t<!std::is_same_v<_U, _T> &&
                                        std::is_convertible_v<_U*, _T*>>>
  PolymorphicValue(PolymorphicValue<U>&& other) {
    if (!other.m_controlBlock) {
      reset();
      return;
    }
    auto cbTmp = std::move(m_controlBlock);
    m_controlBlock = std::make_unique<
        detail::PolymorphicValue::DelegatingControlBlock<T, U>>(
        std::move(other.m_controlBlock));
    m_pointer = m_controlBlock->pointer();
    other.m_pointer = nullptr;
  }

  /// @brief Copy assignment from a concrete value
  ///
  /// @note This only works if the given type @c U can
  /// be converted to the base type of this polymorphic value.
  ///
  /// @tparam U The concrete type to copy-assign
  /// @param u The concrete type
  template <typename U, typename _U = U, typename _T = T,
            typename = std::enable_if_t<!std::is_same_v<_U, _T> &&
                                        std::is_convertible_v<_U*, _T*> &&
                                        !IsPolymorphicValue<_U>::value>>
  PolymorphicValue& operator=(const U& u) {
    auto cbTmp =
        std::move(m_controlBlock);  // extend lifetime until after operation
    m_controlBlock =
        std::make_unique<detail::PolymorphicValue::ControlBlock<T, U>>(
            std::make_unique<U>(u));
    m_pointer = m_controlBlock->pointer();
    return *this;
  }

  /// @brief Move assignment from a concrete value
  ///
  /// @note This only works if the given type @c U can
  /// be converted to the base type of this polymorphic value.
  ///
  /// @tparam U The concrete type to move-assign
  /// @param u The concrete type
  template <typename U, typename _U = U, typename _T = T,
            typename = std::enable_if_t<std::is_convertible_v<_U*, _T*> &&
                                        !IsPolymorphicValue<_U>::value>>
  PolymorphicValue& operator=(U&& u) {
    auto cbTmp = std::move(m_controlBlock);
    m_controlBlock =
        std::make_unique<detail::PolymorphicValue::ControlBlock<T, U>>(
            std::make_unique<U>(std::move(u)));
    m_pointer = m_controlBlock->pointer();
    return *this;
  }

  /// @brief Copy assignment from another @c PolymorphicValue with a different
  /// base type
  ///
  /// @note This only works if the other base type @c U can
  /// be converted to the base type of this polymorphic value.
  ///
  /// @tparam U The base type of the other polymorphic value
  /// @param other The other polymorphic value
  template <typename U, typename _U = U, typename _T = T,
            typename = std::enable_if_t<!std::is_same_v<_U, _T> &&
                                        std::is_convertible_v<_U*, _T*>>>
  PolymorphicValue& operator=(const PolymorphicValue<U>& other) {
    if (!other.m_controlBlock) {
      reset();
      return *this;
    }
    auto cbTmp = std::move(m_controlBlock);
    m_controlBlock = std::make_unique<
        detail::PolymorphicValue::DelegatingControlBlock<T, U>>(
        other.m_controlBlock->clone());
    m_pointer = m_controlBlock->pointer();
    return *this;
  }

  /// @brief Move assignment from another @c PolymorphicValue with a different
  /// base type
  ///
  /// @note This only works if the other base type @c U can
  /// be converted to the base type of this polymorphic value.
  ///
  /// @tparam U The base type of the other polymorphic value
  /// @param other The other polymorphic value
  template <typename U, typename _U = U, typename _T = T,
            typename = std::enable_if_t<!std::is_same_v<_U, _T> &&
                                        std::is_convertible_v<_U*, _T*>>>
  PolymorphicValue& operator=(PolymorphicValue<U>&& other) {
    if (!other.m_controlBlock) {
      reset();
      return *this;
    }
    auto cbTmp = std::move(m_controlBlock);
    m_controlBlock = std::make_unique<
        detail::PolymorphicValue::DelegatingControlBlock<T, U>>(
        std::move(other.m_controlBlock));
    m_pointer = m_controlBlock->pointer();
    other.m_pointer = nullptr;
    return *this;
  }

  /// @brief Copy constructor from another @c PolymorphicValue of the same base
  /// type as this one
  ///
  /// @param other The other polymocphic value
  PolymorphicValue(const PolymorphicValue& other) {
    if (!other.m_controlBlock) {
      reset();
      return;
    }
    auto cbTmp = std::move(m_controlBlock);
    m_controlBlock = other.m_controlBlock->clone();
    m_pointer = m_controlBlock->pointer();
  }

  /// @brief Copy assignment from another @c PolymorphicValue of the same base
  /// type as this one
  ///
  /// @param other The other polymocphic value
  PolymorphicValue& operator=(const PolymorphicValue& other) {
    if (!other.m_controlBlock) {
      reset();
      return *this;
    }
    auto cbTmp = std::move(m_controlBlock);
    m_controlBlock = other.m_controlBlock->clone();
    m_pointer = m_controlBlock->pointer();
    return *this;
  }

  /// @brief Move constructor from another @c PolymorphicValue of the same base
  /// type as this one
  ///
  /// @param other The other polymocphic value
  PolymorphicValue(PolymorphicValue&& other) {
    if (!other.m_controlBlock) {
      reset();
      return;
    }
    auto cbTmp = std::move(m_controlBlock);
    m_controlBlock = std::move(other.m_controlBlock);
    m_pointer = m_controlBlock->pointer();
    other.m_pointer = nullptr;
  }

  /// @brief Move assignment from another @c PolymorphicValue of the same base
  /// type as this one
  ///
  /// @param other The other polymocphic value
  PolymorphicValue& operator=(PolymorphicValue&& other) {
    if (!other.m_controlBlock) {
      reset();
      return *this;
    }
    auto cbTmp = std::move(m_controlBlock);
    m_controlBlock = std::move(other.m_controlBlock);
    m_pointer = m_controlBlock->pointer();
    other.m_pointer = nullptr;
    return *this;
  }

  /// @brief Check whether this @c PolymorphicValue currently contains a value
  ///
  /// @return true If a value is contained
  /// @return false If no value is contained
  operator bool() const { return !!m_controlBlock; }

  /// @brief Release ownership of the contained type
  ///
  /// @return T* The raw pointer to the contained type. This type is now owned
  /// by the caller!
  T* release() {
    // release the value itself first
    T* value = m_controlBlock->releaseValue();
    // now release the control block
    m_controlBlock.reset();
    m_pointer = nullptr;
    return value;
  }

  /// @brief Reset this @c PolymorphicValue to an already existing object in the
  /// form of a raw pointer
  ///
  /// @tparam U The type of the pointer to reset to
  /// @param u The pointer to reset to
  template <typename U, typename _U = U, typename _T = T,
            typename = std::enable_if_t<std::is_copy_constructible_v<_U> &&
                                        std::is_convertible_v<_U*, _T*>>>
  void reset(U* u) {
    auto cbTmp =
        std::move(m_controlBlock);  // extend lifetime until after operation
    m_controlBlock =
        std::make_unique<detail::PolymorphicValue::ControlBlock<T, U>>(
            std::unique_ptr<U>(u));
    m_pointer = m_controlBlock->pointer();
  }

  /// @brief Clear the value of this @c PolymorphicValue
  void reset() {
    m_controlBlock.reset();
    m_pointer = nullptr;
  }

  /// @brief Dereference operator
  ///
  /// @return T* Pointer to the contained value
  T* operator->() {
    assert(m_controlBlock && m_pointer != nullptr);
    return m_pointer;
  }

  /// @brief Dereference operator
  ///
  /// @return T& Reference to the contained value
  T& operator*() {
    assert(m_controlBlock && m_pointer != nullptr);
    return *m_pointer;
  }

  /// @brief Const dereference operator
  ///
  /// @return T* Pointer to the contained value
  const T* operator->() const {
    assert(m_controlBlock && m_pointer != nullptr);
    return m_pointer;
  }

  /// @brief Const dereference operator
  ///
  /// @return T& Reference to the contained value
  const T& operator*() const {
    assert(m_controlBlock && m_pointer != nullptr);
    return *m_pointer;
  }

  /// @brief Getter for the raw pointer contained. Needed for interface
  /// compatibility with other smart pointer types.
  ///
  /// @return T* The raw pointer to the value contained
  T* get() { return m_pointer; }
  const T* get() const { return m_pointer; }

 private:
  std::unique_ptr<detail::PolymorphicValue::ControlBlockBase<T>> m_controlBlock{
      nullptr};
  T* m_pointer{nullptr};
};

/// @cond
/// @brief Factory function for a polymorphic value from constructor arguments
///
/// @tparam T The base type of the returned polymorphic value
/// @tparam U The concrete type to construct in the polymorphic value
/// @tparam Args Template parameter pack to construct from
/// @param args Parameter pack to construct from
/// @return PolymorphicValue<T> The constructed polymorphic value
template <typename T, typename U, typename... Args>
PolymorphicValue<T> makePolymorphicValue(Args&&... args) {
  return PolymorphicValue<T>{std::in_place_type_t<U>(),
                             std::forward<Args>(args)...};
}
/// @endcond

/// @brief Factory function for a polymorphic value from constructor arguments
///
/// @tparam T The type to construct, this will also be the base type of the
/// polymorphic value
/// @tparam Args Template parameter pack to construct from
/// @param args Parameter pack to construct from
/// @return PolymorphicValue<T> The constructed polymorphic value
template <typename T, typename... Args>
PolymorphicValue<T> makePolymorphicValue(Args&&... args) {
  return PolymorphicValue<T>{std::in_place_type_t<T>(),
                             std::forward<Args>(args)...};
}

}  // namespace Acts