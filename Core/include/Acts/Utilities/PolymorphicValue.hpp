// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts {

namespace detail {
template <typename T>
struct ControlBlockBase {
  virtual ~ControlBlockBase() = default;

  virtual std::unique_ptr<ControlBlockBase<T>> clone() const = 0;

  virtual T* pointer() = 0;
  virtual const T* pointer() const = 0;

  virtual T* releaseValue() = 0;
};

template <typename T, typename U>
struct ControlBlock : public ControlBlockBase<T> {
  explicit ControlBlock(std::unique_ptr<U> u) : m_value{std::move(u)} {}

  std::unique_ptr<ControlBlockBase<T>> clone() const override {
    // return std::make_unique<ControlBlock<U>>(*this);
    return std::make_unique<ControlBlock<T, U>>(std::make_unique<U>(*m_value));
  }

  T* pointer() override { return m_value.get(); }
  const T* pointer() const override { return m_value.get(); }

  T* releaseValue() override { return m_value.release(); }

  std::unique_ptr<U> m_value;
};

}  // namespace detail

template <typename T>
class PolymorphicValue {
 public:
  PolymorphicValue() {}

  //   template <typename U>
  //   explicit PolymorphicValue(U value) {}

  template <typename U, typename _U = U, typename _T = T,
            typename = std::enable_if<!std::is_same_v<_U, _T> &&
                                      std::is_copy_constructible_v<_U> &&
                                      std::is_convertible_v<_U*, _T*>>>
  explicit PolymorphicValue(U&& u)
      : m_controlBlock{std::make_unique<detail::ControlBlock<T, U>>(
            std::make_unique<U>(u))} {}

  template <typename U, typename _U = U, typename _T = T,
            typename = std::enable_if<!std::is_same_v<_U, _T> &&
                                      std::is_copy_constructible_v<_U> &&
                                      std::is_convertible_v<_U*, _T*>>,
            typename... Args>
  explicit PolymorphicValue(std::in_place_type_t<U>, Args&&... args)
      : m_controlBlock{std::make_unique<detail::ControlBlock<T, U>>(
            std::make_unique<U>(std::forward<Args>(args)...))} {}

  template <typename U, typename _U = U, typename _T = T,
            typename = std::enable_if<!std::is_same_v<_U, _T> &&
                                      std::is_copy_constructible_v<_U> &&
                                      std::is_convertible_v<_U*, _T*>>>
  PolymorphicValue& operator=(const U& u) {
    auto cbTmp =
        std::move(m_controlBlock);  // extend lifetime until after operation
    m_controlBlock =
        std::make_unique<detail::ControlBlock<T, U>>(std::make_unique<U>(u));
    return *this;
  }

  template <typename U, typename _U = U, typename _T = T,
            typename = std::enable_if<std::is_copy_constructible_v<_U> &&
                                      std::is_convertible_v<_U*, _T*>>>
  PolymorphicValue& operator=(U&& u) {
    auto cbTmp =
        std::move(m_controlBlock);  // extend lifetime until after operation
    m_controlBlock = std::make_unique<detail::ControlBlock<T, U>>(
        std::make_unique<U>(std::move(u)));
    return *this;
  }

  operator bool() const { return !!m_controlBlock; }

  T* release() {
    // release the value itself first
    T* value = m_controlBlock->releaseValue();
    // now release the control block
    m_controlBlock.reset();
    return value;
  }

  template <typename U, typename _U = U, typename _T = T,
            typename = std::enable_if<std::is_copy_constructible_v<_U> &&
                                      std::is_convertible_v<_U*, _T*>>>
  void reset(U* u) {
    auto cbTmp =
        std::move(m_controlBlock);  // extend lifetime until after operation
    m_controlBlock =
        std::make_unique<detail::ControlBlock<T, U>>(std::unique_ptr<U>(u));
  }

  void reset() { m_controlBlock.reset(); }

  T* operator->() {
    assert(m_controlBlock);
    return m_controlBlock->pointer();
  }

  const T* operator->() const {
    assert(m_controlBlock);
    return m_controlBlock->pointer();
  }

 private:
  std::unique_ptr<detail::ControlBlockBase<T>> m_controlBlock{nullptr};
};

}  // namespace Acts