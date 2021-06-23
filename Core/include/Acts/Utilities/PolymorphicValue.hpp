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
};

template <typename T, typename U = T>
struct ControlBlock {
  explicit ControlBlock(std::unique_ptr<U> value) : m_value{std::move(value)} {}

  std::unique_ptr<ControlBlockBase<T>> clone() const override {
    // return std::make_unique<ControlBlock<U>>(*this);
    return std::make_unique<ControlBlock<U>>(*m_value);
  }

  T* pointer() { return &m_value; }

  std::unique_ptr<U> m_value;
};

}  // namespace detail

template <typename T>
class PolymorphicValue {
 public:
  PolymorphicValue() {}

  //   template <typename U>
  //   explicit PolymorphicValue(U value) {}

  template <typename U>
  explicit PolymorphicValue(U&& value)
      : m_controlBlock{
            std::make_unique<detail::ControlBlock<U>>(std::move(value))} {}

  template <typename U>
  PolymorphicValue& operator=(const U& other) {
    return *this;
  }

  template <typename U>
  PolymorphicValue& operator=(U&& other) {
    return *this;
  }

  operator bool() const { return !!m_controlBlock; }

  T* operator->() {
    assert(m_controlBlock);
    return *m_controlBlock->pointer();
  }

  const T* operator->() const {
    assert(m_controlBlock);
    return *m_controlBlock->pointer();
  }

 private:
  std::unique_ptr<detail::ControlBlock<T>> m_controlBlock{nullptr};
};

}  // namespace Acts