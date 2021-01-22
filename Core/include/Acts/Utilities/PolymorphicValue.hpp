// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

namespace Acts {

template <typename T>
class PolymorphicValue {
 public:
  template <typename U, typename... Args>
  static PolymorphicValue<T> make(Args&&... args) {
    static_assert(std::is_base_of_v<T, U>,
                  "Given type is not derived from base type");
    PolymorphicValue<T> pm{std::make_unique<U>(std::forward<Args>(args)...)};

    static Handler<U> handler;
    pm.m_handler = &handler;

    return pm;
  }

  T* operator->() { return m_value.get(); }

  PolymorphicValue& operator=(const PolymorphicValue& other) {
    assert(m_handler != nullptr && "No handler set");
    m_handler->assign(other.m_value.get(), m_value.get());
    return *this;
  }

  PolymorphicValue(const PolymorphicValue& other) {
    assert(m_handler != nullptr && "No handler set");
    m_value = m_handler->clone(other.m_value.get());
  };

 private:
  PolymorphicValue(std::unique_ptr<T> value) : m_value{std::move(value)} {}

  struct HandlerBase {
    virtual std::unique_ptr<T> clone(T* src) const = 0;
    virtual void assign(T* src, T* dest) const = 0;
  };

  template <typename U>
  struct Handler : public HandlerBase {
    std::unique_ptr<T> clone(T* src) const override {
      return std::make_unique<U>(*static_cast<U*>(src));
    }

    void assign(T* src, T* dest) const override {
      U* _src = static_cast<U*>(src);
      U* _dest = static_cast<U*>(dest);
      (*_dest) = (*_src);
    }
  };

  std::unique_ptr<T> m_value{nullptr};
  HandlerBase* m_handler{nullptr};
};

}  // namespace Acts
