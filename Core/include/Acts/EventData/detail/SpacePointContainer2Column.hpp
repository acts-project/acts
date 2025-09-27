// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <memory>
#include <vector>

namespace Acts {

class SpacePointContainer2;
template <typename T, bool read_only>
class SpacePointColumnProxy;
template <typename T>
using MutableSpacePointColumnProxy = SpacePointColumnProxy<T, false>;
template <typename T>
using ConstSpacePointColumnProxy = SpacePointColumnProxy<T, true>;

namespace detail::sp {

// These classes should have gone into `SpacePointContainer2` but a compiler bug
// prevents that

class ColumnHolderBase {
 public:
  virtual ~ColumnHolderBase() = default;

  virtual std::unique_ptr<ColumnHolderBase> copy() const = 0;

  virtual std::size_t size() const = 0;
  virtual void reserve(std::size_t size) = 0;
  virtual void resize(std::size_t size) = 0;
  virtual void clear() = 0;
  virtual void emplace_back() = 0;
};

template <typename T>
class ColumnHolder final : public ColumnHolderBase {
 public:
  using Value = T;
  using Container = std::vector<Value>;
  using MutableProxy = MutableSpacePointColumnProxy<Value>;
  using ConstProxy = ConstSpacePointColumnProxy<Value>;

  ColumnHolder() = default;
  explicit ColumnHolder(Value defaultValue)
      : m_default(std::move(defaultValue)) {}

  MutableProxy proxy(SpacePointContainer2 &container) {
    return MutableProxy(container, m_data);
  }
  ConstProxy proxy(const SpacePointContainer2 &container) const {
    return ConstProxy(container, m_data);
  }

  std::unique_ptr<ColumnHolderBase> copy() const override {
    return std::make_unique<detail::sp::ColumnHolder<T>>(*this);
  }

  std::size_t size() const override { return m_data.size(); }
  void reserve(std::size_t size) override { m_data.reserve(size); }
  void clear() override { m_data.clear(); }
  void resize(std::size_t size) override { m_data.resize(size, m_default); }
  void emplace_back() override { m_data.emplace_back(m_default); }

 private:
  Value m_default{};
  Container m_data;
};

}  // namespace detail::sp
}  // namespace Acts
