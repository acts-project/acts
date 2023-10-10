// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <any>
#include <cassert>
#include <memory>
#include <vector>

namespace Acts::detail {

struct DynamicColumnBase {
  virtual ~DynamicColumnBase() = default;

  virtual std::any get(size_t i) = 0;
  virtual std::any get(size_t i) const = 0;

  virtual void add() = 0;
  virtual void clear() = 0;
  virtual void reserve(size_t size) = 0;
  virtual void erase(size_t i) = 0;
  virtual size_t size() const = 0;
  virtual void copyFrom(size_t dstIdx, const DynamicColumnBase& src,
                        size_t srcIdx) = 0;

  virtual std::unique_ptr<DynamicColumnBase> clone(
      bool empty = false) const = 0;
};

template <typename T>
struct DynamicColumn : public DynamicColumnBase {
  std::any get(size_t i) override {
    assert(i < m_vector.size() && "DynamicColumn out of bounds");
    return &m_vector[i];
  }

  std::any get(size_t i) const override {
    assert(i < m_vector.size() && "DynamicColumn out of bounds");
    return &m_vector[i];
  }

  void add() override { m_vector.emplace_back(); }
  void clear() override { m_vector.clear(); }
  void reserve(size_t size) override { m_vector.reserve(size); }
  void erase(size_t i) override { m_vector.erase(m_vector.begin() + i); }
  size_t size() const override { return m_vector.size(); }

  std::unique_ptr<DynamicColumnBase> clone(bool empty) const override {
    if (empty) {
      return std::make_unique<DynamicColumn<T>>();
    }
    return std::make_unique<DynamicColumn<T>>(*this);
  }

  void copyFrom(size_t dstIdx, const DynamicColumnBase& src,
                size_t srcIdx) override {
    const auto* other = dynamic_cast<const DynamicColumn<T>*>(&src);
    assert(other != nullptr &&
           "Source column is not of same type as destination");
    m_vector.at(dstIdx) = other->m_vector.at(srcIdx);
  }

  std::vector<T> m_vector;
};

template <>
struct DynamicColumn<bool> : public DynamicColumnBase {
  struct Wrapper {
    bool value;
  };

  std::any get(size_t i) override {
    assert(i < m_vector.size() && "DynamicColumn out of bounds");
    return &m_vector[i].value;
  }

  std::any get(size_t i) const override {
    assert(i < m_vector.size() && "DynamicColumn out of bounds");
    return &m_vector[i].value;
  }

  void add() override { m_vector.emplace_back(); }
  void reserve(size_t size) override { m_vector.reserve(size); }
  void clear() override { m_vector.clear(); }
  void erase(size_t i) override { m_vector.erase(m_vector.begin() + i); }
  size_t size() const override { return m_vector.size(); }

  std::unique_ptr<DynamicColumnBase> clone(bool empty) const override {
    if (empty) {
      return std::make_unique<DynamicColumn<bool>>();
    }
    return std::make_unique<DynamicColumn<bool>>(*this);
  }

  void copyFrom(size_t dstIdx, const DynamicColumnBase& src,
                size_t srcIdx) override {
    const auto* other = dynamic_cast<const DynamicColumn<bool>*>(&src);
    assert(other != nullptr &&
           "Source column is not of same type as destination");
    m_vector.at(dstIdx) = other->m_vector.at(srcIdx);
  }

  std::vector<Wrapper> m_vector;
};

}  // namespace Acts::detail
