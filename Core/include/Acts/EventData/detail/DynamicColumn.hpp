// This file is part of the Acts project.
//
// Copyright (C) 2022 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <any>
#include <memory>

namespace Acts::detail {

struct DynamicColumnBase {
  virtual ~DynamicColumnBase() = default;

  virtual std::any get(size_t i) = 0;
  virtual std::any get(size_t i) const = 0;

  virtual void add() = 0;
  virtual void clear() = 0;

  virtual std::unique_ptr<DynamicColumnBase> clone() const = 0;
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

  std::unique_ptr<DynamicColumnBase> clone() const override {
    return std::make_unique<DynamicColumn<T>>(*this);
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
  void clear() override { m_vector.clear(); }

  std::unique_ptr<DynamicColumnBase> clone() const override {
    return std::make_unique<DynamicColumn<bool>>(*this);
  }

  std::vector<Wrapper> m_vector;
};

}  // namespace Acts::detail
