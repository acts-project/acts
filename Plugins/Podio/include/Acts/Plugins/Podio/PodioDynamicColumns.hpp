// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <any>
#include <string>

#include <podio/Frame.h>
#include <podio/UserDataCollection.h>

namespace Acts::podio_detail {

struct ConstDynamicColumnBase {
  ConstDynamicColumnBase(const std::string& name) : m_name{name} {}

  virtual ~ConstDynamicColumnBase() = default;

  virtual std::any get(size_t i) const = 0;

  virtual size_t size() const = 0;

 protected:
  std::string m_name;
};

struct DynamicColumnBase : public ConstDynamicColumnBase {
  DynamicColumnBase(const std::string& name) : ConstDynamicColumnBase{name} {}

  virtual std::any get(size_t i) = 0;
  std::any get(size_t i) const override = 0;

  virtual void add() = 0;
  virtual void clear() = 0;
  virtual void erase(size_t i) = 0;
  virtual void copyFrom(size_t dstIdx, const DynamicColumnBase& src,
                        size_t srcIdx) = 0;

  virtual std::unique_ptr<DynamicColumnBase> clone(
      bool empty = false) const = 0;

  virtual void releaseInto(podio::Frame& frame, const std::string& prefix) = 0;
};

template <typename T>
struct DynamicColumn : public DynamicColumnBase {
  DynamicColumn(const std::string& name,
                podio::UserDataCollection<T> collection = {})
      : DynamicColumnBase(name), m_collection{std::move(collection)} {}

  std::any get(size_t i) override { return &m_collection.vec().at(i); }

  std::any get(size_t i) const override { return &m_collection.vec().at(i); }

  void add() override { m_collection.vec().emplace_back(); }
  void clear() override { m_collection.clear(); }
  void erase(size_t i) override {
    m_collection.vec().erase(m_collection.vec().begin() + i);
  }
  size_t size() const override { return m_collection.size(); }

  std::unique_ptr<DynamicColumnBase> clone(bool empty) const override {
    if (empty) {
      return std::make_unique<DynamicColumn<T>>(m_name);
    }
    podio::UserDataCollection<T> copy;
    copy.vec().reserve(m_collection.size());
    for (const T& v : m_collection) {
      copy.push_back(v);
    }
    return std::make_unique<DynamicColumn<T>>(m_name, std::move(copy));
  }

  void copyFrom(size_t dstIdx, const DynamicColumnBase& src,
                size_t srcIdx) override {
    const auto* other = dynamic_cast<const DynamicColumn<T>*>(&src);
    assert(other != nullptr &&
           "Source column is not of same type as destination");
    m_collection.vec().at(dstIdx) = other->m_collection.vec().at(srcIdx);
  }

  void releaseInto(podio::Frame& frame, const std::string& prefix) override {
    frame.put(std::move(m_collection), prefix + m_name);
  }

  podio::UserDataCollection<T> m_collection;
};

template <typename T>
struct ConstDynamicColumn : public ConstDynamicColumnBase {
  ConstDynamicColumn(const std::string& name,
                     const podio::UserDataCollection<T>& collection)
      : ConstDynamicColumnBase(name), m_collection{collection} {}

  std::any get(size_t i) const override { return &m_collection.vec().at(i); }
  size_t size() const override { return m_collection.size(); }

  const podio::UserDataCollection<T>& m_collection;
};
}  // namespace Acts::podio_detail
