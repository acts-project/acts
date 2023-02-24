// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/ThrowAssert.hpp"
#include "ActsExamples/Framework/IAlgorithm.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <iostream>
#include <typeinfo>

namespace ActsExamples {

class DataHandleBase {
 protected:
  virtual ~DataHandleBase() = default;

  DataHandleBase(IAlgorithm* parent, const std::string& name)
      : m_parent(parent), m_name(name) {}

 public:
  const std::string& key() const { return m_key.value(); }

  virtual const std::type_info& typeInfo() const = 0;

  bool isInitialized() const { return m_key.has_value(); }

  const std::string& name() const { return m_name; }

 protected:
  IAlgorithm* m_parent{nullptr};
  std::string m_name;
  std::optional<std::string> m_key{};
};

template <typename T>
class WriteDataHandle final : public DataHandleBase {
 public:
  WriteDataHandle(IAlgorithm* parent, const std::string& name)
      : DataHandleBase{parent, name} {
    m_parent->registerWriteHandle(*this);
  }

  void operator()(const AlgorithmContext& ctx, T&& value) const {
    throw_assert(isInitialized(), "WriteDataHandle not initialized");
    ctx.eventStore.add(m_key.value(), std::move(value));
  }

  void initialize(const std::string& key) {
    throw_assert(!key.empty(), "Input key cannot be empty");
    m_key = key;
  }

  const std::type_info& typeInfo() const override { return typeid(T); };
};

template <typename T>
class ReadDataHandle final : public DataHandleBase {
 public:
  ReadDataHandle(IAlgorithm* parent, const std::string& name)
      : DataHandleBase{parent, name} {
    m_parent->registerReadHandle(*this);
  }

  void initialize(const std::string& key) {
    throw_assert(!key.empty(), "Input key cannot be empty");
    m_key = key;
  }

  const T& operator()(const AlgorithmContext& ctx) const {
    throw_assert(isInitialized(), "ReadDataHandle not initialized");
    return ctx.eventStore.get<T>(m_key.value());
  }

  const std::type_info& typeInfo() const override { return typeid(T); };
};

}  // namespace ActsExamples
