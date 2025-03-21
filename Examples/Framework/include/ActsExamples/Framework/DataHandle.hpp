// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/SequenceElement.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <stdexcept>
#include <typeinfo>
#include <unordered_map>

namespace Acts {
class Logger;
}

namespace ActsExamples {

/// Base class for all data handles.
///
/// Provides common functionality for tracking the parent sequence element
/// and key name. The key is optional until explicitly initialized.
class DataHandleBase {
 private:
  struct StringHash {
    using is_transparent = void;  // Enables heterogeneous operations.

    std::size_t operator()(std::string_view sv) const {
      std::hash<std::string_view> hasher;
      return hasher(sv);
    }
  };

 protected:
  DataHandleBase(SequenceElement* parent, const std::string& name)
      : m_parent(parent), m_name(name) {}

  // We can't change addresses after construction
  DataHandleBase(const DataHandleBase&) = delete;
  DataHandleBase(DataHandleBase&&) = default;

 public:
  virtual ~DataHandleBase() = default;

  const std::string& key() const { return m_key.value(); }

  virtual const std::type_info& typeInfo() const = 0;

  bool isInitialized() const { return m_key.has_value(); }

  const std::string& name() const { return m_name; }

  void maybeInitialize(std::string_view key);

  virtual bool isCompatible(const DataHandleBase& other) const = 0;

  using StateMapType = std::unordered_map<std::string, const DataHandleBase*,
                                          StringHash, std::equal_to<>>;

  virtual void emulate(StateMapType& state, WhiteBoard::AliasMapType& aliases,
                       const Acts::Logger& logger) const = 0;

  std::string fullName() const { return m_parent->name() + "." + name(); }

 protected:
  void registerAsWriteHandle();
  void registerAsReadHandle();

  // Trampoline functions to avoid having the WhiteBoard as a friend
  template <typename T>
  void add(WhiteBoard& wb, T&& object) const {
    wb.add(m_key.value(), std::forward<T>(object));
  }

  template <typename T>
  const T& get(const WhiteBoard& wb) const {
    return wb.get<T>(m_key.value());
  }

  template <typename T>
  T pop(WhiteBoard& wb) const {
    return wb.pop<T>(m_key.value());
  }

  SequenceElement* m_parent{nullptr};
  std::string m_name;
  std::optional<std::string> m_key{};
};

/// Base class for write data handles.
///
/// Write handles are used to store data in the WhiteBoard. They ensure that:
/// - Each key can only be written once
/// - The key must be non-empty
/// - The data type is consistent for each key
class WriteDataHandleBase : public DataHandleBase {
 protected:
  WriteDataHandleBase(SequenceElement* parent, const std::string& name)
      : DataHandleBase{parent, name} {}

 public:
  void initialize(std::string_view key);

  bool isCompatible(const DataHandleBase& other) const final;

  void emulate(StateMapType& state, WhiteBoard::AliasMapType& aliases,
               const Acts::Logger& logger) const final;
};

/// Base class for read data handles.
///
/// Read handles are used to access data from the WhiteBoard. They ensure that:
/// - The data exists before reading
/// - The data type matches the expected type
/// - The data can be read multiple times
class ReadDataHandleBase : public DataHandleBase {
 protected:
  using DataHandleBase::DataHandleBase;

 public:
  void initialize(std::string_view key);

  bool isCompatible(const DataHandleBase& other) const final;

  void emulate(StateMapType& state, WhiteBoard::AliasMapType& aliases,
               const Acts::Logger& logger) const override;
};

/// Base class for consume data handles.
///
/// Consume handles are used to take ownership of data from the WhiteBoard.
/// They ensure that:
/// - The data exists before consuming
/// - The data type matches the expected type
/// - The data can only be consumed once
/// - The data is removed from the WhiteBoard after consumption
class ConsumeDataHandleBase : public ReadDataHandleBase {
 protected:
  using ReadDataHandleBase::ReadDataHandleBase;

 public:
  void emulate(StateMapType& state, WhiteBoard::AliasMapType& aliases,
               const Acts::Logger& logger) const override;
};

/// A write handle for storing data in the WhiteBoard.
///
/// @tparam T The type of data to store
///
/// Example usage:
/// @code
/// WriteDataHandle<int> handle(parent, "my_data");
/// handle.initialize("my_key");
/// handle(wb, 42);  // Store value
/// @endcode
template <typename T>
class WriteDataHandle final : public WriteDataHandleBase {
 public:
  WriteDataHandle(SequenceElement* parent, const std::string& name)
      : WriteDataHandleBase{parent, name} {
    registerAsWriteHandle();
  }

  void operator()(const AlgorithmContext& ctx, T&& value) const {
    (*this)(ctx.eventStore, std::move(value));
  }

  void operator()(WhiteBoard& wb, T&& value) const {
    if (!isInitialized()) {
      throw std::runtime_error{"WriteDataHandle '" + fullName() +
                               "' not initialized"};
    }
    add(wb, std::move(value));
  }

  const std::type_info& typeInfo() const override { return typeid(T); };
};

/// A read handle for accessing data from the WhiteBoard.
///
/// @tparam T The type of data to read
///
/// Example usage:
/// @code
/// ReadDataHandle<int> handle(parent, "my_data");
/// handle.initialize("my_key");
/// const auto& value = handle(wb);  // Access value
/// @endcode
template <typename T>
class ReadDataHandle final : public ReadDataHandleBase {
 public:
  ReadDataHandle(SequenceElement* parent, const std::string& name)
      : ReadDataHandleBase{parent, name} {
    registerAsReadHandle();
  }

  const T& operator()(const AlgorithmContext& ctx) const {
    return (*this)(ctx.eventStore);
  }

  const T& operator()(const WhiteBoard& wb) const {
    if (!isInitialized()) {
      throw std::runtime_error{"ReadDataHandle '" + fullName() +
                               "' not initialized"};
    }
    return get<T>(wb);
  }

  const std::type_info& typeInfo() const override { return typeid(T); };
};

/// A consume handle for taking ownership of data from the WhiteBoard.
///
/// @tparam T The type of data to consume
///
/// Example usage:
/// @code
/// ConsumeDataHandle<int> handle(parent, "my_data");
/// handle.initialize("my_key");
/// auto value = handle(wb);  // Take ownership of value
/// // value is no longer in WhiteBoard
/// @endcode
template <typename T>
class ConsumeDataHandle final : public ConsumeDataHandleBase {
 public:
  ConsumeDataHandle(SequenceElement* parent, const std::string& name)
      : ConsumeDataHandleBase{parent, name} {
    registerAsReadHandle();
  }

  T operator()(const AlgorithmContext& ctx) const {
    return (*this)(ctx.eventStore);
  }

  T operator()(WhiteBoard& wb) const {
    if (!isInitialized()) {
      throw std::runtime_error{"ConsumeDataHandle '" + fullName() +
                               "' not initialized"};
    }
    return pop<T>(wb);
  }

  const std::type_info& typeInfo() const override { return typeid(T); };
};

}  // namespace ActsExamples
