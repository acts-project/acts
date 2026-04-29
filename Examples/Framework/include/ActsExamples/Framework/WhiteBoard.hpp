// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Any.hpp"
#include "Acts/Utilities/HashedString.hpp"
#include "Acts/Utilities/Logger.hpp"

#include <algorithm>
#include <cstddef>
#include <memory>
#include <sstream>
#include <stdexcept>
#include <string>
#include <string_view>
#include <typeinfo>
#include <unordered_map>
#include <utility>
#include <vector>

namespace ActsExamples {

/// A container to store arbitrary objects with ownership transfer.
///
/// This is an append-only container that takes ownership of the objects
/// added to it. Once an object has been added, it can only be read but not
/// be modified. Trying to replace an existing object is considered an error.
/// Its lifetime is bound to the lifetime of the white board.
class WhiteBoard {
 public:
  struct StringHash {
    using is_transparent = void;  // Enables heterogeneous operations.

    std::size_t operator()(std::string_view sv) const {
      std::hash<std::string_view> hasher;
      return hasher(sv);
    }
  };

  using StoreValue =
      std::pair<std::shared_ptr<Acts::AnyMoveOnly>, std::uint64_t>;
  using StoreMapType =
      std::unordered_map<std::string, StoreValue, StringHash, std::equal_to<>>;
  using AliasMapType = std::unordered_multimap<std::string, std::string,
                                               StringHash, std::equal_to<>>;

  explicit WhiteBoard(std::unique_ptr<const Acts::Logger> logger =
                          Acts::getDefaultLogger("WhiteBoard",
                                                 Acts::Logging::INFO),
                      AliasMapType objectAliases = {});

  WhiteBoard(const WhiteBoard& other) = delete;
  WhiteBoard& operator=(const WhiteBoard&) = delete;

  WhiteBoard(WhiteBoard&& other) = default;
  WhiteBoard& operator=(WhiteBoard&& other) = default;

  bool exists(const std::string& name) const;

  /// Copies key from another whiteboard to this whiteboard.
  /// This is a low overhead operation, since the data holders are
  /// shared pointers.
  /// Throws an exception if this whiteboard already contains one of
  /// the keys in the other whiteboard.
  void copyFrom(const WhiteBoard& other);

  std::vector<std::string> getKeys() const;

 private:
  /// Find similar names for suggestions with levenshtein-distance
  std::vector<std::string_view> similarNames(const std::string_view& name,
                                             int distThreshold,
                                             std::size_t maxNumber) const;

  /// Store a value on the white board.
  ///
  /// @param name Non-empty identifier to store it under
  /// @param holder The value to store (shared ownership for copyFrom)
  /// @param typeHash Hash of the stored type for runtime verification
  /// @throws std::invalid_argument on empty or duplicate name
  void addHolder(const std::string& name,
                 const std::shared_ptr<Acts::AnyMoveOnly>& holder,
                 std::uint64_t typeHash);

  /// Store a value on the white board (from unique ownership).
  void addHolder(const std::string& name,
                 std::unique_ptr<Acts::AnyMoveOnly> holder,
                 std::uint64_t typeHash);

  /// Store an object on the white board and transfer ownership.
  ///
  /// @param name Non-empty identifier to store it under
  /// @param object Movable reference to the transferable object
  template <typename T>
  const T& add(const std::string& name, T&& object) {
    addHolder(name,
              std::make_shared<Acts::AnyMoveOnly>(std::forward<T>(object)),
              Acts::typeHash<T>());
    return get<T>(name);
  }

  /// Get access to a stored object.
  ///
  /// @param[in] name Identifier for the object
  /// @return reference to the stored object
  /// @throws std::out_of_range if no object is stored under the requested name
  template <typename T>
  const T& get(const std::string& name) const;

  template <typename T>
  Acts::AnyMoveOnly* getHolder(const std::string& name) const;

  /// Returns (pointer to stored value, type hash). Throws if not found.
  std::pair<Acts::AnyMoveOnly*, std::uint64_t> getHolder(
      const std::string& name) const;

  template <typename T>
  T pop(const std::string& name);

  std::unique_ptr<const Acts::Logger> m_logger;

  StoreMapType m_store;

  AliasMapType m_objectAliases;

  const Acts::Logger& logger() const { return *m_logger; }

  static std::string typeMismatchMessage(const std::string& name,
                                         const char* req, const char* act);

  friend class DataHandleBase;
};

inline WhiteBoard::WhiteBoard(std::unique_ptr<const Acts::Logger> logger,
                              AliasMapType objectAliases)
    : m_logger(std::move(logger)), m_objectAliases(std::move(objectAliases)) {}

template <typename T>
Acts::AnyMoveOnly* WhiteBoard::getHolder(const std::string& name) const {
  auto it = m_store.find(name);
  if (it == m_store.end()) {
    const auto names = similarNames(name, 10, 3);

    std::stringstream ss;
    if (!names.empty()) {
      ss << ", similar ones are: [ ";
      for (std::size_t i = 0; i < std::min(3ul, names.size()); ++i) {
        ss << "'" << names[i] << "' ";
      }
      ss << "]";
    }

    throw std::out_of_range("Object '" + name + "' does not exists" + ss.str());
  }

  auto& [holder, storedTypeHash] = it->second;
  if (storedTypeHash != Acts::typeHash<T>()) {
    const char* holderTypeName =
        holder->typeInfo() ? holder->typeInfo()->name() : "unknown";
    std::string msg =
        typeMismatchMessage(name, typeid(T).name(), holderTypeName);
    throw std::out_of_range(msg.c_str());
  }

  return holder.get();
}

template <typename T>
inline const T& WhiteBoard::get(const std::string& name) const {
  ACTS_VERBOSE("Attempt to get object '" << name << "' of type "
                                         << typeid(T).name());
  ACTS_VERBOSE("Retrieved object '" << name << "'");
  auto* holder = getHolder<T>(name);
  return holder->template as<T>();
}

template <typename T>
T WhiteBoard::pop(const std::string& name) {
  ACTS_VERBOSE("Pop object '" << name << "'");
  (void)getHolder<T>(name);  // validates type and existence
  auto node = m_store.extract(name);
  return node.mapped().first->template take<T>();
}

inline bool WhiteBoard::exists(const std::string& name) const {
  // TODO remove this function?
  return m_store.contains(name);
}

}  // namespace ActsExamples
