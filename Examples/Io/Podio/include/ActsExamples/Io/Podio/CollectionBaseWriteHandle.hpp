// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/DataHandle.hpp"

#include <cstdint>
#include <memory>
#include <type_traits>

#include <podio/CollectionBase.h>

namespace ActsExamples {

namespace detail {
class CollectionBaseTypedReadHandle {
 public:
  virtual ~CollectionBaseTypedReadHandle() = default;
  virtual std::uint64_t collectionTypeHash() const = 0;
};
}  // namespace detail

/// A write handle for Podio collections that provides type-safe storage and
/// retrieval. This handle manages the lifecycle of Podio collections and
/// ensures they are properly stored in the event store.
class CollectionBaseWriteHandle : public WriteDataHandleBase {
 public:
  /// Construct a new collection write handle
  /// @param parent The sequence element that owns this handle
  /// @param name The name of the handle
  CollectionBaseWriteHandle(SequenceElement* parent, const std::string& name);

  /// Store a collection in the event store
  /// @param ctx The algorithm context containing the event store
  /// @param collection The collection to store
  void operator()(const AlgorithmContext& ctx,
                  std::unique_ptr<podio::CollectionBase> collection) const {
    validateCollectionPointer(collection.get());
    store(ctx.eventStore, std::move(collection));
  }

  /// Store a collection directly in a whiteboard
  /// @param wb The whiteboard to store the collection in
  /// @param collection The collection to store
  void operator()(WhiteBoard& wb,
                  std::unique_ptr<podio::CollectionBase> collection) const {
    validateCollectionPointer(collection.get());
    store(wb, std::move(collection));
  }

  /// Store a derived Podio collection in the event store
  /// @param ctx The algorithm context containing the event store
  /// @param collection The collection to store
  template <typename T>
  void operator()(const AlgorithmContext& ctx, T collection) const
    requires(!std::is_same_v<std::remove_cvref_t<T>, podio::CollectionBase> &&
             std::is_base_of_v<podio::CollectionBase, std::remove_cvref_t<T>>)
  {
    using Collection = std::remove_cvref_t<T>;
    validateConcreteType<Collection>();
    store(ctx.eventStore, std::make_unique<Collection>(std::move(collection)));
  }

  /// Store a derived Podio collection directly in a whiteboard
  /// @param wb The whiteboard to store the collection in
  /// @param collection The collection to store
  template <typename T>
  void operator()(WhiteBoard& wb, T collection) const
    requires(!std::is_same_v<std::remove_cvref_t<T>, podio::CollectionBase> &&
             std::is_base_of_v<podio::CollectionBase, std::remove_cvref_t<T>>)
  {
    using Collection = std::remove_cvref_t<T>;
    validateConcreteType<Collection>();
    store(wb, std::make_unique<Collection>(std::move(collection)));
  }

  /// Declare the concrete collection type for sequencer type checking.
  ///
  /// Data is still stored in the whiteboard as
  /// std::unique_ptr<podio::CollectionBase>, but reads through
  /// CollectionBaseReadHandle<T> can be validated against this concrete type.
  template <typename T>
  void declareConcreteType()
    requires std::is_base_of_v<podio::CollectionBase, std::remove_cvref_t<T>>
  {
    using Collection = std::remove_cvref_t<T>;
    m_declaredTypeInfo = &typeid(std::unique_ptr<Collection>);
    m_declaredTypeHash = Acts::typeHash<std::unique_ptr<Collection>>();
    if constexpr (std::is_same_v<Collection, podio::CollectionBase>) {
      m_typeValidator = nullptr;
    } else {
      m_typeValidator = [](const podio::CollectionBase* collection) {
        return collection == nullptr ||
               dynamic_cast<const Collection*>(collection) != nullptr;
      };
    }
  }

  /// Internal method to store a collection in a whiteboard
  /// @param wb The whiteboard to store the collection in
  /// @param collection The collection to store
  void store(WhiteBoard& wb,
             std::unique_ptr<podio::CollectionBase> collection) const;

  bool isCompatible(const DataHandleBase& other) const override;

  /// Get the type info for this handle
  /// @return The currently declared concrete collection type
  const std::type_info& typeInfo() const override;

  std::uint64_t typeHash() const override;

 private:
  template <typename T>
  void validateConcreteType() const {
    static_assert(std::is_base_of_v<podio::CollectionBase, T>);
    if (m_declaredTypeHash == Acts::typeHash<std::unique_ptr<podio::CollectionBase>>()) {
      return;
    }
    if (m_declaredTypeHash != Acts::typeHash<std::unique_ptr<T>>()) {
      throw std::invalid_argument{
          "CollectionBaseWriteHandle '" + fullName() +
          "' received a collection with a different concrete type than "
          "declared"};
    }
  }

  void validateCollectionPointer(const podio::CollectionBase* collection) const {
    if (m_typeValidator != nullptr && !m_typeValidator(collection)) {
      throw std::invalid_argument{
          "CollectionBaseWriteHandle '" + fullName() +
          "' received a collection with a different concrete type than "
          "declared"};
    }
  }

  const std::type_info* m_declaredTypeInfo =
      &typeid(std::unique_ptr<podio::CollectionBase>);
  std::uint64_t m_declaredTypeHash =
      Acts::typeHash<std::unique_ptr<podio::CollectionBase>>();
  using TypeValidator = bool (*)(const podio::CollectionBase*);
  TypeValidator m_typeValidator = nullptr;
};

/// A read handle for concrete Podio collection types.
///
/// This handle reads from whiteboard entries stored as
/// std::unique_ptr<podio::CollectionBase> and performs a checked cast to T.
template <typename T>
class CollectionBaseReadHandle final : public ReadDataHandleBase,
                                       public detail::CollectionBaseTypedReadHandle {
 public:
  CollectionBaseReadHandle(SequenceElement* parent, const std::string& name)
      : ReadDataHandleBase{parent, name} {
    static_assert(std::is_base_of_v<podio::CollectionBase, T>);
    registerAsReadHandle();
  }

  const T& operator()(const AlgorithmContext& ctx) const {
    return (*this)(ctx.eventStore);
  }

  const T& operator()(const WhiteBoard& wb) const {
    if (!isInitialized()) {
      throw std::runtime_error{"CollectionBaseReadHandle '" + fullName() +
                               "' not initialized"};
    }
    const auto& collection = get<std::unique_ptr<podio::CollectionBase>>(wb);
    if (collection == nullptr) {
      throw std::runtime_error{"CollectionBaseReadHandle '" + fullName() +
                               "' found null collection"};
    }
    const auto* typedCollection = dynamic_cast<const T*>(collection.get());
    if (typedCollection == nullptr) {
      throw std::out_of_range{"CollectionBaseReadHandle '" + fullName() +
                              "' could not cast collection to requested "
                              "concrete type"};
    }
    return *typedCollection;
  }

  const std::type_info& typeInfo() const override {
    return typeid(std::unique_ptr<T>);
  }

  std::uint64_t typeHash() const override {
    return Acts::typeHash<std::unique_ptr<T>>();
  }

  std::uint64_t collectionTypeHash() const override { return typeHash(); }
};

}  // namespace ActsExamples
