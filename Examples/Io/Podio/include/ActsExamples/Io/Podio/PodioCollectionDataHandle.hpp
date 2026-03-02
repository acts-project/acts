// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/DataHandle.hpp"

#include <concepts>
#include <cstdint>
#include <memory>
#include <stdexcept>

#include <podio/CollectionBase.h>

namespace ActsExamples {

namespace detail {
/// Base class for typed Podio collection read handles.
///
/// Extends @ref ReadDataHandleBase with @c collectionTypeHash() for compatibility
/// checking against @ref PodioCollectionWriteHandle. This enables single inheritance
/// for PodioCollectionReadHandle<T>.
class PodioCollectionTypedReadHandle : public ReadDataHandleBase {
 protected:
  using ReadDataHandleBase::ReadDataHandleBase;

 public:
  virtual ~PodioCollectionTypedReadHandle() = default;

  /// Hash of the concrete collection type for compatibility checking.
  /// Used by @ref PodioCollectionWriteHandle::isCompatible to validate that
  /// typed read handles match the declared write handle type.
  virtual std::uint64_t collectionTypeHash() const = 0;
};
}  // namespace detail

/// Type-safe PODIO collection write handle.
///
/// This is the handle intended for algorithm code. It enforces and advertises
/// the concrete collection type for sequencer type checking while storing as
/// `std::unique_ptr<podio::CollectionBase>` in the whiteboard.
///
/// @tparam T The PODIO collection type (e.g. `edm4hep::TrackCollection`).
///           Must inherit from `podio::CollectionBase`.
template <typename T>
  requires std::derived_from<T, podio::CollectionBase>
class PodioCollectionWriteHandle final : public WriteDataHandleBase {
 public:
  /// Construct a write handle.
  /// @param parent The sequence element that owns this handle
  /// @param name The handle name (used for configuration and error messages)
  PodioCollectionWriteHandle(SequenceElement* parent, const std::string& name)
      : WriteDataHandleBase{parent, name} {
    registerAsWriteHandle();
  }

  /// Store a collection in the event store.
  /// @param ctx The algorithm context (uses `ctx.eventStore`)
  /// @param collection The collection to store (moved from function argument)
  void operator()(const AlgorithmContext& ctx, T&& collection) const {
    (*this)(ctx.eventStore, std::move(collection));
  }

  /// Store a collection in the given whiteboard.
  /// @param wb The whiteboard to store in
  /// @param collection The collection to store
  void operator()(WhiteBoard& wb, T&& collection) const {
    store(wb, std::make_unique<T>(std::move(collection)));
  }

  /// Store an existing unique_ptr to a collection.
  /// @param ctx The algorithm context (uses `ctx.eventStore`)
  /// @param collection The collection to store (ownership transferred)
  void operator()(const AlgorithmContext& ctx,
                  std::unique_ptr<T> collection) const {
    (*this)(ctx.eventStore, std::move(collection));
  }

  /// Store an existing unique_ptr to a collection.
  /// @param wb The whiteboard to store in
  /// @param collection The collection to store (ownership transferred)
  void operator()(WhiteBoard& wb, std::unique_ptr<T> collection) const {
    store(wb, std::move(collection));
  }

  /// Return the type information of the collection.
  /// @return The type information of the collection
  const std::type_info& typeInfo() const override {
    return typeid(std::unique_ptr<T>);
  }

  /// Return the type hash of the collection.
  /// @return The type hash of the collection
  std::uint64_t typeHash() const override {
    return Acts::typeHash<std::unique_ptr<T>>();
  }

  /// Check if the collection is compatible with another data handle.
  /// @param other The other data handle to check compatibility with
  /// @return True if the collection is compatible, false otherwise
  bool isCompatible(const DataHandleBase& other) const override {
    const auto baseCollectionHash =
        Acts::typeHash<std::unique_ptr<podio::CollectionBase>>();
    if (other.typeHash() == baseCollectionHash) {
      return true;
    }
    const auto* typedReadHandle =
        dynamic_cast<const detail::PodioCollectionTypedReadHandle*>(&other);
    return typedReadHandle != nullptr &&
           typedReadHandle->collectionTypeHash() == typeHash();
  }

 private:
  /// Internal helper to add a collection to the whiteboard.
  /// @param wb The whiteboard to add the collection to
  /// @param collection The collection to add
  void store(WhiteBoard& wb,
             std::unique_ptr<podio::CollectionBase> collection) const {
    if (!isInitialized()) {
      throw std::runtime_error{"PodioCollectionWriteHandle '" + fullName() +
                               "' not initialized"};
    }
    add(wb, std::move(collection));
  }
};

/// A read handle for concrete Podio collection types.
///
/// This handle reads from whiteboard entries stored as
/// std::unique_ptr<podio::CollectionBase> and performs a checked cast to T.
/// Throws if the stored collection is null or has an incompatible type.
///
/// @tparam T The PODIO collection type to read (e.g. edm4hep::TrackCollection).
///           Must inherit from podio::CollectionBase.
template <typename T>
  requires std::derived_from<T, podio::CollectionBase>
class PodioCollectionReadHandle final
    : public detail::PodioCollectionTypedReadHandle {
 public:
  /// Construct a read handle.
  /// @param parent The sequence element that owns this handle
  /// @param name The handle name (used for configuration and error messages)
  PodioCollectionReadHandle(SequenceElement* parent, const std::string& name)
      : PodioCollectionTypedReadHandle{parent, name} {
    registerAsReadHandle();
  }

  /// Read the collection from the event store.
  /// @param ctx The algorithm context (uses ctx.eventStore)
  /// @return Reference to the stored collection
  /// @throws std::runtime_error if not initialized or collection is null
  /// @throws std::out_of_range if type cast fails
  const T& operator()(const AlgorithmContext& ctx) const {
    return (*this)(ctx.eventStore);
  }

  /// Read the collection from the given whiteboard.
  /// @param wb The whiteboard to read from
  /// @return Reference to the stored collection
  /// @throws std::runtime_error if not initialized or collection is null
  /// @throws std::out_of_range if type cast fails
  const T& operator()(const WhiteBoard& wb) const {
    if (!isInitialized()) {
      throw std::runtime_error{"PodioCollectionReadHandle '" + fullName() +
                               "' not initialized"};
    }
    const auto& collection = get<std::unique_ptr<podio::CollectionBase>>(wb);
    if (collection == nullptr) {
      throw std::runtime_error{"PodioCollectionReadHandle '" + fullName() +
                               "' found null collection"};
    }
    const auto* typedCollection = dynamic_cast<const T*>(collection.get());
    if (typedCollection == nullptr) {
      throw std::out_of_range{"PodioCollectionReadHandle '" + fullName() +
                              "' could not cast collection to requested "
                              "concrete type"};
    }
    return *typedCollection;
  }

  /// Return the type information of the collection.
  /// @return The type information of the collection
  const std::type_info& typeInfo() const override {
    return typeid(std::unique_ptr<T>);
  }

  /// Return the type hash of the collection.
  /// @return The type hash of the collection
  std::uint64_t typeHash() const override {
    return Acts::typeHash<std::unique_ptr<T>>();
  }

  /// Return the type hash of the collection.
  /// @return The type hash of the collection
  std::uint64_t collectionTypeHash() const override { return typeHash(); }
};

}  // namespace ActsExamples
