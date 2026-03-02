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
#include <stdexcept>
#include <type_traits>

#include <podio/CollectionBase.h>

namespace ActsExamples {

namespace detail {
/// Base class for typed Podio collection read handles.
///
/// Extends ReadDataHandleBase with collectionTypeHash() for compatibility
/// checking against PodioCollectionWriteHandle. This enables single inheritance
/// for PodioCollectionReadHandle<T>.
class PodioCollectionTypedReadHandle : public ReadDataHandleBase {
 protected:
  using ReadDataHandleBase::ReadDataHandleBase;

 public:
  virtual ~PodioCollectionTypedReadHandle() = default;
  virtual std::uint64_t collectionTypeHash() const = 0;
};
}  // namespace detail

/// Type-safe PODIO collection write handle.
///
/// This is the handle intended for algorithm code. It enforces and advertises
/// the concrete collection type for sequencer type checking while storing as
/// std::unique_ptr<podio::CollectionBase> in the whiteboard.
template <typename T>
class PodioCollectionWriteHandle final : public WriteDataHandleBase {
 public:
  PodioCollectionWriteHandle(SequenceElement* parent, const std::string& name)
      : WriteDataHandleBase{parent, name} {
    registerAsWriteHandle();
  }

  void operator()(const AlgorithmContext& ctx, T collection) const {
    (*this)(ctx.eventStore, std::move(collection));
  }

  void operator()(WhiteBoard& wb, T collection) const {
    store(wb, std::make_unique<T>(std::move(collection)));
  }

  void operator()(const AlgorithmContext& ctx,
                  std::unique_ptr<T> collection) const {
    (*this)(ctx.eventStore, std::move(collection));
  }

  void operator()(WhiteBoard& wb, std::unique_ptr<T> collection) const {
    store(wb, std::move(collection));
  }

  const std::type_info& typeInfo() const override {
    return typeid(std::unique_ptr<T>);
  }

  std::uint64_t typeHash() const override {
    return Acts::typeHash<std::unique_ptr<T>>();
  }

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
template <typename T>
class PodioCollectionReadHandle final
    : public detail::PodioCollectionTypedReadHandle {
 public:
  PodioCollectionReadHandle(SequenceElement* parent, const std::string& name)
      : PodioCollectionTypedReadHandle{parent, name} {
    registerAsReadHandle();
  }

  const T& operator()(const AlgorithmContext& ctx) const {
    return (*this)(ctx.eventStore);
  }

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

  const std::type_info& typeInfo() const override {
    return typeid(std::unique_ptr<T>);
  }

  std::uint64_t typeHash() const override {
    return Acts::typeHash<std::unique_ptr<T>>();
  }

  std::uint64_t collectionTypeHash() const override { return typeHash(); }
};

}  // namespace ActsExamples
