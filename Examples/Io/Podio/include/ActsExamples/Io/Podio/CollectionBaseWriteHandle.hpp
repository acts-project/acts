// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "ActsExamples/Framework/DataHandle.hpp"

#include <memory>

namespace podio {
class CollectionBase;
}

namespace ActsExamples {

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
    store(ctx.eventStore, std::move(collection));
  }

  /// Store a collection directly in a whiteboard
  /// @param wb The whiteboard to store the collection in
  /// @param collection The collection to store
  void operator()(WhiteBoard& wb,
                  std::unique_ptr<podio::CollectionBase> collection) const {
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
    store(ctx.eventStore, std::make_unique<T>(std::move(collection)));
  }

  /// Store a derived Podio collection directly in a whiteboard
  /// @param wb The whiteboard to store the collection in
  /// @param collection The collection to store
  template <typename T>
  void operator()(WhiteBoard& wb, T collection) const
    requires(!std::is_same_v<std::remove_cvref_t<T>, podio::CollectionBase> &&
             std::is_base_of_v<podio::CollectionBase, std::remove_cvref_t<T>>)
  {
    store(wb, std::make_unique<T>(std::move(collection)));
  }

  /// Internal method to store a collection in a whiteboard
  /// @param wb The whiteboard to store the collection in
  /// @param collection The collection to store
  void store(WhiteBoard& wb,
             std::unique_ptr<podio::CollectionBase> collection) const;

  /// Get the type info for this handle
  /// @return The type info for std::unique_ptr<podio::CollectionBase>
  const std::type_info& typeInfo() const override;
};

}  // namespace ActsExamples
