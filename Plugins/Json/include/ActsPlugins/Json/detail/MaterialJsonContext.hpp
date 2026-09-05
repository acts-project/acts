// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/MaterialSlab.hpp"

#include <cstddef>
#include <memory>
#include <stdexcept>
#include <string>
#include <vector>

namespace Acts::detail {

/// A store of material slabs that several index grids can share
using MaterialSlabStore = std::shared_ptr<std::vector<MaterialSlab>>;

/// Context threaded through the material encoders of a single document.
///
/// Material grids that index into a shared slab store register that store
/// here instead of inlining it, so the document holds one copy per store.
/// A default constructed context has no store table, which makes the
/// encoders fall back to inlining and keeps standalone payloads
/// self-contained.
class MaterialJsonEncodeContext {
 public:
  /// Default construction, the store table is disabled
  MaterialJsonEncodeContext() = default;

  /// Create a context that collects slab stores in a table
  ///
  /// @return a context with the store table enabled
  static MaterialJsonEncodeContext withStoreTable() {
    MaterialJsonEncodeContext ctx;
    ctx.m_storeTable = true;
    return ctx;
  }

  /// @return whether encoders should reference the store table
  bool storeTableEnabled() const { return m_storeTable; }

  /// Look up, and if needed assign, the table id of a slab store
  ///
  /// Stores are keyed on pointer identity, never on content, so two stores
  /// that happen to hold the same slabs stay distinct.
  ///
  /// @param store the slab store to register
  ///
  /// @return the id of the store in the table
  std::size_t storeId(const MaterialSlabStore& store) {
    if (!m_storeTable) {
      throw std::logic_error(
          "MaterialJsonEncodeContext: store table is not enabled");
    }
    if (store == nullptr) {
      throw std::invalid_argument(
          "MaterialJsonEncodeContext: cannot register a null slab store");
    }
    for (std::size_t is = 0; is < m_stores.size(); ++is) {
      if (m_stores[is] == store) {
        return is;
      }
    }
    m_stores.push_back(store);
    return m_stores.size() - 1;
  }

  /// @return the collected slab stores, in encounter order
  const std::vector<MaterialSlabStore>& stores() const { return m_stores; }

 private:
  bool m_storeTable = false;
  std::vector<MaterialSlabStore> m_stores;
};

/// Context threaded through the material decoders of a single document.
///
/// It carries the slab store table read from the document so that all
/// grids referencing the same id end up sharing one allocation.
class MaterialJsonDecodeContext {
 public:
  /// Default construction, no store table is available
  MaterialJsonDecodeContext() = default;

  /// Install the slab store table of the document
  ///
  /// @param stores the stores, indexed by their document id
  void setStores(std::vector<MaterialSlabStore> stores) {
    m_stores = std::move(stores);
    m_storeTable = true;
  }

  /// @return whether a store table is available
  bool storeTableEnabled() const { return m_storeTable; }

  /// Resolve a slab store by its document id
  ///
  /// @param id the id of the store in the table
  ///
  /// @return the shared slab store
  MaterialSlabStore store(std::size_t id) const {
    if (!m_storeTable) {
      throw std::invalid_argument(
          "MaterialJsonDecodeContext: the payload references slab store " +
          std::to_string(id) + " but no store table was read");
    }
    if (id >= m_stores.size()) {
      throw std::invalid_argument("MaterialJsonDecodeContext: slab store id " +
                                  std::to_string(id) +
                                  " is out of range, the table holds " +
                                  std::to_string(m_stores.size()) + " stores");
    }
    return m_stores[id];
  }

  /// @return the number of stores in the table
  std::size_t size() const { return m_stores.size(); }

 private:
  bool m_storeTable = false;
  std::vector<MaterialSlabStore> m_stores;
};

}  // namespace Acts::detail
