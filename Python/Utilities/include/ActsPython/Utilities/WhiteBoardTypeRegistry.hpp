// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/HashedString.hpp"
#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <functional>
#include <typeindex>
#include <unordered_map>

#include <pybind11/pybind11.h>

namespace ActsPython {

class WhiteBoardRegistry {
 public:
  using DowncastFunction = std::function<pybind11::object(
      const void* data, const pybind11::object& wbPy)>;

  /// Register a pybind11-bound type T for WhiteBoard read access.
  /// Call this immediately after the py::class_<T> definition.
  template <typename... Ts>
  static void registerType(const pybind11::class_<Ts...>& pyType) {
    namespace py = pybind11;

    using type = pybind11::class_<Ts...>::type;

    instance()[pyType.ptr()] = {
        .fn = [](const void* ptr, const py::object& wbPy) -> py::object {
          // wb py seems to be needed to ensure correct lifetime
          return py::cast(*static_cast<const type*>(ptr),
                          py::return_value_policy::reference_internal, wbPy);
        },
        .typeinfo = &typeid(type),
        .typeHash = Acts::typeHash<type>(),
    };
  }

  struct RegistryEntry {
    DowncastFunction fn{nullptr};
    const std::type_info* typeinfo{nullptr};
    std::uint64_t typeHash{0};
  };

  static RegistryEntry* find(const pybind11::object& pyType) {
    if (auto it = instance().find(pyType.ptr()); it != instance().end()) {
      return &it->second;
    }
    return nullptr;
  }

 private:
  WhiteBoardRegistry() = default;

  static inline std::unordered_map<PyObject*, RegistryEntry>& instance() {
    static std::unordered_map<PyObject*, RegistryEntry> map;
    return map;
  }
};

}  // namespace ActsPython
