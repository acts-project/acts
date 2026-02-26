// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/HashedString.hpp"
#include "Acts/Utilities/TypeTag.hpp"

#include <functional>
#include <unordered_map>

#include <pybind11/pybind11.h>

namespace ActsPython {

/// The WhiteBoard is an event-store container that holds arbitrary C++ objects
/// by name. Python algorithms need to read these objects through pybind11, but
/// the WhiteBoard stores them in a type-erased form (`void*`). This registry
/// bridges that gap by mapping Python types to their C++ counterparts and
/// providing a downcast function that converts the stored pointer into a
/// properly reference-managed pybind11 object.
///
/// Usage:
/// 1. When defining pybind11 bindings for a type T that will be stored on the
///    WhiteBoard, call `WhiteBoardRegistry::registerClass(pyClass)` or
///    `WhiteBoardRegistry::registerType<T>(pyType)` immediately after the
///    py::class_<T> definition.
/// 2. When a Python algorithm creates a `ReadDataHandle` for that type, the
///    registry is consulted to find the type info and downcast function for
///    safe retrieval from the `WhiteBoard`.
class WhiteBoardRegistry {
 public:
  /// Function that converts a type-erased pointer from the WhiteBoard into a
  /// pybind11 object. The wbPy argument is used for reference_internal
  /// lifetime.
  using DowncastFunction = std::function<pybind11::object(
      const void* data, const pybind11::object& wbPy)>;

  /// Register a pybind11-bound type T for WhiteBoard read access.
  /// Call this after the `py::class_<T>` definition.
  /// @tparam Ts The types to register.
  /// @param pyClass The pybind11 class object to register.
  template <typename... Ts>
  static void registerClass(const pybind11::class_<Ts...>& pyClass) {
    namespace py = pybind11;
    using type = pybind11::class_<Ts...>::type;
    registerType<type>(pyClass);
  }

  /// Register a C++ type `~T` with its pybind11 Python type for WhiteBoard
  /// access. Use when the `py::class_<T>` type cannot be deduced (e.g. for
  /// template types).
  /// @tparam T The C++ type to register.
  /// @param pyType The pybind11 Python type object to register.
  template <typename T>
  static void registerType(const pybind11::object& pyType) {
    namespace py = pybind11;

    using type = T;

    instance()[pyType.ptr()] = {
        .fn = [](const void* ptr, const py::object& wbPy) -> py::object {
          // wb needed to ensure correct lifetime
          return py::cast(*static_cast<const type*>(ptr),
                          py::return_value_policy::reference_internal, wbPy);
        },
        .typeinfo = &typeid(type),
        .typeHash = Acts::typeHash<type>(),
    };
  }

  /// Per-type registry entry: downcast function and type metadata for lookups.
  struct RegistryEntry {
    DowncastFunction fn{
        nullptr};  ///< Converts `void*` + `WhiteBoard` -> `py::object`
    const std::type_info* typeinfo{nullptr};  ///< C++ type for type checking
    std::uint64_t typeHash{0};  ///< Hash for runtime type verification
  };

  /// Look up a registered type by its pybind11 Python type object.
  /// @param pyType The pybind11 Python type object to look up.
  /// @return Pointer to the `RegistryEntry`, or `nullptr` if not registered
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
