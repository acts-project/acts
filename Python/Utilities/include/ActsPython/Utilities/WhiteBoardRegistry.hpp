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

#include <concepts>
#include <cstdint>
#include <functional>
#include <memory>
#include <typeinfo>

#include <pybind11/pybind11.h>

namespace ActsPython {

/// Concept: pybind11 class must use smart_holder for WhiteBoard ownership.
template <typename... Ts>
concept PyClassWithSmartHolder =
    std::same_as<typename pybind11::class_<Ts...>::holder_type,
                 pybind11::smart_holder>;

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
 private:
  /// Function that converts a type-erased pointer from the WhiteBoard into a
  /// pybind11 object. The wbPy argument is used for reference_internal
  /// lifetime.
  using ToPythonFunction = std::function<pybind11::object(
      const Acts::AnyMoveOnly& any, const pybind11::object& wbPy)>;

  /// Function that converts a pybind11 object into a type-erased pointer for
  /// the WhiteBoard.
  using FromPythonFunction = std::function<std::unique_ptr<Acts::AnyMoveOnly>(
      const pybind11::object& obj)>;

 public:
  /// Opaque handle to an internal registry entry.
  struct EntryHandle {
   protected:
    EntryHandle() = default;
    ~EntryHandle() = default;
  };

 private:
  static void registerTypeImpl(PyObject* pyType, ToPythonFunction toPython,
                               FromPythonFunction fromPython,
                               const std::type_info* typeinfo,
                               std::uint64_t typeHash);

 public:
  /// Register a C++ type T with its pybind11 Python type for WhiteBoard
  /// access. Use when the `py::class_<T>` type cannot be deduced (e.g. for
  /// template types).
  /// @tparam T The C++ type to register.
  template <typename T>
  static void registerType() {
    namespace py = pybind11;
    return registerType<T>(py::type::of<T>());
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

    registerTypeImpl(
        pyType.ptr(),
        [](const Acts::AnyMoveOnly& any, const py::object& wbPy) -> py::object {
          // wb needed to ensure correct lifetime
          return py::cast(any.as<type>(),
                          py::return_value_policy::reference_internal, wbPy);
        },
        [](const py::object& obj) {
          // This communicates to pybind11's smart_holder that the object is
          // consumed in C++
          auto up = py::cast<std::unique_ptr<T>>(obj);
          return std::make_unique<Acts::AnyMoveOnly>(std::move(*up));
        },
        &typeid(type), Acts::typeHash<type>());
  }

  /// Register a pybind11-bound type T for WhiteBoard read access.
  /// Call this after the `py::class_<T>` definition.
  /// @tparam Ts The types to register.
  /// @param pyClass The pybind11 class object to register.
  template <typename... Ts>
  static void registerClass(const pybind11::class_<Ts...>& pyClass)
    requires PyClassWithSmartHolder<Ts...>
  {
    using type = pybind11::class_<Ts...>::type;
    registerType<type>(pyClass);
  }

  /// Look up a registered type by Python type object.
  /// @param pyType The Python type object to look up.
  /// @return Opaque handle to the internal entry, or nullptr if not registered.
  static const EntryHandle* find(PyObject* pyType) noexcept;

  /// Access the registered C++ type metadata.
  static const std::type_info* typeInfo(const EntryHandle* entry) noexcept;

  /// Access the registered hash for runtime type verification.
  static std::uint64_t typeHash(const EntryHandle* entry) noexcept;

  /// Convert a WhiteBoard object into a Python object.
  /// @return New reference (`PyObject*`) owned by the caller.
  static PyObject* toPython(const EntryHandle* entry,
                            const Acts::AnyMoveOnly& any, PyObject* wbPy);

  /// Convert a Python object into a WhiteBoard-storable holder.
  static std::unique_ptr<Acts::AnyMoveOnly> fromPython(const EntryHandle* entry,
                                                       PyObject* obj);

 private:
  WhiteBoardRegistry() = default;
};

}  // namespace ActsPython
