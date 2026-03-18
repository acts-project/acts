// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsPython/Utilities/WhiteBoardRegistry.hpp"

#include <unordered_map>

namespace py = pybind11;

namespace ActsPython {

namespace {

struct RegistryEntry : WhiteBoardRegistry::EntryHandle {
  std::function<py::object(const Acts::AnyMoveOnly& any,
                           const py::object& wbPy)>
      toPython{nullptr};
  std::function<std::unique_ptr<Acts::AnyMoveOnly>(const py::object& obj)>
      fromPython{nullptr};
  const std::type_info* typeinfo{nullptr};
  std::uint64_t typeHash{0};
};

std::unordered_map<PyObject*, RegistryEntry>& instance() {
  static std::unordered_map<PyObject*, RegistryEntry> map;
  return map;
}

}  // namespace

void WhiteBoardRegistry::registerTypeImpl(PyObject* pyType,
                                          ToPythonFunction toPython,
                                          FromPythonFunction fromPython,
                                          const std::type_info* typeinfo,
                                          std::uint64_t typeHash) {
  instance()[pyType] = {
      .toPython = std::move(toPython),
      .fromPython = std::move(fromPython),
      .typeinfo = typeinfo,
      .typeHash = typeHash,
  };
}

const WhiteBoardRegistry::EntryHandle* WhiteBoardRegistry::find(
    PyObject* pyType) noexcept {
  if (auto it = instance().find(pyType); it != instance().end()) {
    return &it->second;
  }
  return nullptr;
}

const std::type_info* WhiteBoardRegistry::typeInfo(
    const EntryHandle* entry) noexcept {
  auto* e = static_cast<const RegistryEntry*>(entry);
  return (e != nullptr) ? e->typeinfo : nullptr;
}

std::uint64_t WhiteBoardRegistry::typeHash(const EntryHandle* entry) noexcept {
  auto* e = static_cast<const RegistryEntry*>(entry);
  return (e != nullptr) ? e->typeHash : 0;
}

PyObject* WhiteBoardRegistry::toPython(const EntryHandle* entry,
                                       const Acts::AnyMoveOnly& any,
                                       PyObject* wbPy) {
  if (entry == nullptr) {
    throw py::type_error("Type is not registered for WhiteBoard access");
  }
  auto* e = static_cast<const RegistryEntry*>(entry);
  py::object pyWb = py::reinterpret_borrow<py::object>(wbPy);
  py::object obj = e->toPython(any, pyWb);
  return obj.release().ptr();
}

std::unique_ptr<Acts::AnyMoveOnly> WhiteBoardRegistry::fromPython(
    const EntryHandle* entry, PyObject* obj) {
  if (entry == nullptr) {
    throw py::type_error("Type is not registered for WhiteBoard access");
  }
  auto* e = static_cast<const RegistryEntry*>(entry);
  py::object pyObj = py::reinterpret_borrow<py::object>(obj);
  return e->fromPython(pyObj);
}

}  // namespace ActsPython
