// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <cstddef>
#include <functional>

#include <pybind11/pybind11.h>

namespace ActsPython {

/// Several EventData proxies/iterators hold a raw `Container*` + index and are
/// handed to Python. When the backing container is transferred to the
/// whiteboard, pybind11's smart_holder disowns the Python wrapper and the C++
/// object is moved out and freed, leaving the proxy dangling (a later access is
/// a use-after-free -> SIGSEGV).
///
/// `keep_alive` cannot help: it tracks the Python object's lifetime, but after
/// the transfer the wrapper is a live-but-disowned shell. The data genuinely
/// leaves, so the correct behaviour is to fail loud. pybind11 already maintains
/// the authoritative "is it still owned" signal (set when the transfer does
/// `py::cast<std::unique_ptr<T>>`), and we consult it by tethering each
/// proxy/iterator to its parent container's `py::object` and re-validating on
/// access. Holding the `py::object` also keeps the container alive while the
/// proxy lives, which subsumes `keep_alive` for the garbage-collection path.

/// Return whether `owner`'s wrapped C++ container is still owned (not disowned
/// by a whiteboard transfer). Implemented by attempting to re-acquire the C++
/// reference: this throws iff the instance was disowned. Any pending Python
/// error from the failed cast is cleared so the caller can raise a clean one.
template <typename Container>
bool ownerAlive(const pybind11::object& owner) {
  try {
    owner.cast<const Container&>();
    return true;
  } catch (...) {
    PyErr_Clear();
    return false;
  }
}

/// Tethers a C++ proxy to its parent container's Python object. The disown
/// check is type-erased via a function pointer (`&ownerAlive<Container>`)
/// rather than a Container template parameter, so a single bound proxy type can
/// be produced by several different container types (e.g. a measurement proxy
/// returned by both MeasurementContainer and MeasurementSubset).
template <typename Proxy, typename Owner = void>
struct ProxyTether {
  using AliveFn = bool (*)(const pybind11::object&);

  pybind11::object owner;
  Proxy proxy;
  AliveFn aliveFn;

  /// Single-owner constructor: alive check is derived from Owner at compile
  /// time. Use when there is exactly one container type that can own the proxy.
  ProxyTether(pybind11::object owner_, Proxy proxy_)
    requires(!std::is_void_v<Owner>)
      : owner{std::move(owner_)},
        proxy{std::move(proxy_)},
        aliveFn{&ownerAlive<Owner>} {}

  /// Multi-owner constructor: alive check is supplied explicitly. Use when the
  /// same bound proxy type can be produced by more than one container type
  /// (e.g. MeasurementContainer and MeasurementSubset both yield the same
  /// ConstVariableBoundMeasurementProxy).
  ProxyTether(pybind11::object owner_, Proxy proxy_, AliveFn aliveFn_)
      : owner{std::move(owner_)}, proxy{std::move(proxy_)}, aliveFn{aliveFn_} {}

  /// Validate the container is still owned, then return the proxy. Throws
  /// pybind11::value_error (a Python ValueError) if it was consumed.
  const Proxy& checked() const {
    validate();
    return proxy;
  }
  Proxy& checked() {
    validate();
    return proxy;
  }

 private:
  void validate() const {
    if (!aliveFn(owner)) {
      throw pybind11::value_error(
          "proxy is no longer valid: its container was consumed (moved to the "
          "whiteboard)");
    }
  }
};

/// Build a read accessor for a tethered proxy: validate the tether, then call
/// `access` with the inner proxy. `Tether` is given explicitly at the call
/// site.
template <typename Tether, typename Access>
auto tetheredRead(Access access) {
  return [access = std::move(access)](const Tether& self) {
    return access(self.checked());
  };
}

/// Build a write accessor for a tethered proxy. `V` is the concrete value type,
/// given explicitly at the call site so pybind11 sees a concrete setter
/// signature (never `auto`).
template <typename Tether, typename V, typename Assign>
auto tetheredWrite(Assign assign) {
  return [assign = std::move(assign)](Tether& self, const V& value) {
    assign(self.checked(), value);
  };
}

/// A Python iterator that gates an underlying pybind11 iterator on the owner
/// still being alive. Used for containers whose iteration yields values/pairs
/// whose Python conversion is still registered (e.g. the flat multimaps).
template <typename Container>
struct IteratorTether {
  pybind11::object owner;
  pybind11::iterator inner;

  pybind11::object next() {
    if (!ownerAlive<Container>(owner)) {
      throw pybind11::value_error(
          "iterator is no longer valid: its container was consumed (moved to "
          "the whiteboard)");
    }
    return inner.attr("__next__")();
  }
};

/// Register a IteratorTether type. MUST be called before any container's
/// `__iter__` that returns it, or pybind11 raises an unregistered-type error.
template <typename Container>
void bindIteratorTether(pybind11::module_& m, const char* name) {
  pybind11::class_<IteratorTether<Container>>(m, name)
      .def("__iter__", [](pybind11::object self) { return self; })
      .def("__next__", &IteratorTether<Container>::next);
}

/// A Python iterator over an index-accessible container that yields freshly
/// tethered proxies. Used for containers whose proxy type is itself bound as a
/// ProxyTether (so a plain make_iterator over raw proxies would hit an
/// unregistered type). `maker(owner, container, index)` builds the tethered
/// proxy for each element.
template <typename Container>
struct IndexIteratorTether {
  using Maker = std::function<pybind11::object(const pybind11::object&,
                                               Container&, std::size_t)>;
  pybind11::object owner;
  std::size_t index{0};
  Maker maker;

  pybind11::object next() {
    Container* container = nullptr;
    try {
      container = &owner.cast<Container&>();
    } catch (...) {
      PyErr_Clear();
      throw pybind11::value_error(
          "iterator is no longer valid: its container was consumed (moved to "
          "the whiteboard)");
    }
    if (index >= container->size()) {
      throw pybind11::stop_iteration();
    }
    return maker(owner, *container, index++);
  }
};

/// Register a IndexIteratorTether type. MUST be called before any container's
/// `__iter__` that returns it.
template <typename Container>
void bindIndexIteratorTether(pybind11::module_& m, const char* name) {
  pybind11::class_<IndexIteratorTether<Container>>(m, name)
      .def("__iter__", [](pybind11::object self) { return self; })
      .def("__next__", &IndexIteratorTether<Container>::next);
}

}  // namespace ActsPython
