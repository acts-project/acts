// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Utilities/Logger.hpp"

#include <concepts>

#include <boost/preprocessor/seq/for_each.hpp>
#include <boost/preprocessor/variadic/to_seq.hpp>
#include <pybind11/pybind11.h>

namespace ActsExamples {
struct AlgorithmContext;
enum class ProcessCode;
}  // namespace ActsExamples

namespace ActsPython::Concepts {
template <typename T>
concept has_write_method =
    requires(T& t, const ActsExamples::AlgorithmContext& ctx) {
      { t.write(ctx) } -> std::same_as<ActsExamples::ProcessCode>;
    };
}  // namespace ActsPython::Concepts

#define ACTS_PYTHON_MEMBER(name) \
  _binding_instance.def_readwrite(#name, &_struct_type::name)

#define ACTS_PYTHON_STRUCT_BEGIN(object)               \
  {                                                    \
    [[maybe_unused]] auto& _binding_instance = object; \
    using _struct_type = decltype(object)::type;       \
    do {                                               \
    } while (0)

#define ACTS_PYTHON_STRUCT_END() \
  }                              \
  do {                           \
  } while (0)

/// This macro is needed to use the BOOST_PP_SEQ_FOR_EACH loop macro
#define ACTS_PYTHON_MEMBER_LOOP(r, data, elem) ACTS_PYTHON_MEMBER(elem);

/// Macro that accepts a variadic set of member names that are to be registered
/// into an object as read-write fields
#define ACTS_PYTHON_STRUCT(object, ...)                          \
  do {                                                           \
    ACTS_PYTHON_STRUCT_BEGIN(object);                            \
    BOOST_PP_SEQ_FOR_EACH(ACTS_PYTHON_MEMBER_LOOP, _,            \
                          BOOST_PP_VARIADIC_TO_SEQ(__VA_ARGS__)) \
    ACTS_PYTHON_STRUCT_END();                                    \
  } while (0)

template <typename A, typename B>
auto declareAlgorithm(pybind11::module_& m, const char* name) {
  using Config = typename A::Config;
  namespace py = pybind11;
  auto alg =
      py::class_<A, B, std::shared_ptr<A>>(m, name)
          .def(py::init([name](const Config& cfg, Acts::Logging::Level level) {
                 return std::make_shared<A>(
                     cfg, Acts::getDefaultLogger(name, level));
               }),
               py::arg("config"), py::arg("level"))
          .def(py::init([](const Config& cfg,
                           std::unique_ptr<const Acts::Logger> logger) {
                 return std::make_shared<A>(cfg, std::move(logger));
               }),
               py::arg("config"), py::arg("logger"))
          .def_property_readonly("config", &A::config);
  auto c = py::class_<Config>(alg, "Config");
  if constexpr (std::is_default_constructible_v<Config>) {
    c.def(py::init<>());
  }
  return std::tuple{alg, c};
}

/// A macro that uses Boost.Preprocessor to create the python binding for and
/// algorithm and the additional config struct.
#define ACTS_PYTHON_DECLARE_ALGORITHM(algorithm, mod, name, ...)          \
  do {                                                                    \
    auto [alg, c] =                                                       \
        declareAlgorithm<algorithm, ActsExamples::IAlgorithm>(mod, name); \
    ACTS_PYTHON_STRUCT(c, __VA_ARGS__);                                   \
  } while (0)

/// Similar as above for writers
#define ACTS_PYTHON_DECLARE_WRITER(writer, mod, name, ...)                  \
  do {                                                                      \
    using Writer = writer;                                                  \
    using Config = Writer::Config;                                          \
    auto w =                                                                \
        py::class_<Writer, ActsExamples::IWriter, std::shared_ptr<Writer>>( \
            mod, name)                                                      \
            .def(py::init<const Config&, Acts::Logging::Level>(),           \
                 py::arg("config"), py::arg("level"))                       \
            .def_property_readonly("config", &Writer::config);              \
                                                                            \
    constexpr bool has_write_method =                                       \
        ActsPython::Concepts::has_write_method<Writer>;                     \
                                                                            \
    if constexpr (has_write_method) {                                       \
      w.def("write", &Writer::write);                                       \
    }                                                                       \
    auto c = py::class_<Config>(w, "Config").def(py::init<>());             \
    ACTS_PYTHON_STRUCT(c, __VA_ARGS__);                                     \
  } while (0)

/// Similar as above for readers
#define ACTS_PYTHON_DECLARE_READER(reader, mod, name, ...)                  \
  do {                                                                      \
    using Reader = reader;                                                  \
    using Config = Reader::Config;                                          \
    auto r =                                                                \
        py::class_<Reader, ActsExamples::IReader, std::shared_ptr<Reader>>( \
            mod, name)                                                      \
            .def(py::init<const Config&, Acts::Logging::Level>(),           \
                 py::arg("config"), py::arg("level"))                       \
            .def_property_readonly("config", &Reader::config);              \
                                                                            \
    auto c = py::class_<Config>(r, "Config").def(py::init<>());             \
    ACTS_PYTHON_STRUCT(c, __VA_ARGS__);                                     \
  } while (0)
