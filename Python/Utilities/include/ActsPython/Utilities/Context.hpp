// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <string>
#include <tuple>
#include <unordered_map>

#include <pybind11/pybind11.h>

namespace ActsPython {

/// A Context class that holds the Python modules used in ACTS
struct Context {
  std::unordered_map<std::string, pybind11::module_,
                     std::hash<std::string_view>, std::equal_to<>>
      modules;

  /// @brief Retrieve a module by name
  /// @param name the name of the module you retrieve
  /// @return a reference to the module
  pybind11::module_& get(const std::string& name) { return modules.at(name); }

  /// @brief Retrieve multiple modules
  /// @tparam ...Args parameter pack of module names
  /// @param ...args
  /// @return
  template <typename... Args>
  auto get(Args&&... args)
    requires(sizeof...(Args) >= 2)
  {
    return std::make_tuple((modules.at(args))...);
  }
};

}  // namespace ActsPython
