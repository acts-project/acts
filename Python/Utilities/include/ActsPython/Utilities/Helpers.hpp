// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <algorithm>

#include <pybind11/pybind11.h>

namespace ActsPython {

/// This method calls the acts adapter to patch the classes with a config object
///
/// @param m the module to patch
inline void patchClassesWithConfig(pybind11::module_& m) {
  pybind11::module::import("acts._adapter").attr("_patch_config")(m);
}

/// @brief This method patches a class with kwargs arguments
/// @tparam T the config or class object type
/// @param c the config or class object to patch
template <typename T>
void patchKwargsConstructor(T& c) {
  pybind11::module::import("acts._adapter").attr("_patchKwargsConstructor")(c);
}

/// @brief  This sets a range property on a class
template <typename T, typename Ur, typename Ut>
void pythonRangeProperty(T& obj, const std::string& name, Ur Ut::*begin,
                         Ur Ut::*end) {
  obj.def_property(
      name.c_str(), [=](Ut& self) { return std::pair{self.*begin, self.*end}; },
      [=](Ut& self, std::pair<Ur, Ur> p) {
        self.*begin = p.first;
        self.*end = p.second;
      });
}

}  // namespace ActsPython
