// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <iostream>
#include <type_traits>
#include <typeinfo>

#include <boost/core/demangle.hpp>

namespace Acts {

template <typename T, typename U>
T checked_cast(U* ptr)
  requires std::is_pointer_v<T>
{
  if (ptr == nullptr) {
    return nullptr;
  }

  const auto& typeA = typeid(std::remove_cv_t<std::remove_pointer_t<T>>);
  const auto& typeB = typeid(*ptr);

  if (typeA.hash_code() != typeB.hash_code()) {
    return nullptr;
  }

  return static_cast<T>(ptr);
}

template <typename T, typename U>
T checked_cast(U& ref)
  requires std::is_reference_v<T>
{
  const auto& typeA = typeid(std::remove_cv_t<std::remove_reference_t<T>>);
  const auto& typeB = typeid(ref);

  if (typeA.hash_code() != typeB.hash_code()) {
    throw std::bad_cast();
  }

  return static_cast<T>(ref);
}

}  // namespace Acts
