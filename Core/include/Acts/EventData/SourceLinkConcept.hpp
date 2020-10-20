// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <type_traits>

namespace Acts {

class GeometryIdentifier;

namespace Concepts {

template <typename T>
struct SourceLinkConceptImpl {
  // these manual checks are roughly equivalent to the C++20 `regular` concept:
  // a type the behaves similar to built-in types, i.e. is default
  // initializable, copyable, and equality comparable.
  constexpr static bool isDefaultInitializable =
      std::is_default_constructible_v<T>;
  // this is not testing all neccesary functionality but should be sufficient
  constexpr static bool isCopyable = std::is_copy_constructible_v<T>;
  constexpr static bool isEqualityComparable =
      std::is_same_v<decltype(std::declval<T>() == std::declval<T>()), bool> and
      std::is_same_v<decltype(std::declval<T>() != std::declval<T>()), bool>;

  constexpr static bool hasGeometryIdAccessor =
      std::is_same_v<std::decay_t<decltype(std::declval<T>().geometryId())>,
                     GeometryIdentifier>;

  // provide meaningful error messages in case of non-compliance
  static_assert(isDefaultInitializable,
                "Source link must be default initializable");
  static_assert(isCopyable, "Source link must be copyable");
  static_assert(isEqualityComparable,
                "Source link must be equality comparable");
  static_assert(hasGeometryIdAccessor,
                "Source link must provide a `geometryId()` accessor");

  constexpr static bool value = isDefaultInitializable and isCopyable and
                                isEqualityComparable and hasGeometryIdAccessor;
};

}  // namespace Concepts

template <typename T>
constexpr bool SourceLinkConcept = Concepts::SourceLinkConceptImpl<T>::value;

}  // namespace Acts
