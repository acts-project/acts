// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.


namespace Acts {

/// TODO This probably belongs to Acts/Core/Utilities
/// Helper object to allow the overload pattern like shown e.g. in
/// https://www.modernescpp.com/index.php/visiting-a-std-variant-with-the-overload-pattern
/// This is especially useful for using a std::variant with std::visit.
template <typename... Ts>
struct Overload : Ts... {
  using Ts::operator()...;
};

/// Deduction guide for the Overload class
template <class... Ts>
Overload(Ts...) -> Overload<Ts...>;

}  // namespace Acts
