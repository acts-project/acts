// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/HashedString.hpp"

namespace Acts::detail {
namespace {
struct TypeHashProbeA {};
struct TypeHashProbeB {};
}  // namespace

// Fail the build on a compiler that omits the template arguments from
// std::source_location::function_name(), rather than silently collide hashes.
static_assert(typeHash<TypeHashProbeA>() != typeHash<TypeHashProbeB>(),
              "typeHash<T> does not distinguish types on this compiler: "
              "std::source_location::function_name() apparently omits template "
              "arguments. Fall back to __PRETTY_FUNCTION__ / __FUNCSIG__.");

}  // namespace Acts::detail
