// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <system_error>
#include <type_traits>

class GeoFullPhysVol;

namespace ActsPlugins {

/// @ingroup errors
enum class GeoModelConversionError {
  // ensure all values are non-zero
  /// Wrong shape provided for this converter
  WrongShapeForConverter = 1,
  /// Shape parameters can not be converted to Surface representation
  InvalidShapeParameters,
  /// Unknown Shape provided, no converter available
  UnkownShape,
  /// No logical volume found for the shape
  MissingLogicalVolume
};

std::error_code make_error_code(GeoModelConversionError e);

}  // namespace ActsPlugins

namespace std {
// register with STL
template <>
struct is_error_code_enum<ActsPlugins::GeoModelConversionError>
    : std::true_type {};
}  // namespace std
