// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"

#include <ostream>
#include <sstream>

namespace Acts {

/// Enum which determines how the axis handle its outer boundaries
/// possible values values
enum class AxisBoundaryType {
  /// Default behaviour: out of bounds
  /// positions are filled into the over or underflow bins
  Open,
  /// Out-of-bounds positions resolve to first/last bin
  /// respectively
  Bound,
  /// Out-of-bounds positions resolve to the outermost
  /// bin on the opposite side
  Closed,
};

/// Tag helper type for Axis constructors with class template deduction
/// @tparam bdt the boundary type
template <AxisBoundaryType bdt>
struct AxisBoundaryTypeTag {};

/// Convenience typedefs for AxisBoundaryTypeTag
constexpr auto AxisOpen = AxisBoundaryTypeTag<AxisBoundaryType::Open>{};
constexpr auto AxisBound = AxisBoundaryTypeTag<AxisBoundaryType::Bound>{};
constexpr auto AxisClosed = AxisBoundaryTypeTag<AxisBoundaryType::Closed>{};

inline std::ostream& operator<<(std::ostream& os, AxisBoundaryType bdt) {
  switch (bdt) {
    case AxisBoundaryType::Open:
      os << "Open";
      break;
    case AxisBoundaryType::Bound:
      os << "Bound";
      break;
    case AxisBoundaryType::Closed:
      os << "Closed";
      break;
  }
  return os;
}

inline std::string axisBoundaryTypeToString(AxisBoundaryType bdt) {
  std::stringstream ss;
  ss << bdt;
  return ss.str();
}

/// Enum which determines the binning type of the axis
enum class AxisType {
  /// An axis where all bins have the same size
  Equidistant,
  /// An axis where bins can have different sizes
  Variable,
};

inline std::ostream& operator<<(std::ostream& os, AxisType type) {
  switch (type) {
    case AxisType::Equidistant:
      os << "Equidistant";
      break;
    case AxisType::Variable:
      os << "Variable";
      break;
  }
  return os;
}

inline std::string axisTypeToString(AxisType type) {
  std::stringstream ss;
  ss << type;
  return ss.str();
}

/// @brief calculate bin indices from a given binning structure
///
/// This class provides some basic functionality for calculating bin indices
/// for a given binning configuration. Both equidistant as well as variable
/// binning structures are supported.
///
/// Bin intervals are defined such that the lower bound is closed and the
/// upper bound is open.
///
/// @tparam equidistant flag whether binning is equidistant (@c true)
///                     or not (@c false)
template <AxisType type, AxisBoundaryType bdt = AxisBoundaryType::Open>
class Axis;

Axis(double min, double max,
     std::size_t bins) -> Axis<AxisType::Equidistant, AxisBoundaryType::Open>;

template <AxisBoundaryType bdt>
Axis(AxisBoundaryTypeTag<bdt> /*bdt*/, double min, double max,
     std::size_t bins) -> Axis<AxisType::Equidistant, bdt>;

Axis(std::vector<double> bins)
    -> Axis<AxisType::Variable, AxisBoundaryType::Open>;

template <AxisBoundaryType bdt>
Axis(AxisBoundaryTypeTag<bdt> /*bdt*/,
     std::vector<double> bins) -> Axis<AxisType::Variable, bdt>;

/// @enum AxisDirection defines which parameter is described by this
/// axis or coordinate dimension. This can be either local or global
/// coordinates, the actual meaning is defined by the context in which it is
/// called.
enum class AxisDirection : int {
  /// The x direction
  AxisX = 0,
  /// The y direction
  AxisY = 1,
  /// The z direction
  AxisZ = 2,
  /// The radial direction
  AxisR = 3,
  /// The azimuthal direction
  AxisPhi = 4,
  /// The radial-azimuthal direction
  AxisRPhi = 5,
  /// The theta direction
  AxisTheta = 6,
  /// The eta direction
  AxisEta = 7,
  /// The magnitude direction
  AxisMag = 8
};

/// Returns the total number of binningvalues
/// @return the number of binning values
constexpr std::size_t numAxisDirections() {
  return 9;
}

/// @brief convert an AxisDirection from a string
/// @param adName the name of the AxisDirection
/// @note temporarily backwards compatibility enabled
inline AxisDirection axisDirectionFromName(const std::string& adName) {
  if (adName == "AxisX" || adName == "binX")
    return AxisDirection::AxisX;
  if (adName == "AxisY" || adName == "binY")
    return AxisDirection::AxisY;
  if (adName == "AxisZ" || adName == "binZ")
    return AxisDirection::AxisZ;
  if (adName == "AxisR" || adName == "binR")
    return AxisDirection::AxisR;
  if (adName == "AxisPhi" || adName == "binPhi")
    return AxisDirection::AxisPhi;
  if (adName == "AxisRPhi" || adName == "binRPhi")
    return AxisDirection::AxisRPhi;
  if (adName == "AxisTheta" || adName == "binH")
    return AxisDirection::AxisTheta;
  if (adName == "AxisEta" || adName == "binEta")
    return AxisDirection::AxisEta;
  if (adName == "AxisMag" || adName == "binMag")
    return AxisDirection::AxisMag;
  throw std::invalid_argument("Unknown AxisDirection name: " + adName);
}

inline std::ostream& operator<<(std::ostream& os, AxisDirection descr) {
  switch (descr) {
    case AxisDirection::AxisX:
      os << "AxisX";
      break;
    case AxisDirection::AxisY:
      os << "AxisY";
      break;
    case AxisDirection::AxisZ:
      os << "AxisZ";
      break;
    case AxisDirection::AxisR:
      os << "AxisR";
      break;
    case AxisDirection::AxisPhi:
      os << "AxisPhi";
      break;
    case AxisDirection::AxisRPhi:
      os << "AxisRPhi";
      break;
    case AxisDirection::AxisTheta:
      os << "AxisTheta";
      break;
    case AxisDirection::AxisEta:
      os << "AxisEta";
      break;
    case AxisDirection::AxisMag:
      os << "AxisMag";
      break;
  }
  return os;
}

inline std::string axisDirectionToString(AxisDirection descr) {
  std::stringstream ss;
  ss << descr;
  return ss.str();
}

inline std::vector<AxisDirection> allAxisDirections() {
  return {
      AxisDirection::AxisX,     AxisDirection::AxisY,   AxisDirection::AxisZ,
      AxisDirection::AxisR,     AxisDirection::AxisPhi, AxisDirection::AxisRPhi,
      AxisDirection::AxisTheta, AxisDirection::AxisEta, AxisDirection::AxisMag};
}

}  // namespace Acts
