// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Utilities/MathHelpers.hpp"

namespace Acts::detail {
/// @brief Auxiliary class to coherently sort Transforms and Vectors such that
///        a strict ordering between two objects may be defined and hence
///        a sorted map or set of these objects becomes possible. In the case of
///        Vectors, the sorter compares the differences of each component. As
///        soon as a difference exceeds the predefined translational tolerance
///        threshold the < operator is assigned based on this difference. In the
///        case of rotations, the three Euler angles of the rotation matrix are
///        evaluated and compared in an analogous way
class TransformComparator {
 public:
  /// @brief Default constructor setting the predefined tolerance values
  ///        for translation & rotation
  TransformComparator() = default;
  /// @brief Constructor with an adaption of the translation & rotTolerance
  /// @param transTolerance: Tolerance value within the difference of the
  ///         i-th component between two vectors is considered to be 0
  /// @param rotTolerance: Tolerance value within the difference of the
  ///        i-th Euler angle between two Rotations is considered to be 0
  TransformComparator(const double transTolerance, const double rotTolerance);
  /// @brief Generic comparison function between two kSize-dimensional vectors
  ///        Returns 0 if all components agree within the translational
  ///        tolerance. As soon as one component differs, 1 is returned if the
  ///        component from vector a is larger and otherwise -1
  /// @param a: Reference to the first vector to compare
  /// @param b: Reference to the second vector to compare
  template <unsigned int kSize>
  int compare(const Acts::Vector<kSize>& a,
              const Acts::Vector<kSize>& b) const {
    for (unsigned int i = 0; i < kSize; ++i) {
      const double diff = a[i] - b[i];
      if (std::abs(diff) > m_tolTrans) {
        return copySign(1, diff);
      }
    }
    return 0;
  }
  /// @brief Brief compares the three Euler angles of the two rotation matrices
  ///        If all angles agree within the rotational tolerance, 0 is returned.
  ///        Otherwise, -1 or 1 is returned depending on whether the first
  ///        differing i-th angle from a or b is larger.
  /// @param a: Reference to the first rotation matrix to compare
  /// @param b: Reference to the second rotation matrix to compare
  int compare(const Acts::RotationMatrix3& a,
              const Acts::RotationMatrix3& b) const;
  /// @brief Compares two transforms. First, it's checked whether the translational
  ///        components between the two differ and if not the comparison between
  ///        the two rotation matrices is returned
  /// @param a: Reference to the first transform to compare
  /// @param b: Reference to the second transform to compare
  int compare(const Acts::Transform3& a, const Acts::Transform3& b) const;

  /// @brief Implementation of the < operator for Transforms
  bool operator()(const Acts::Transform3& a, const Acts::Transform3& b) const;
  /// @brief Implementation of the < operator for RotationMatrices
  bool operator()(const Acts::RotationMatrix3& a,
                  const Acts::RotationMatrix3& b) const;
  /// @brief Implementation of the < operator for 3-vectors
  bool operator()(const Acts::Vector3& a, const Acts::Vector3& b) const;
  /// @brief Implementation of the < operator for 2-vectors
  bool operator()(const Acts::Vector2& a, const Acts::Vector2& b) const;

 private:
  /** @brief Maximum tolerance per translational vector component */
  double m_tolTrans{0.1 * Acts::UnitConstants::um};
  /** @brief Maximum tolerance per euler angle */
  double m_tolRot{0.01 * Acts::UnitConstants::mrad};
};
}  // namespace Acts::detail
