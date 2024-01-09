// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

#include <optional>
#include <type_traits>

namespace Acts {

class Surface;

namespace Concepts {

// nested types that must be available
template <typename T>
using TypeScalar = typename T::Scalar;
template <typename T>
using TypeParametersVector = typename T::ParametersVector;
template <typename T>
using TypeCovarianceMatrix = typename T::CovarianceMatrix;

template <typename T>
using ReturnTypeParameters = decltype(std::declval<T>().parameters());
template <typename T>
using ReturnTypeCovariance = decltype(std::declval<T>().covariance());
template <typename T>
using ReturnTypeFourPositionFromContext =
    decltype(std::declval<T>().fourPosition(std::declval<GeometryContext>()));
template <typename T>
using ReturnTypeFourPosition = decltype(std::declval<T>().fourPosition());
template <typename T>
using ReturnTypePositionFromContext =
    decltype(std::declval<T>().position(std::declval<GeometryContext>()));
template <typename T>
using ReturnTypePosition = decltype(std::declval<T>().position());
template <typename T>
using ReturnTypeTime = decltype(std::declval<T>().time());
template <typename T>
using ReturnTypeDirection = decltype(std::declval<T>().direction());
template <typename T>
using ReturnTypeAbsoluteMomentum =
    decltype(std::declval<T>().absoluteMomentum());
template <typename T>
using ReturnTypeCharge = decltype(std::declval<T>().charge());
template <typename T>
using ReturnTypeReferenceSurface =
    decltype(std::declval<T>().referenceSurface());

template <typename T>
struct BoundTrackParametersConceptImpl {
  // check for required nested types
  constexpr static bool hasTypeScalar = exists<TypeScalar, const T>;
  constexpr static bool hasTypeParametersVector =
      exists<TypeParametersVector, const T>;
  constexpr static bool hasTypeCovarianceMatrix =
      exists<TypeCovarianceMatrix, const T>;

  // check for required methods
  constexpr static bool hasMethodParameters =
      std::is_convertible_v<ReturnTypeParameters<T>, BoundVector>;
  constexpr static bool hasMethodCovariance =
      std::is_convertible_v<ReturnTypeCovariance<T>,
                            std::optional<BoundSquareMatrix>>;
  constexpr static bool hasMethodFourPositionFromContext =
      identical_to<Vector4, ReturnTypeFourPositionFromContext, const T>;
  constexpr static bool hasMethodPositionFromContext =
      identical_to<Vector3, ReturnTypePositionFromContext, const T>;
  constexpr static bool hasMethodTime =
      identical_to<TypeScalar<T>, ReturnTypeTime, const T>;
  constexpr static bool hasMethodDirection =
      identical_to<Vector3, ReturnTypeDirection, const T>;
  constexpr static bool hasMethodAbsoluteMomentum =
      identical_to<TypeScalar<T>, ReturnTypeAbsoluteMomentum, const T>;
  constexpr static bool hasMethodCharge =
      identical_to<TypeScalar<T>, ReturnTypeCharge, const T>;
  constexpr static bool hasMethodReferenceSurface =
      identical_to<const Surface&, ReturnTypeReferenceSurface, const T>;

  // provide meaningful error messages in case of non-compliance
  static_assert(hasTypeScalar, "Scalar type is missing");
  static_assert(hasTypeParametersVector, "Parameters vector type is missing");
  static_assert(hasTypeCovarianceMatrix, "Covariance matrix type is missing");
  static_assert(hasMethodParameters, "Missing or invalid 'parameters' method");
  static_assert(hasMethodCovariance, "Missing or invalid 'covariance' method");
  static_assert(hasMethodFourPositionFromContext,
                "Missing or invalid 'fourPosition' method");
  static_assert(hasMethodPositionFromContext,
                "Missing or invalid 'position' method");
  static_assert(hasMethodTime, "Missing or invalid 'time' method");
  static_assert(hasMethodDirection, "Missing or invalid 'direction' method");
  static_assert(hasMethodAbsoluteMomentum,
                "Missing or invalid 'absoluteMomentum' method");
  static_assert(hasMethodCharge, "Missing or invalid 'charge' method");
  static_assert(hasMethodReferenceSurface,
                "Missing or invalid 'referenceSurface' method");

  constexpr static bool value =
      require<hasTypeScalar, hasTypeParametersVector, hasTypeCovarianceMatrix,
              hasMethodParameters, hasMethodCovariance,
              hasMethodFourPositionFromContext, hasMethodPositionFromContext,
              hasMethodTime, hasMethodDirection, hasMethodAbsoluteMomentum,
              hasMethodCharge, hasMethodReferenceSurface>;
};

template <typename T>
struct FreeTrackParametersConceptImpl {
  // check for required nested types
  constexpr static bool hasTypeScalar = exists<TypeScalar, const T>;
  constexpr static bool hasTypeParametersVector =
      exists<TypeParametersVector, const T>;
  constexpr static bool hasTypeCovarianceMatrix =
      exists<TypeCovarianceMatrix, const T>;

  // check for required methods
  constexpr static bool hasMethodParameters =
      std::is_convertible_v<ReturnTypeParameters<T>, FreeVector>;
  constexpr static bool hasMethodCovariance =
      std::is_convertible_v<ReturnTypeCovariance<T>,
                            std::optional<FreeSquareMatrix>>;
  constexpr static bool hasMethodFourPosition =
      identical_to<Vector4, ReturnTypeFourPosition, const T>;
  constexpr static bool hasMethodPosition =
      identical_to<Vector3, ReturnTypePosition, const T>;
  constexpr static bool hasMethodTime =
      identical_to<TypeScalar<T>, ReturnTypeTime, const T>;
  constexpr static bool hasMethodDirection =
      identical_to<Vector3, ReturnTypeDirection, const T>;
  constexpr static bool hasMethodAbsoluteMomentum =
      identical_to<TypeScalar<T>, ReturnTypeAbsoluteMomentum, const T>;
  constexpr static bool hasMethodCharge =
      identical_to<TypeScalar<T>, ReturnTypeCharge, const T>;

  // provide meaningful error messages in case of non-compliance
  static_assert(hasTypeScalar, "Scalar type is missing");
  static_assert(hasTypeParametersVector, "Parameters vector type is missing");
  static_assert(hasTypeCovarianceMatrix, "Covariance matrix type is missing");
  static_assert(hasMethodParameters, "Missing or invalid 'parameters' method");
  static_assert(hasMethodCovariance, "Missing or invalid 'covariance' method");
  static_assert(hasMethodFourPosition,
                "Missing or invalid 'fourPosition' method");
  static_assert(hasMethodPosition, "Missing or invalid 'position' method");
  static_assert(hasMethodTime, "Missing or invalid 'time' method");
  static_assert(hasMethodDirection, "Missing or invalid 'direction' method");
  static_assert(hasMethodAbsoluteMomentum,
                "Missing or invalid 'absoluteMomentum' method");
  static_assert(hasMethodCharge, "Missing or invalid 'charge' method");

  constexpr static bool value =
      require<hasTypeScalar, hasTypeParametersVector, hasTypeCovarianceMatrix,
              hasMethodParameters, hasMethodCovariance, hasMethodFourPosition,
              hasMethodPosition, hasMethodTime, hasMethodDirection,
              hasMethodAbsoluteMomentum, hasMethodCharge>;
};

template <typename parameters_t>
constexpr bool BoundTrackParametersConcept =
    BoundTrackParametersConceptImpl<parameters_t>::value;

template <typename parameters_t>
constexpr bool FreeTrackParametersConcept =
    FreeTrackParametersConceptImpl<parameters_t>::value;

}  // namespace Concepts
}  // namespace Acts
