// This file is part of the Acts project.
//
// Copyright (C) 2019-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/TypeTraits.hpp"

namespace Acts {

class Surface;

namespace Concepts {
using namespace Acts::concept;

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
using ReturnTypePosition4FromContext =
    decltype(std::declval<T>().position4(GeometryContext()));
template <typename T>
using ReturnTypePosition4 = decltype(std::declval<T>().position4());
template <typename T>
using ReturnTypePositionFromContext =
    decltype(std::declval<T>().position(GeometryContext()));
template <typename T>
using ReturnTypePosition = decltype(std::declval<T>().position());
template <typename T>
using ReturnTypeTime = decltype(std::declval<T>().time());
template <typename T>
using ReturnTypeDirectionUnit = decltype(std::declval<T>().directionUnit());
template <typename T>
using ReturnTypeAbsMomentum = decltype(std::declval<T>().absMomentum());
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
      identical_to<BoundVector, ReturnTypeParameters, const T>;
  constexpr static bool hasMethodCovariance =
      identical_to<const std::optional<BoundSymMatrix>&, ReturnTypeCovariance,
                   const T>;
  constexpr static bool hasMethodPosition4FromContext =
      identical_to<Vector4D, ReturnTypePosition4FromContext, const T>;
  constexpr static bool hasMethodPositionFromContext =
      identical_to<Vector3D, ReturnTypePositionFromContext, const T>;
  constexpr static bool hasMethodTime =
      identical_to<double, ReturnTypeTime, const T>;
  constexpr static bool hasMethodDirectionUnit =
      identical_to<Vector3D, ReturnTypeDirectionUnit, const T>;
  constexpr static bool hasMethodAbsMomentum =
      identical_to<double, ReturnTypeAbsMomentum, const T>;
  constexpr static bool hasMethodCharge =
      identical_to<double, ReturnTypeCharge, const T>;
  constexpr static bool hasMethodReferenceSurface =
      identical_to<const Surface&, ReturnTypeReferenceSurface, const T>;

  // provide meaningful error messages in case of non-compliance
  static_assert(hasTypeScalar, "Scalar type is missing");
  static_assert(hasTypeParametersVector, "Parameters vector type is missing");
  static_assert(hasTypeCovarianceMatrix, "Covariance matrix type is missing");
  static_assert(hasMethodParameters, "Missing or invalid 'parameters' method");
  static_assert(hasMethodCovariance, "Missing or invalid 'covariance' method");
  static_assert(hasMethodPositionFromContext,
                "Missing or invalid 'position' method");
  static_assert(hasMethodTime, "Missing or invalid 'time' method");
  static_assert(hasMethodDirectionUnit,
                "Missing or invalid 'directionUnit' method");
  static_assert(hasMethodAbsMomentum,
                "Missing or invalid 'absMomentum' method");
  static_assert(hasMethodCharge, "Missing or invalid 'charge' method");
  static_assert(hasMethodReferenceSurface,
                "Missing or invalid 'referenceSurface' method");

  constexpr static bool value =
      require<hasTypeScalar, hasTypeParametersVector, hasTypeCovarianceMatrix,
              hasMethodParameters, hasMethodCovariance,
              hasMethodPosition4FromContext, hasMethodPositionFromContext,
              hasMethodTime, hasMethodDirectionUnit, hasMethodAbsMomentum,
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
      identical_to<FreeVector, ReturnTypeParameters, const T>;
  constexpr static bool hasMethodCovariance =
      identical_to<const std::optional<FreeSymMatrix>&, ReturnTypeCovariance,
                   const T>;
  constexpr static bool hasMethodPosition4 =
      identical_to<Vector4D, ReturnTypePosition4, const T>;
  constexpr static bool hasMethodPosition =
      identical_to<Vector3D, ReturnTypePosition, const T>;
  constexpr static bool hasMethodTime =
      identical_to<double, ReturnTypeTime, const T>;
  constexpr static bool hasMethodDirectionUnit =
      identical_to<Vector3D, ReturnTypeDirectionUnit, const T>;
  constexpr static bool hasMethodAbsMomentum =
      identical_to<double, ReturnTypeAbsMomentum, const T>;
  constexpr static bool hasMethodCharge =
      identical_to<double, ReturnTypeCharge, const T>;

  // provide meaningful error messages in case of non-compliance
  static_assert(hasTypeScalar, "Scalar type is missing");
  static_assert(hasTypeParametersVector, "Parameters vector type is missing");
  static_assert(hasTypeCovarianceMatrix, "Covariance matrix type is missing");
  static_assert(hasMethodParameters, "Missing or invalid 'parameters' method");
  static_assert(hasMethodCovariance, "Missing or invalid 'covariance' method");
  static_assert(hasMethodPosition4, "Missing or invalid 'position4' method");
  static_assert(hasMethodPosition, "Missing or invalid 'position' method");
  static_assert(hasMethodTime, "Missing or invalid 'time' method");
  static_assert(hasMethodDirectionUnit,
                "Missing or invalid 'directionUnit' method");
  static_assert(hasMethodAbsMomentum,
                "Missing or invalid 'absMomentum' method");
  static_assert(hasMethodCharge, "Missing or invalid 'charge' method");

  constexpr static bool value =
      require<hasTypeScalar, hasTypeParametersVector, hasTypeCovarianceMatrix,
              hasMethodParameters, hasMethodCovariance, hasMethodPosition4,
              hasMethodPosition, hasMethodTime, hasMethodDirectionUnit,
              hasMethodAbsMomentum, hasMethodCharge>;
};

template <typename parameters_t>
constexpr bool BoundTrackParametersConcept =
    BoundTrackParametersConceptImpl<parameters_t>::value;

template <typename parameters_t>
constexpr bool FreeTrackParametersConcept =
    FreeTrackParametersConceptImpl<parameters_t>::value;

}  // namespace Concepts
}  // namespace Acts
