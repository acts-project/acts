// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"

#include <algorithm>
#include <optional>
#include <tuple>
#include <variant>

namespace Acts {

using VariantCovariance = std::variant<BoundSquareMatrix, FreeSquareMatrix>;

using VariantTransportJacobian =
    std::variant<BoundMatrix, BoundToFreeMatrix, FreeToBoundMatrix, FreeMatrix>;

/// Helper struct holding the necessary cache for the Covariance
/// transport in various parametrisations.
struct CovarianceCache {
  /// Internal cache state, indicates correct setup
  bool applyTransport = false;

  /// Internal cache state, at surface indication
  const Surface* atSurface = nullptr;

  /// Internal cache state, at position
  Vector3 atPosition = Vector3::Zero();

  /// Variant: the currently held covariance
  VariantCovariance covariance;

  /// Optional for starting from bound or curvilinear
  std::optional<BoundToFreeMatrix> boundToFreeJacobian = std::nullopt;
  /// Options for starting from free
  std::optional<ActsMatrix<8, 7>> anglesToDirectionJacobian = std::nullopt;
  std::optional<ActsMatrix<7, 8>> directionToAnglesJacobian = std::nullopt;

  /// Non-variant: the free transport jacobian
  FreeMatrix freeTransportJacobian = FreeMatrix::Identity();
  /// Non-variant: the free derivatives
  FreeVector freeToPathDerivatives = FreeVector::Zero();

  /// Defaulted constructor, gives invalid cache
  CovarianceCache() = default;

  /// Constructor from bound & surface
  ///
  /// @param gctx The current geometry context
  /// @param surface The surface of the bound representation
  /// @param position The position of the representation
  /// @param boundParameters The bound parameters at the surface
  /// @param boundCovariance The bound covariance to be propagated
  ///
  /// This constructor will set the variant covariance type to
  /// a bound matrix, remember the surface & establish the
  /// jacobian between bound and free parametrisation.
  CovarianceCache(const GeometryContext& gctx, const Surface& surface,
                  const Vector3& position, const BoundVector& boundParameters,
                  const BoundSquareMatrix& boundCovariance);

  /// Constructor from curvilinear
  ///
  /// @param position The position of the representation
  /// @param direction The direction of at the representation
  /// @param boundCovariance The bound covariance to be propagated
  ///
  /// This constructor will set the variant covariance type to
  /// a bound matrix, remember the surface & establish the
  /// jacobian between bound and free parametrisation.
  CovarianceCache(const Vector3& position, const Vector3& direction,
                  const BoundSquareMatrix& boundCovariance);

  /// Construction from free
  ///
  /// @param freeParameters The free parameters
  /// @param freeCovariance The free covariance to be propagated
  ///
  CovarianceCache(const FreeVector& freeParameters,
                  const FreeSquareMatrix& freeCovariance);
};

/// Transport the covariance to a new (surface) bound definition
///
/// @param gctx [in] The current geometry context
/// @param surface [in] The surface of the bound representation
/// @param parameters [in] The free parametrisation on the surface
/// @param cCache [in,out] the covariance cache (to be modified)
///
/// @return a variant transport jacobian
std::tuple<VariantCovariance, VariantTransportJacobian>
transportCovarianceToBound(const GeometryContext& gctx, const Surface& surface,
                           const FreeVector& parameters,
                           CovarianceCache& cCache);

/// Transport the covariance to a new (curvilinear) bound definition
///
/// @param direction [in] The track direction at the new curvilinear definition
/// @param cCache [in,out]  the covariance cache (to be modified)
std::tuple<VariantCovariance, VariantTransportJacobian>
transportCovarianceToCurvilinear(const Vector3& direction,
                                 CovarianceCache& cCache);

/// Transport the covariance to a new free definition
///
/// @param cCache [in,out] the covariance cache (to be modified)
std::tuple<VariantCovariance, VariantTransportJacobian>
transportCovarianceToFree(CovarianceCache& cCache);

}  // namespace Acts
