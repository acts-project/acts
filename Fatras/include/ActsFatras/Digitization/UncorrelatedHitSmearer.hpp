// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/ParameterSet.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Surfaces/SurfaceError.hpp"
#include "Acts/Utilities/Helpers.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"
#include "Acts/Utilities/Result.hpp"
#include "ActsFatras/Digitization/detail/ParametersCast.hpp"
#include "ActsFatras/Digitization/detail/ParametersSmearer.hpp"
#include "ActsFatras/EventData/Hit.hpp"

namespace ActsFatras {

/// Smearing functions definition:
/// - it takes the unsmeared parameter
/// - it returns the smeared parameter and a covariance
using SmearFunction =
    std::function<Acts::Result<std::pair<double, double>>(double)>;

/// Parameter smearer for fast digitisation that produces
///
///
/// The logic of this smearer is the following:
/// - Input is a single simulated ActsFatras::Hit which is wrapped into
/// - It returns smeared parameter vector and a covariance matrix
///
/// @note This smearer only supports uncorrelated smearing of parameters
///
class UncorrelatedHitSmearer {
 public:
  /// Generic implementation of a smearing meathod for bound parameters
  ///
  /// @tparam kParameters parameter pack describing the parameters to smear
  ///
  /// @param gctx The Geometry context (alignment, etc.)
  /// @param hit The Fatras hit to be smeared (will be source linked)
  /// @param sf The Surface to which this hit is attached
  /// @param sFunctions The smearing functions that are applied
  ///
  /// @return Smeared bound parameter vector and covariance matrix
  template <Acts::BoundIndices... kParameters>
  Acts::Result<Acts::ParameterSet<Acts::BoundIndices, kParameters...>>
  smearedParameterSet(const Acts::GeometryContext& gctx, const Hit& hit,
                      const Acts::Surface& sf,
                      const std::array<SmearFunction, sizeof...(kParameters)>&
                          sFunctions) const {
    using ParSet = Acts::ParameterSet<Acts::BoundIndices, kParameters...>;
    using Result = Acts::Result<ParSet>;
    using ParametersSmearer =
        detail::ParametersSmearer<Acts::BoundIndices, kParameters...>;

    auto dir = hit.unitDirection();
    auto gltResult = sf.globalToLocal(gctx, hit.position(), dir);
    if (not gltResult.ok()) {
      return Result(Acts::SurfaceError::GlobalPositionNotOnSurface);
    }
    const auto& lPosition = gltResult.value();

    typename ParSet::ParametersVector sParameters;
    typename ParSet::CovarianceMatrix sCovariance;
    sCovariance.setZero();
    fillBoundParameters(sParameters, lPosition, dir, hit.time(),
                        kParameters...);

    auto smearResult =
        ParametersSmearer::run(sParameters, sCovariance, sFunctions);
    if (not smearResult.ok()) {
      return Result(smearResult.error());
    }
    return ParSet{sCovariance, sParameters};
  }

  /// Generic implementation of a smearing meathod for free parameters
  ///
  /// @tparam kParameters parameter pack describing the parameters to smear
  ///
  /// @param gctx The Geometry context (alignment, etc.)
  /// @param hit The Fatras hit to be smeared (will be source linked)
  /// @param sFunctions The smearing functions that are applied
  ///
  /// @note uncorrelated smearing of the direction using the components
  ///       is not recommended
  ///
  /// @return Smeared free parameter vector and covariance matrix
  template <Acts::FreeIndices... kParameters>
  Acts::Result<Acts::ParameterSet<Acts::FreeIndices, kParameters...>>
  smearedParameterSet(const Hit& hit,
                      const std::array<SmearFunction, sizeof...(kParameters)>&
                          sFunctions) const {
    using ParSet = Acts::ParameterSet<Acts::FreeIndices, kParameters...>;
    using Result = Acts::Result<ParSet>;
    using ParametersSmearer =
        detail::ParametersSmearer<Acts::FreeIndices, kParameters...>;

    typename ParSet::FullParametersVector fParameters;
    fParameters.template segment<4>(0) = hit.position4();
    fParameters.template segment<3>(4) = hit.unitDirection();

    typename ParSet::ParametersVector sParameters =
        ParSet::projector() * fParameters;
    typename ParSet::CovarianceMatrix sCovariance;
    sCovariance.setZero();

    auto smearResult =
        ParametersSmearer::run(sParameters, sCovariance, sFunctions);
    if (not smearResult.ok()) {
      return Result(smearResult.error());
    }
    return ParSet{sCovariance, sParameters};
  }

 private:
  /// Apply the recursive smearing on the parameter pack,
  /// the base local position or time is used as unsmeared
  ///
  /// @tparam vec_t The parameter vector, static size known at compile time
  /// @tparam param_t The first of the pack
  /// @tparam paramts_t The remaining pack
  ///
  /// @param vector[in,out] The to be filled parameters vector
  /// @param lhit The local hit parameters
  /// @param dir The direction of the hit
  /// @param ctime The current time
  /// @param first The first of the pack
  /// @param rest The remaining pack
  ///
  template <typename vector_t, typename param_t, typename... kParameters_t>
  void fillBoundParameters(vector_t& vector, const Acts::Vector2D& lhit,
                           const Acts::Vector3D& dir, double ctime,
                           const param_t& first,
                           kParameters_t const&... rest) const {
    const size_t par = vector.rows() - sizeof...(rest) - 1;
    vector[par] = detail::BoundCasts[first](lhit, dir, ctime);
    if constexpr (sizeof...(rest) > 0) {
      fillBoundParameters(vector, lhit, dir, ctime, rest...);
    }
  }
};

}  // namespace ActsFatras