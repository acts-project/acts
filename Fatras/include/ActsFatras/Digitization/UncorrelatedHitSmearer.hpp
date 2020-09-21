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
#include "ActsFatras/Digitization/detail/ParametersSmearer.hpp"
#include "ActsFatras/EventData/Hit.hpp"

#include <functional>

namespace ActsFatras {

/// Smearing functions definition:
/// - it takes the unsmeared parameter
/// - it returns the smeared parameter and a covariance
using SmearFunction =
    std::function<Acts::Result<std::pair<double, double>>(double)>;

/// Smearing input to be used by the smearers
/// - this struct helps to harmonize the interface between
///   free and bound smearers
struct SmearInput {
  /// Only valid constructor, wraps the @param hit_,
  /// the  and optionally the @param surface_
  SmearInput(std::reference_wrapper<const Hit> hit_,
             std::reference_wrapper<const Acts::GeometryContext> geoContext_,
             std::shared_ptr<const Acts::Surface> surface_ = nullptr)
      : hit(hit_), geoContext(geoContext_), surface(surface_) {}

  SmearInput() = delete;

  std::reference_wrapper<const Hit> hit;
  std::reference_wrapper<const Acts::GeometryContext> geoContext;
  std::shared_ptr<const Acts::Surface> surface = nullptr;
};

/// Parameter smearer for fast digitisation for bound parameters
///
///
/// The logic of this smearer is the following:
/// - Input is a single simulated ActsFatras::Hit which is wrapped into
/// - It returns smeared parameter vector and a covariance matrix
///
/// @note This smearer only supports uncorrelated smearing of parameters
///
template <Acts::BoundIndices... kParameters>
struct BoundParametersSmearer {
  /// Generic implementation of a smearing meathod for bound parameters
  ///
  /// @tparam kParameters parameter pack describing the parameters to smear
  ///
  /// @param sInput The smearing input struct: surface and simulated hit
  /// @param sFunctions The smearing functions that are applied
  ///
  /// @return Smeared bound parameter set wrapped in a Result<...> object
  Acts::Result<Acts::ParameterSet<Acts::BoundIndices, kParameters...>>
  operator()(const SmearInput& sInput,
             const std::array<SmearFunction, sizeof...(kParameters)>&
                 sFunctions) const {
    using ParSet = Acts::ParameterSet<Acts::BoundIndices, kParameters...>;
    using Result = Acts::Result<ParSet>;
    using ParametersSmearer =
        detail::ParametersSmearer<Acts::BoundIndices, kParameters...>;

    if (sInput.surface == nullptr) {
      return Result(ActsFatras::DigitizationError::NoSurfaceDefined);
    }

    const auto& hit = sInput.hit.get();

    auto dir = hit.unitDirection();
    auto gltResult =
        sInput.surface->globalToLocal(sInput.geoContext, hit.position(), dir);
    if (not gltResult.ok()) {
      return Result(Acts::SurfaceError::GlobalPositionNotOnSurface);
    }
    const auto& lPosition = gltResult.value();

    typename ParSet::FullParametersVector fParameters;
    fParameters.setZero();
    fParameters.template segment<2>(0) = lPosition;
    fParameters[Acts::eBoundPhi] = Acts::VectorHelpers::phi(dir);
    fParameters[Acts::eBoundTheta] = Acts::VectorHelpers::theta(dir);
    fParameters[Acts::eBoundTime] = hit.time();
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
};

/// Generic implementation of a smearing meathod for free parameters
///
/// @tparam kParameters parameter pack describing the parameters to smear
template <Acts::FreeIndices... kParameters>
struct FreeParametersSmearer {
  /// Smearing function
  ///
  /// @param sInput The smearing input struct with the simulated hit
  /// @param sFunctions The smearing functions that are applied
  ///
  /// @note uncorrelated smearing of the direction using the components
  ///       is not recommended
  ///
  /// @return Smeared free parameter set wrapped in a Result<...> object
  Acts::Result<Acts::ParameterSet<Acts::FreeIndices, kParameters...>>
  operator()(const SmearInput& sInput,
             const std::array<SmearFunction, sizeof...(kParameters)>&
                 sFunctions) const {
    using ParSet = Acts::ParameterSet<Acts::FreeIndices, kParameters...>;
    using Result = Acts::Result<ParSet>;
    using ParametersSmearer =
        detail::ParametersSmearer<Acts::FreeIndices, kParameters...>;

    const auto& hit = sInput.hit.get();

    typename ParSet::FullParametersVector fParameters;
    fParameters.setZero();
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
};

}  // namespace ActsFatras
