// This file is part of the Acts project.
//
// Copyright (C) 2016-2019 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include <cmath>
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "Acts/Utilities/ParameterDefinitions.hpp"

#ifdef ACTS_COORDINATE_TRANSFORM_PLUGIN

#include ACTS_COORDINATE_TRANSFORM_PLUGIN

#else

namespace Acts {
/// @cond detail
namespace detail {

/// @brief helper structure summarizing coordinate transformations
///
struct coordinate_transformation {
  using ParVector_t = BoundVector;

  /// @brief static method to transform the local information in
  /// the track parameterisation to a global position
  ///
  /// This transformation uses the surface and hence needs a context
  /// object to guarantee the local to global transformation is done
  /// within the right (alginment/conditions) context
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param pars the parameter vector
  /// @param s the surface for the local to global transform
  ///
  /// @return position in the global frame
  static Vector3D parameters2globalPosition(const GeometryContext& gctx,
                                            const ParVector_t& pars,
                                            const Surface& s) {
    Vector3D globalPosition;
    s.localToGlobal(gctx, Vector2D(pars(Acts::eLOC_0), pars(Acts::eLOC_1)),
                    parameters2globalMomentum(pars), globalPosition);
    return globalPosition;
  }

  /// @brief static method to transform the momentum parameterisation
  /// into a global momentum vector
  ///
  /// This transformation does not use the surface and hence no context
  /// object is needed
  ///
  /// @param pars the parameter vector
  ///
  /// @return momentum in the global frame
  static Vector3D parameters2globalMomentum(const ParVector_t& pars) {
    Vector3D momentum;
    double p = std::abs(1. / pars(Acts::eQOP));
    double phi = pars(Acts::ePHI);
    double theta = pars(Acts::eTHETA);
    momentum << p * sin(theta) * cos(phi), p * sin(theta) * sin(phi),
        p * cos(theta);

    return momentum;
  }

  /// @brief static method to transform a global representation into
  /// a curvilinear represenation
  ///
  /// This transformation does not use the surface and hence no context
  /// object is needed
  ///
  /// @param pos - ignored
  /// @param mom the global momentum parameters
  /// @param charge of the particle/track
  ///
  /// @return curvilinear parameter representation
  static ParVector_t global2curvilinear(const Vector3D& /*pos*/,
                                        const Vector3D& mom, double charge,
                                        double time) {
    using VectorHelpers::phi;
    using VectorHelpers::theta;
    ParVector_t parameters;
    parameters << 0, 0, phi(mom), theta(mom),
        ((std::abs(charge) < 1e-4) ? 1. : charge) / mom.norm(), time;

    return parameters;
  }

  /// @brief static method to transform the global information into
  /// the track parameterisation
  ///
  /// This transformation uses the surface and hence needs a context
  /// object to guarantee the local to global transformation is done
  /// within the right (alginment/conditions) context
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  /// @param pos position of the parameterisation in global
  /// @param mom position of the parameterisation in global
  /// @param charge of the particle/track
  /// @param s the surface for the global to local transform
  ///
  /// @return the track parameterisation
  static ParVector_t global2parameters(const GeometryContext& gctx,
                                       const Vector3D& pos, const Vector3D& mom,
                                       double charge, double time,
                                       const Surface& s) {
    using VectorHelpers::phi;
    using VectorHelpers::theta;
    Vector2D localPosition;
    s.globalToLocal(gctx, pos, mom, localPosition);
    ParVector_t result;
    result << localPosition(0), localPosition(1), phi(mom), theta(mom),
        ((std::abs(charge) < 1e-4) ? 1. : charge) / mom.norm(), time;
    return result;
  }

  /// @brief static calculate the charge from the track parameterisation
  ///
  /// @return the charge as a double
  static double parameters2charge(const ParVector_t& pars) {
    return (pars(Acts::eQOP) > 0) ? 1. : -1.;
  }
};
}  // namespace detail
/// @endcond
}  // namespace Acts

#endif  // ACTS_COORDINATE_TRANSFORM_PLUGIN
