// This file is part of the Acts project.
//
// Copyright (C) 2016-2017 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ApproachDescriptor.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

#include <vector>
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/Intersection.hpp"

namespace Acts {

class Surface;
class Layer;
class BoundaryCheck;

/// @class ApproachDescriptor
///
/// Virtual base class to decide and return which approaching surface to be
/// taken,
/// the surfaces are std::shared_ptr, as they can be the boundary surfaces of
/// the
/// representingVolume of the Layer
///
///
class ApproachDescriptor
{
public:
  /// Default constructor
  ApproachDescriptor() {}

  /// Virtual destructor
  virtual ~ApproachDescriptor() {}

  /// @brief Register Layer
  /// this gives the approach surfaces the link to the layer
  ///
  /// @param lay is the layer to be assigned
  virtual void
  registerLayer(const Layer& lay)
      = 0;

  /// @brief Get the surface on approach
  ///
  /// @tparam parameters_t Type of the Parameters for the approach search
  /// @tparam options_t Type of the Navigation options for the search
  /// @tparam corrector_t Type of the Corrector for the approach search
  ///
  /// @param parameters The actual parameters object
  /// @param options are the steering options for the search
  template <typename parameters_t,
            typename options_t,
            typename corrector_t = VoidCorrector>
  ObjectIntersection<Surface>
  approachSurface(const parameters_t& parameters,
                  const options_t&    options,
                  const corrector_t&  correct = corrector_t()) const;

  /// @brief Get the surface on approach
  ///
  /// @param pos is the position from start of the search
  /// @param mom is the momentum
  /// @param bcheck is the boundary check directive
  /// @param correciton is the function pointer to a corrector
  ///
  /// @return is a surface intersection
  virtual ObjectIntersection<Surface>
  approachSurface(const Vector3D&      pos,
                  const Vector3D&      mom,
                  NavigationDirection  navDir,
                  const BoundaryCheck& bcheck,
                  CorrFnc              correct = nullptr) const = 0;

  /// Tet to all the contained surfaces
  /// @return all contained surfaces of this approach descriptor
  virtual const std::vector<const Surface*>&
  containedSurfaces() const = 0;

  /// Non-const version
  virtual std::vector<const Surface*>&
  containedSurfaces()
      = 0;
};

template <typename parameters_t, typename options_t, typename corrector_t>
ObjectIntersection<Surface>
ApproachDescriptor::approachSurface(const parameters_t& parameters,
                                    const options_t&    options,
                                    const corrector_t&  corrfnc) const
{
  // calculate the actual intersection
  return approachSurface(parameters.position(),
                         parameters.momentum(),
                         options.navDir,
                         options.boundaryCheck,
                         corrfnc);
}

}  // end of namespace Acts
