// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ILayerCreator.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYINTERFACES_ILAYERCREATOR_H
#define ACTS_GEOMETRYINTERFACES_ILAYERCREATOR_H 1

// STL
#include <memory>
#include <string>
#include <vector>

namespace Acts {

class Surface;
class Layer;
typedef std::shared_ptr<const Layer> LayerPtr;

/// @class ILayerCreator
///
/// Interface class for LayerCreator from DetectorElements
///
class ILayerCreator
{
public:
  /// Virtual destructor
  virtual ~ILayerCreator() {}
  /// ILayerCreator interface method - returning a cylindrical layer
  /// @param surfaces is the vector of sensitive surfaces represented by this
  /// layer
  /// @param envelopeR is the additional envelope applied in R
  /// @param envelopeZ is the additional envelope applied in z
  /// @param binsRPhi is number of bins the sensitive surfaces are ordered in
  /// phi
  /// @param binsZ is number of bins the sensitive surfaces are ordered in Z
  virtual LayerPtr
  cylinderLayer(const std::vector<const Surface*>& surfaces,
                double                             envelopeR,
                double                             envelopeZ,
                size_t                             binsRPhi,
                size_t                             binsZ) const = 0;

  /// ILayerCreator interface method - returning a disc layer
  /// @param surfaces is the vector of sensitive surfaces represented by this
  /// layer
  /// @param envelopeMinR is the additional envelope applied in R at Rmin
  /// @param envelopeMaxR is the additional envelope applied in R in Rmax
  /// @param envelopeZ is the additional envelope applied in z
  /// @param binsR is number of bins the sensitive surfaces are ordered in R
  /// @param binsPhi is number of bins the sensitive surfaces are ordered in Phi
  virtual LayerPtr
  discLayer(const std::vector<const Surface*>& surfaces,
            double                             envelopeMinR,
            double                             envelopeMaxR,
            double                             envelopeZ,
            size_t                             binsR,
            size_t                             binsPhi) const = 0;

  /// ILayerCreator interface method - returning a plane layer
  /// @param surfaces is the vector of sensitive surfaces represented by this
  /// layer
  /// @param envelopeXY is the additional envelope applied in XY
  /// @param envelopeZ is the additional envelope applied in Z
  /// @param binsX is number of bins the sensitive surfaces are ordered in X
  /// @param binsY is number of bins the sensitive surfaces are ordered in Y
  virtual LayerPtr
  planeLayer(const std::vector<const Surface*>& surfaces,
             double                             envelopeXY,
             double                             envelopeZ,
             size_t                             binsX,
             size_t                             binsY) const = 0;
};

}  // end of namespace

#endif  // ACTS_GEOMETRYINTERFACES_ILAYERCREATOR_H
