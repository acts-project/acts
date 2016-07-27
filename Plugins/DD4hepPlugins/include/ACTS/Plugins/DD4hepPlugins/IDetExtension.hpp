// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// IDetExtension.h, ACTS project, DD4hepDetector plugin
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DD4HEPDETECTORELEMENT_IDETEXTENSION_H
#define ACTS_DD4HEPDETECTORELEMENT_IDETEXTENSION_H 1

// Algebra
#include "ACTS/Utilities/Definitions.hpp"
#include <vector>

namespace DD4hep {
namespace Geometry {
  class DetElement;
  class Segmentation;
}
}

namespace Acts {

/// @class IDetExtension
///
/// @brief Interface class for making extensions to the DD4hep::DetElement class
///
/// Interface class for making extensions to the DD4hep::DetElement class,
/// needed
/// for the translation from the DD4hep geometry into the tracking geometry of
/// the
/// ATS package.
/// In this way, the segmentation of the sensitive detector elements can be
/// directly accessed from DD4hep to ensure consistency between the full and the
/// tracking geometry.
/// Since in DD4hep volumes used as a cylinder (detector layers are binned in r
/// and
/// z, e.g. central barrel volume) and discs (detector layers are binned in r
/// and
/// phi, e.g. end caps) are both described as a ROOT TGeoConeSeg one needs to
/// distinguish between these volume types by setting the shape.
/// @TODO find replacement for Gaudi exeption and message stream

class Module;

    /// @enum The ShapeType enumeration is used to distinguish between a volume
    /// having a disc and a volume having a cylinder shape, since in DD4hep they
    /// are both described as TGeoConeSeg.
enum ShapeType {

  None     = 0,
  Cylinder = 1,
  Disc     = 2
};
    /// @enum The LayerMaterialPos enumeration is foreseen to mark on which surface
    /// the two dimensional material of the layer sits. Either on the inner or
    /// the other boundary surface of the layer or at the representing (center)
    /// surface of the layer
enum LayerMaterialPos{
        inner   = 0,
        central = 1,
        outer   = 2,
};

class IDetExtension
{
public:
  /// Virtual destructor
  virtual ~IDetExtension() {}
  /// Possibility to set shape of a volume to distinguish between disc and
  /// cylinder volume
  /// @param shape The type of the shape defined in IDetExtension can be either
  /// disc or cylinder
  virtual void
  setShape(ShapeType shape)
      = 0;
  /// Access shape type of a volume to distinguish between disc and cylinder
  /// volume
  /// @return shape The type of the shape defined in IDetExtension can be either
  /// disc or cylinder
  virtual ShapeType
  shape() const = 0;
  /// Method to set the DD4hep segmentation for the readout
  /// @param segmentation DD4hep segmentation for the readout
  virtual void
  setSegmentation(const DD4hep::Geometry::Segmentation segmentation)
      = 0;
  /// Method to access the DD4hep segmentation for the readout
  /// @return segmentation DD4hep segmentation for the readout
  virtual const DD4hep::Geometry::Segmentation
  segmentation() const = 0;
  /// possibility to mark layer to have support material
  /// @note the numer of bins determines the granularity of the material
  /// map of the layer
  /// @param bins1 The number of bins in first direction of the layer
  /// which is phi for both, cylinder and disc layers.
  /// @param bins2 The number of bins in second direction of the layer
  /// which is r in case of a disc layer and z in case of a cylinder layer
  /// @param layerMatPos States if the material should be mapped on the inner,
  /// the center or the outer surface of the layer
  virtual void
  supportMaterial(size_t bins1, size_t bins2, LayerMaterialPos layerMatPos)
      = 0;
  /// Bool returning true if the layers should carry material
  virtual bool
  hasSupportMaterial() const = 0;
  /// Access to the two bin numbers determining the granularity of the two dimensional grid
  /// on which the material of the layer should be mapped on
  virtual std::pair<size_t,size_t> materialBins() const = 0;
  /// returns states if the material should be mapped on the inner,
  /// the center or the outer surface of the layer
  virtual Acts::LayerMaterialPos layerMaterialPos() const = 0;
  /// Possibility to set contained detector modules of a layer
  /// @param mod Possible sensitive modules contained by a layer
  virtual void
  setModules(std::vector<DD4hep::Geometry::DetElement> mod)
      = 0;
  /// Access modules detector module contained by a layer
  /// @return mod Possible sensitive modules contained by a layer
  virtual std::vector<DD4hep::Geometry::DetElement>
  modules() const = 0;

protected:
  /// Protected constructor
  IDetExtension() {}
};
}

#endif  // ACTS_DD4HEPDETECTORELEMENT_DET_IDETEXTENSION_H
