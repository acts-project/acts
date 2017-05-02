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

#ifndef ACTS_DD4HEPDETECTORELEMENT_IACTSEXTENSION_H
#define ACTS_DD4HEPDETECTORELEMENT_IACTSEXTENSION_H 1

// Algebra
#include <memory>
#include <vector>
#include "ACTS/Utilities/Definitions.hpp"

namespace DD4hep {
namespace Geometry {
  class DetElement;
  class Segmentation;
  class Volume;
}
}

namespace Acts {

/// @class IActsExtension
///
/// @brief Interface for Extension of the DD4hep \a DetElement needed for
/// translation into ACTS tracking geometry
///
/// Interface class for making extensions to the DD4hep \a DetElement class,
/// needed  for the translation from the DD4hep geometry into the tracking
/// geometry of the ACTS package.
///
/// This extensions are necessary in order to distinguish in the translation if
/// a DD4hep \a DetElement is
/// 	- the beampipe
///		- a barrel volume
///		- an endcap volume
/// 	- a layer
/// and to hand over needed parameters.
///
/// Every DD4hep \a DetElement containing sensitive DD4hep \a DetElements has to
/// be
/// declared as a layer. However the layer does not need to be the direct mother
/// of these sensitive DD4hep \a DetElements - they can also be nested in other
/// DD4hep \a DetElement substructures. Moreover every DD4hep \a DetElement
/// layer
/// which should carry material should also be declared as a layer and the
/// specific parameters needed for material mapping should be handed over.
/// In case the sensitive modules/components contained by a layer have a
/// different orientation in respect to the local tracking frame of ACTS, the
/// axes orientation of these modules can be set for the layer.
/// In case the segmentation is set for the module (by using the second
/// constructor setting the segmentation or using the setSegmentation()
/// function), the axes do not need to be set for the layer, which containes
/// the modules.
///
/// In DD4hep cylinder and disc volumes are both described with the underlying
/// ROOT TGeoConeSeg class. In ACTS one needs to distinguish these two volume
/// types. Therefore volumes which are endcaps or barrels should be indicated as
/// these.
/// The ActsExtension should also be used to indicate that a \a DetElement is
/// the
/// beampipe.
///
/// @note see Acts::ActsExtension for implementation.

class DigitizationModule;

/// @enum LayerMaterialPos The LayerMaterialPos enumeration is foreseen to mark
/// on which surface
/// the two dimensional material of the layer sits. Either on the inner or
/// the other boundary surface of the layer or at the representing (center)
/// surface of the layer
enum LayerMaterialPos {
  inner   = 0,
  central = 1,
  outer   = 2,
};

class IActsExtension
{
public:
  /// Virtual destructor
  virtual ~IActsExtension() {}
  /// Indicates if the DD4hep::DetElement is the beampipe
  virtual bool
  isBeampipe() const = 0;
  /// Indicates that the DD4hep::DetElement is a barrel
  virtual bool
  isBarrel() const = 0;
  /// Indicates that the DD4hep::DetElement is an endcap
  virtual bool
  isEndcap() const = 0;
  /// Indicates that the DD4hep::DetElement is a layer
  virtual bool
  isLayer() const = 0;
  /// Bool returning true if the layers should carry material
  /// @note automatically set when the material bins are set
  virtual bool
  hasSupportMaterial() const = 0;
  /// Access to the two bin numbers determining the granularity of the two
  /// dimensional grid on which the material of the layer should be mapped on
  /// @return std::pair with the number of bins in th first and the second
  /// direction
  virtual std::pair<size_t, size_t>
  materialBins() const = 0;
  /// @return states if the material should be mapped on the inner,
  /// the center or the outer surface of the layer
  virtual Acts::LayerMaterialPos
  layerMaterialPosition() const = 0;
  /// Access the orientation of the module in respect to the tracking frame
  /// @return string describing the orientation of the axes
  virtual const std::string
  axes() const = 0;
  /// @return states if the geometrical boundaries of the current object should
  /// be built automatically by adding given tolerances to the expansion of the
  /// contained modules
  virtual bool
  buildEnvelope() const = 0;
  /// @return the tolerance which should be added in r to the geometrical
  /// expansion of the contained surfaces (sensituive DetElements) of this
  /// DetElement to automatically create the layer envelope
  virtual double
  envelopeR() const = 0;
  /// @return the tolerance which should be added in z to the geometrical
  /// expansion of the contained surfaces (sensituive DetElements) of this
  /// DetElement to automatically create the layer envelope
  virtual double
  envelopeZ() const = 0;
  /// In case several sensitive modules have the same segmentation function can
  /// to set the segmentation of this extension, which can than be attached to
  /// the different modules. In this way the Acts::DigitizationModule will
  /// be shared amongst these modules which saves memory.
  /// In case this function is used for the module, the axes do not need
  /// to be set for the layer, which containes this module.
  /// @param segmentation The DD4hep segmentation object
  /// @param volume The DD4hep logical volume
  /// @param axes The orientation of the axes in respect to the tracking frame
  ///  A different orientation can occur because in TGeo (which
  /// is the underlying geometry model of %DD4hep) all shapes are 3D volumes
  /// also the sensitive components of a detector. In the ACTS tracking
  /// geometry these sensitive components are described as 2D surfaces, which
  /// have their local 2D coordinate system. Therefore one needs to know which
  /// coordinates should be taken as the local coordinates.
  /// A string of the three characters x, y and z (standing for the
  /// three axes) needs to be handed over. There is a distinction between
  /// capital and lower case
  /// characters :
  /// 	- capital      -> positive orientation of the axis
  ///		- lower case   -> negative oriantation of the axis
  ///
  ///
  /// Example options are:
  /// 	- "XYZ" -> identical frame definition (default value)
  /// 	- "YZX" -> node y axis is tracking x axis, etc.
  ///		- "XzY" -> negative node z axis is tracking y axis, etc.
  virtual void
  setSegmentation(DD4hep::Geometry::Segmentation segmentation,
                  DD4hep::Geometry::Volume       volume,
                  std::string                    axes = "XYZ")
      = 0;
  /// @return the Acts::DigitizationModule
  virtual std::shared_ptr<const DigitizationModule>
  digitizationModule() const = 0;

protected:
  /// Protected constructor
  IActsExtension() {}
};
}

#endif  // ACTS_DD4HEPDETECTORELEMENT_DET_IACTSEXTENSION_H
