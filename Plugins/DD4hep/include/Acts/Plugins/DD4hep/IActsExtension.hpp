// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// IDetExtension.h, Acts project, DD4hepDetector plugin
///////////////////////////////////////////////////////////////////

#pragma once
// Algebra
#include <memory>
#include <vector>
#include "Acts/Utilities/Definitions.hpp"

namespace dd4hep {
class DetElement;
class Segmentation;
class Volume;
class Material;
}  // namespace dd4hep

namespace Acts {

/// forwared declaration of DigitzationModule is enough
class DigitizationModule;

/// @class IActsExtension
///
/// @brief Interface for Extension of the DD4hep \a DetElement needed for
/// translation into Acts tracking geometry
///
/// Interface class for making extensions to the DD4hep \a DetElement class,
/// needed  for the translation from the DD4hep geometry into the tracking
/// geometry of the Acts package.
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
/// be declared as a layer.
/// However, the layer does not need to be the direct mother
/// of these sensitive DD4hep \a DetElements - they can also be nested in other
/// DD4hep \a DetElement substructures. Moreover every DD4hep \a DetElement
/// layer which should carry material should also be declared as layer and the
/// specific parameters needed for material mapping should be handed over.
/// In case the sensitive modules/components contained by a layer have a
/// different orientation in respect to the local tracking frame of Acts, the
/// axes orientation of these modules can be set for the layer.
/// In case the segmentation is set for the module (by using the second
/// constructor setting the segmentation or using the setSegmentation()
/// function), the axes do not need to be set for the layer, which containes
/// the modules.
///
/// In DD4hep cylinder and disc volumes are both described with the underlying
/// ROOT TGeoConeSeg class. In Acts one needs to distinguish these two volume
/// types. Therefore volumes which are endcaps or barrels should be indicated as
/// these.
/// The ActsExtension should also be used to indicate that a \a DetElement is
/// the
/// beampipe.
///
/// @note see Acts::ActsExtension for implementation.

class ISurfaceMaterial;

/// @enum LayerMaterialPos The LayerMaterialPos enumeration is foreseen to
/// mark
/// on which surface
/// the two dimensional material of the layer sits. Either on the inner or
/// the other boundary surface of the layer or at the representing (center)
/// surface of the layer
enum LayerMaterialPos {
  inner = 0,
  central = 1,
  outer = 2,
};

class IActsExtension {
 public:
  /// Virtual destructor
  virtual ~IActsExtension() = default;
  /// Indicates if the DD4hep::DetElement is the beampipe
  virtual bool isBeampipe() const = 0;
  /// Indicates that the DD4hep::DetElement is a sensitive barrel
  virtual bool isBarrel() const = 0;
  /// Indicates that the DD4hep::DetElement is an sensitive endcap
  virtual bool isEndcap() const = 0;
  /// Indicates that the DD4hep::DetElement is a layer
  virtual bool isLayer() const = 0;
  /// Indicates that the DD4hep::DetElement is a passive cylinder layer
  virtual bool isPassiveCylinder() const = 0;
  /// Bool returning true if the layers should carry material using material
  /// mapping
  /// @note automatically set when the material bins are set
  virtual bool hasSupportMaterial() const = 0;
  /// Access to the two bin numbers determining the granularity of the two
  /// dimensional grid on which the material of the layer should be mapped on
  /// @return std::pair with the number of bins in th first and the second
  /// direction
  virtual std::pair<size_t, size_t> materialBins() const = 0;
  /// @return states if the material should be mapped on the inner,
  /// the center or the outer surface of the layer
  virtual Acts::LayerMaterialPos layerMaterialPosition() const = 0;
  /// Access the orientation of the module in respect to the tracking frame
  /// @return string describing the orientation of the axes
  virtual const std::string axes() const = 0;
  /// @return states if the geometrical boundaries of the current object should
  /// be built automatically by adding given tolerances to the expansion of the
  /// contained modules
  virtual bool buildEnvelope() const = 0;
  /// @return the tolerance which should be added in r to the geometrical
  /// expansion of the contained surfaces (sensituive DetElements) of this
  /// DetElement to automatically create the layer envelope
  virtual double envelopeR() const = 0;
  /// @return the tolerance which should be added in z to the geometrical
  /// expansion of the contained surfaces (sensituive DetElements) of this
  /// DetElement to automatically create the layer envelope
  virtual double envelopeZ() const = 0;
  /// @return The SurfaceMaterial
  virtual std::shared_ptr<const Acts::ISurfaceMaterial> material() const = 0;
  /// @return the shared pointer to the digitization module
  virtual std::shared_ptr<const DigitizationModule> digitizationModule()
      const = 0;

 protected:
  /// Protected constructor
  IActsExtension() = default;
};
}  // namespace Acts
