// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DetExtension.h, ACTS project, DD4hepDetector plugin
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DD4HEPDETECTORELEMENT_DETEXTENSION_H
#define ACTS_DD4HEPDETECTORELEMENT_DETEXTENSION_H 1

#include "ACTS/Plugins/DD4hepPlugins/IDetExtension.hpp"
// DD4hep
#include "DD4hep/Detector.h"

namespace Acts {

/// @class DetExtension
///
/// Implementation of the IDetExtension class, which uses the extension
/// mechanism
/// of DD4hep, needed for the translation from the DD4hep geometry into the
/// tracking geometry of the ACTS package.
/// In this way, the segmentation of the sensitive detector elements can be
/// directly accessed from DD4hep to ensure consistency between the full and the
/// tracking geometry.
/// Since in DD4hep volumes used as a cylinder (detector layers are binned in r
/// and
/// z, e.g. central barrel volume) and discs (detector layers are binned in r
/// and
/// phi, e.g. end caps) are both described as a ROOT TGeoConeSeg one needs to
/// distinguish between these volume types by setting the shape.

class DetExtension : virtual public IDetExtension
{
public:
  /// Constructor
  DetExtension();
  /// Constructor for sensitive detector element
  /// @param segmentation DD4hep segmentation for the readout
  DetExtension(const DD4hep::Geometry::Segmentation segmentation);
  /// Constructor for volume with shape to distinguish between disc and cylinder
  /// volume
  /// @param shape The type of the shape defined in IDetExtension can be either
  /// disc or cylinder
  DetExtension(ShapeType shape);
  /// Constructor for layer with support structure
  /// @param support Possible support structure of the layer
  DetExtension(const DD4hep::Geometry::DetElement support);
  /// Constructor for layer with modules
  /// @param mod Possible sensitive modules contained by a layer
  DetExtension(std::vector<DD4hep::Geometry::DetElement> mod);
  /// Constructor for layer with support structure and modules
  /// @param support Possible support structure of the layer
  /// @param mod Possible sensitive modules contained by a layer
  DetExtension(const DD4hep::Geometry::DetElement        support,
               std::vector<DD4hep::Geometry::DetElement> mod);
  /// Copy constructor
  DetExtension(const DetExtension&, const DD4hep::Geometry::DetElement&);
  /// Virtual destructor
  virtual ~DetExtension() = default;
  /// Possibility to set shape of a volume to distinguish between disc and
  /// cylinder volume
  /// @param shape The type of the shape defined in IDetExtension can be either
  /// disc or cylinder
  void
  setShape(ShapeType shape) override;
  /// Access shape type of a volume to distinguish between disc and cylinder
  /// volume
  /// @return shape The type of the shape defined in IDetExtension can be either
  /// disc or cylinder
  ShapeType
  shape() const override;
  /// Method to set the DD4hep segmentation for the readout
  /// @param segmentation DD4hep segmentation for the readout
  void
  setSegmentation(const DD4hep::Geometry::Segmentation segmentation) override;
  /// Method to access the DD4hep segmentation for the readout
  /// @return segmentation DD4hep segmentation for the readout
  const DD4hep::Geometry::Segmentation
  segmentation() const override;
  /// possibility to hand over supporte structure of a layer
  /// @param support Possible support structure of the layer
  void
  setSupportStructure(const DD4hep::Geometry::DetElement support) override;
  /// Access supporting structure of a layer
  /// @return support Possible support structure of the layer
  const DD4hep::Geometry::DetElement&
  supportStructure() const override;
  /// Possibility to set contained detector modules of a layer
  /// @param mod Possible sensitive modules contained by a layer
  void
  setModules(std::vector<DD4hep::Geometry::DetElement> mod) override;
  /// Access modules detector module contained by a layer
  /// @return mod Possible sensitive modules contained by a layer
  std::vector<DD4hep::Geometry::DetElement>
  modules() const override;

private:
  /// segmentation of a sensitive detector module
  DD4hep::Geometry::Segmentation m_segmentation;
  /// shape type of a volume defined in IDetExtension can be either disc or
  /// cylinder
  ShapeType m_shape;
  /// possible support structure of a layer
  DD4hep::Geometry::DetElement m_supportStructure;
  /// possible contained modules of a layer
  std::vector<DD4hep::Geometry::DetElement> m_modules;
};
}

inline void
Acts::DetExtension::setShape(Acts::ShapeType type)
{
  m_shape = std::move(type);
}

inline Acts::ShapeType
Acts::DetExtension::shape() const
{
  return m_shape;
}

inline void
Acts::DetExtension::setSegmentation(const DD4hep::Geometry::Segmentation seg)
{
  m_segmentation = std::move(seg);
}

inline const DD4hep::Geometry::Segmentation
Acts::DetExtension::segmentation() const
{
  return m_segmentation;
}

inline void
Acts::DetExtension::setSupportStructure(
    const DD4hep::Geometry::DetElement support)
{
  m_supportStructure = std::move(support);
}

inline const DD4hep::Geometry::DetElement&
Acts::DetExtension::supportStructure() const
{
  return m_supportStructure;
}

inline void
Acts::DetExtension::setModules(std::vector<DD4hep::Geometry::DetElement> mod)
{
  m_modules = std::move(mod);
}

inline std::vector<DD4hep::Geometry::DetElement>
Acts::DetExtension::modules() const
{
  return m_modules;
}

#endif  // ACTS_DD4HEPDETECTORELEMENT_DETEXTENSION_H
