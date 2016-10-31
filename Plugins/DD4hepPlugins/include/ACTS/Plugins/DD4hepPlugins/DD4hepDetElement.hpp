// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// DD4hepDetElement.h, ACTS project, DD4hepDetector plugin
///////////////////////////////////////////////////////////////////

#ifndef ACTS_DD4HEPDETECTORELEMENT_DD4HEPDETELEMENT_H
#define ACTS_DD4HEPDETECTORELEMENT_DD4HEPDETELEMENT_H 1

#include "ACTS/Plugins/TGeoPlugins/TGeoDetectorElement.hpp"
#include "DD4hep/Detector.h"

namespace Acts {

/// @class DD4hepDetElement
///
/// @brief DetectorElement class implementation for DD4hep geometry
///
/// DetectorElement plugin for DD4hep detector elements. DD4hep is based on TGeo
/// shapes, therefore the DD4hepDetElement inherits from TGeoDetectorElement.
/// The
/// full geometrical information is provided by the TGeoDetectorElement. The
/// DD4hepDetrElement ectends the TGeoDetectorElement by containing a
/// segmentation
/// for the readout.

/// @todo what if shape conversion failes? add implementation of more than one
/// surface per module, implementing also for other shapes->Cone,ConeSeg,Tube?
/// what
/// if not used with DD4hep?
/// @todo segmentation

class DD4hepDetElement : public TGeoDetectorElement
{
public:
  /// Constructor
  /// @param detElement The DD4hep DetElement which should be linked to a
  /// surface
  DD4hepDetElement(const DD4hep::Geometry::DetElement detElement,
                   const std::string&                 axes   = "XYZ",
                   double                             scalor = 1.);
  /// Desctructor
  virtual ~DD4hepDetElement() = default;

private:
  /// DD4hep detector element
  DD4hep::Geometry::DetElement m_detElement;
  /// DD4hep segmentation
  DD4hep::Geometry::Segmentation m_segmentation;
};
}

#endif  // ACTS_DD4HEPDETECTORELEMENT_DD4HEPDETELEMENT_H
