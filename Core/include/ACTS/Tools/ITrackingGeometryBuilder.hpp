// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// ITrackingGeometryBuilder.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_GEOMETRYINTERFACES_ITRACKINGGEOMETRYBUILDER_H
#define ACTS_GEOMETRYINTERFACES_ITRACKINGGEOMETRYBUILDER_H 1

#include<memory>

namespace Acts
{
  class TrackingGeometry;

  /** @class ITrackingGeometryBuilder
    
    Interface class for the TrackingGeometry building,
    this is used by the TrackingGeometrySvc to build the geoemtry.
  
    The TrackingGeometry is written to the detector store and thus not created
    as a std::shared_ptr.
  
    The TrackingGeometry is returned as a non-const object in order to recreate
    from conditions callback if necessary.
      
    */
  class ITrackingGeometryBuilder
  {
  public:
    /**Virtual destructor*/
    virtual ~ITrackingGeometryBuilder() = default;
      
    /** TrackingGeometry Interface methode */
    virtual std::unique_ptr<TrackingGeometry> trackingGeometry() const = 0;      
  };
} // end of namespace

#endif // ACTS_GEOMETRYINTERFACES_IGEOMETRYBUILDER_H
