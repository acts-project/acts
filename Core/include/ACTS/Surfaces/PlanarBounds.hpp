// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// PlanarBounds.h, ACTS project
///////////////////////////////////////////////////////////////////

#ifndef ACTS_SURFACES_PLANARBOUNDS_H
#define ACTS_SURFACES_PLANARBOUNDS_H 1

#include "ACTS/Surfaces/SurfaceBounds.hpp"

namespace Acts {

  /**
   @class PlanarBounds
   
   common base class for all bounds that are in a local x/y cartesian frame
    - simply introduced to avoid wrong bound assigments to surfaces
    
   @author Andreas.Salzburger@cern.ch
   */
      
  class PlanarBounds : public SurfaceBounds {

    public:

      /** Default Constructor */
      PlanarBounds() : SurfaceBounds() {}
      
      /** Destructor */
      virtual ~PlanarBounds(){}
      
      /** Virtual Constructor */
      virtual PlanarBounds* clone() const = 0;
      
      /** Return the vertices - or, the points of the extremas */
      virtual const std::vector< Vector2D > vertices() const = 0;
      
  };
                              
} // end of namespace

#endif // ACTS_SURFACES_PLANARBOUNDS_H
