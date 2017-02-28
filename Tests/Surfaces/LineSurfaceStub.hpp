// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef ACTS_TESTS_SURFACES_LINESURFACESTUB
#define ACTS_TESTS_SURFACES_LINESURFACESTUB 1

//
#include "ACTS/Surfaces/LineSurface.hpp"
#include "ACTS/Utilities/Definitions.hpp"
//
//
#include <limits>


namespace Acts {
class LineSurfaceStub: public LineSurface{
  public:
  LineSurfaceStub() = delete;
  //
  LineSurfaceStub(std::shared_ptr<Transform3D> htrans, double radius, double halfz):
    LineSurface(htrans, radius, halfz){ /* nop */}
  //
  LineSurfaceStub(std::shared_ptr<Transform3D> htrans, 
    std::shared_ptr<const LineBounds> lbounds=nullptr):
    LineSurface(htrans, lbounds){ /*nop */ }
  //
  LineSurfaceStub(std::shared_ptr<const LineBounds> lbounds, 
    const DetectorElementBase & detelement,
    const Identifier & identifier=Identifier()): 
    LineSurface(lbounds, detelement, identifier){ /* nop */ }
  //
  LineSurfaceStub(const LineSurfaceStub & ls):LineSurface(ls){ /* nop */}
  //
  LineSurfaceStub(const LineSurfaceStub & ls, const Transform3D &t):
    LineSurface(ls,t){ /* nop */}
  ///pure virtual functions of baseclass implemented here
  Surface*
  clone(const Transform3D* shift = nullptr) const final{
    return nullptr;
  }
  /// Return method for the Surface type to avoid dynamic casts
  SurfaceType
  type() const final
  {
    return Surface::Straw;
  }
  
  /// Simply return true to show object exists and is callable
  bool constructedOk() const {
    return true;
  }
    
};
}//end of ns
#endif

