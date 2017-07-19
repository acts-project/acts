// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#ifndef ACTS_TEST_LAYERSTUB
#define ACTS_TEST_LAYERSTUB 1

#include "../Surfaces/SurfaceStub.hpp"
#include "ACTS/Layers/Layer.hpp"

namespace Acts {
/// Layer derived class stub
class LayerStub : virtual public SurfaceStub, public Layer
{
public:
  /// constructor (deleted in Surface baseclass)
  LayerStub() = delete;
  /// copy constructor (deleted in Surface baseclass)
  LayerStub(const LayerStub& otherLayer) = delete;
  /// constructor with pointer to SurfaceArray (protected in Layer baseclass)
  LayerStub(std::unique_ptr<SurfaceArray>       surfaceArray,
            double                              thickness = 0,
            std::unique_ptr<ApproachDescriptor> ad        = nullptr,
            LayerType                           ltype     = passive)
    : SurfaceStub()
    , Layer(std::move(surfaceArray), thickness, std::move(ad), ltype)
  {
  }

  /// Destructor
  virtual ~LayerStub() {}

  /// Assignment is deleted in the Layer baseclass
  LayerStub&
  operator=(const LayerStub& lay)
      = delete;

  /// surfaceRepresentation is pure virtual in baseclass
  const Surface&
  surfaceRepresentation() const
  {
    return (*this);
  }

  Surface&
  surfaceRepresentation()
  {
    return (*this);
  }

  /// Other methods have implementation in baseclass
  /// templated 'onLayer()' from baseclass ?
};
}
#endif
