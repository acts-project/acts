// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.
#pragma once
//
#include "Acts/Surfaces/LineSurface.hpp"
#include "Acts/Utilities/Definitions.hpp"
#include "Acts/Utilities/VariantDataFwd.hpp"
//
//
#include <limits>

namespace Acts {
namespace Test {

  class LineSurfaceStub : public LineSurface
  {
  public:
    LineSurfaceStub() = delete;
    //
    LineSurfaceStub(std::shared_ptr<const Transform3D> htrans,
                    double                             radius,
                    double                             halfz)
      : GeometryObject(), LineSurface(htrans, radius, halfz)
    { /* nop */
    }
    //
    LineSurfaceStub(std::shared_ptr<const Transform3D> htrans,
                    std::shared_ptr<const LineBounds>  lbounds = nullptr)
      : GeometryObject(), LineSurface(htrans, lbounds)
    { /*nop */
    }
    //
    LineSurfaceStub(std::shared_ptr<const LineBounds> lbounds,
                    const DetectorElementBase&        detelement)
      : GeometryObject(), LineSurface(lbounds, detelement)
    { /* nop */
    }

    LineSurfaceStub(const variant_data& data)
      : GeometryObject(), LineSurface(data)
    { /* nop */
    }

    //
    LineSurfaceStub(const LineSurfaceStub& ls)
      : GeometryObject(), LineSurface(ls)
    { /* nop */
    }
    //
    LineSurfaceStub(const LineSurfaceStub& ls, const Transform3D& t)
      : GeometryObject(), LineSurface(ls, t)
    { /* nop */
    }
    /// pure virtual functions of baseclass implemented here
    std::shared_ptr<LineSurfaceStub>
    clone(const Transform3D* /*shift = nullptr*/) const
    {
      return nullptr;
    }

    /// Return method for the Surface type to avoid dynamic casts
    SurfaceType
    type() const final
    {
      return Surface::Straw;
    }

    /// Simply return true to show object exists and is callable
    bool
    constructedOk() const
    {
      return true;
    }

    using Surface::normal;

  private:
    Surface*
    clone_impl(const Transform3D* /*shift*/) const
    {
      return nullptr;
    }
  };
}  // end of ns
}
