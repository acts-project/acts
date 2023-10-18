// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Detector/detail/Blueprint.hpp"
#include "Acts/Detector/CylindricalContainerBuilder.hpp"
#include "Acts/Detector/GeometryIdGenerator.hpp"
#include "Acts/Surfaces/CylinderSurface.hpp"
#include "Acts/Surfaces/DiscSurface.hpp"
#include "Acts/Detector/interface/IInternalStructureBuilder.hpp"


template<typename surface_type>
class SurfaceBuilder : public Acts::Experimental::IInternalStructureBuilder {
    public :
        SurfaceBuilder(const Acts::Transform3& trf, Acts::ActsScalar p0, ActsScalar p1) :
          m_surface(Surface::makeShared<surface_type>(trf, p0, p1)) 
        {}


}