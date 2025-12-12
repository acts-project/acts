// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Propagator/NavigationTarget.hpp"

#include "Acts/Geometry/BoundarySurfaceT.hpp"
#include "Acts/Geometry/Layer.hpp"
#include "Acts/Geometry/Portal.hpp"
#include "Acts/Surfaces/Surface.hpp"

namespace Acts {
void NavigationTarget::print(std::ostream& ostr) const {
  std::visit(
      [&](const auto& target) {
        using Target_t = std::decay_t<decltype(target)>;
        if constexpr (std::is_same_v<Target_t, const Surface*>) {
          ostr << target->type() << " surface " << target->geometryId();
        } else if constexpr (std::is_same_v<Target_t, const Layer*>) {
          ostr << target->layerType() << " layer " << target->geometryId();
        } else if constexpr (std::is_same_v<Target_t, const BoundarySurface*>) {
          const auto& surf{target->surfaceRepresentation()};
          ostr << "Boundary " << surf.type() << " surface "
               << surf.geometryId();
        } else if constexpr (std::is_same_v<Target_t, const Portal*>) {
          const auto& surf{target->surface()};
          ostr << surf.type() << " portal " << surf.geometryId();
        }
      },
      m_target);
}
}  // namespace Acts
