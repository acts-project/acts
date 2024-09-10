// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Material/ISurfaceMaterial.hpp"
#include "Acts/Material/MaterialSlab.hpp"
#include "Acts/Plugins/Detray/DetrayConversionUtils.hpp"
#include "Acts/Utilities/BinningData.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Navigation/InternalNavigation.hpp"
#include "Acts/Utilities/TypeList.hpp"
#include "Acts/Utilities/IAxis.hpp"



#include <detray/io/frontend/payloads.hpp>

#include "detray/builders/detector_builder.hpp"
#include "detray/core/detector.hpp"
#include "detray/definitions/geometry.hpp"
#include "detray/io/common/geometry_reader.hpp"
#include "detray/io/common/surface_grid_reader.hpp"
#include "detray/io/frontend/detector_writer.hpp"
#include "detray/io/frontend/payloads.hpp"
#include "detray/utils/consistency_checker.hpp"



namespace Acts {


class Surface;
class SurfaceBounds;
class IAxis;

namespace Experimental {
class Detector;
}

namespace DetraySurfaceGridsConverter {
    //aadd code here

    //convertAxis
    detray::io::axis_payload convertAxis(const IAxis& ia);

    //convertGrid
    template <typename grid_type>
    detray::io::grid_payload<std::size_t, detray::io::accel_id> convertGrid(
        const grid_type& grid, bool swapAxis = false);

    //convertImpl -> probs avoidable
    template <typename index_grid>
    detray::io::grid_payload<std::size_t, detray::io::accel_id> convertImpl(
        const index_grid& indexGrid);

    template <typename instance_type>
    std::optional<detray::io::grid_payload<std::size_t, detray::io::accel_id>> convert(
                const Experimental::InternalNavigationDelegate& delegate,
                [[maybe_unused]] const instance_type& refInstance);

    template <typename... Args>
    std::vector<detray::io::grid_payload<std::size_t, detray::io::accel_id>> unrollConvert(
        const Acts::Experimental::InternalNavigationDelegate& delegate,
        Acts::TypeList<Args...>);

    detray::io::detector_grids_payload<std::size_t, detray::io::accel_id> convertSurfaceGrids(
        const Acts::Experimental::Detector& detector);

}  // namespace DetrayMaterialConverter

}  // namespace Acts
