// This file is part of the Acts project.
//
// Copyright (C) 2024 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/SourceLink.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Utilities/detail/Axis.hpp"
#include "Acts/Utilities/Grid.hpp"

namespace Acts {

/// @brief Interface class for binning source links
/// into a lookup table -- a grid of source links
///
/// @tparam axis_t Type of the axis
///
template <typename axis_t>
class ISourceLinkGrid {
    using slGrid = 
        Grid<
            std::vector<SourceLink>, 
            axis_t, 
            axis_t>;

    public:
        /// @brief Virtual destructor
        virtual ~ISourceLinkGrid() = default;

        /// @brief Interface function to sort 
        /// the source links into the grid
        virtual void initialize(
            const GeometryContext& gctx, 
            std::vector<SourceLink> sourceLinks) = 0;

        /// @brief Interface function to get the 
        /// grid of source links for a 
        /// given geometry identifier
        virtual slGrid getSourceLinkTable(
            const GeometryIdentifier& geoId) const = 0;
};

}  // namespace Acts
