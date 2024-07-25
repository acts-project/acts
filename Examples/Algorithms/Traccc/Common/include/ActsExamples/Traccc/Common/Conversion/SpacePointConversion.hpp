// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Acts include(s)
#include "Acts/Definitions/Algebra.hpp"
#include "Acts/EventData/SourceLink.hpp"

// Acts Examples include(s)
#include "ActsExamples/EventData/SimSpacePoint.hpp"
#include "ActsExamples/EventData/Measurement.hpp"
#include "ActsExamples/Traccc/Common/Util/MapUtil.hpp"

// Traccc include(s)
#include "traccc/edm/spacepoint.hpp"

// System include(s)
#include <cstdint>
#include <cstdlib>
#include <map>
#include <variant>
#include <utility>

namespace ActsExamples::Traccc::Common::Conversion {

template <typename map_t>
SimSpacePoint convertSpacePoint(traccc::spacepoint& spacePoint, const map_t& measurementMap){
    using Scalar = Acts::ActsScalar;
    const Acts::Vector3 globalPos(
        spacePoint.x(), 
        spacePoint.y(), 
        spacePoint.z()
    );
    const std::optional<Scalar> t = std::nullopt;
    const Scalar varRho = 0;
    const Scalar varZ = 0;
    const std::optional<Scalar> varT = std::nullopt;
    const Acts::SourceLink sourceLink = std::visit(
        [](auto m) {return m.sourceLink();},
        measurementMap.at(spacePoint.meas)
    );
    boost::container::static_vector<Acts::SourceLink, 2> sourceLinks = {std::move(sourceLink)};

    return SimSpacePoint(
        globalPos,
        t,
        varRho,
        varZ,
        varT,
        sourceLinks
    );
}

template <typename allocator_t, typename map_t>
auto convertSpacePoints(std::vector<traccc::spacepoint, allocator_t>& spacePoints, const map_t& measurementMap){
    auto fn = [&measurementMap](traccc::spacepoint& spacePoint){
        return convertSpacePoint(spacePoint, measurementMap);
    };
    return Util::convert<traccc::spacepoint, SimSpacePoint>(spacePoints, fn);
}

}  // namespace Acts::TracccPlugin
