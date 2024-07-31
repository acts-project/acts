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

template <typename T>
SimSpacePoint convertSpacePoint(traccc::spacepoint& spacePoint, T& measurementConv){
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
        [](auto& m) {return m.sourceLink();},
        measurementConv.valueToValue(spacePoint.meas)
    );

    boost::container::static_vector<Acts::SourceLink, 2> sourceLinks = {std::move(sourceLink)};

    return SimSpacePoint(
        globalPos,
        t,
        varRho,
        varZ,
        varT,
        std::move(sourceLinks)
    );
}

template <typename T, typename allocator_t, typename output_container_t>
auto convertSpacePoints(std::vector<traccc::spacepoint, allocator_t>& spacePoints, T& measurementConv, output_container_t& outputContainer){
    auto fn = [&measurementConv](traccc::spacepoint& spacePoint){
        return convertSpacePoint(spacePoint, measurementConv);
    };
    return Util::convert<SimSpacePoint>(spacePoints, fn, outputContainer);
}
/*
template <typename T>
traccc::spacepoint convertSpacePoint(SimSpacePoint spacePoint,T& measurementConv){
    using Scalar = typename traccc::point3::value_type;
    auto idx = spacePoint.sourceLinks()[0].get<ActsExamples::IndexSourceLink>().index();
    return traccc::spacepoint{
        traccc::point3{
            static_cast<Scalar>(spacePoint.x()), 
            static_cast<Scalar>(spacePoint.y()), 
            static_cast<Scalar>(spacePoint.z())
            },
        measurementConv.indexToValue(idx)
    };
}

template <typename T>
auto convertSpacePoints(SimSpacePointContainer& spacePoints, T& measurementConv){
    auto fn = [&measurementConv](SimSpacePoint& spacePoint){
        return convertSpacePoint(spacePoint, measurementConv);
    };
    return Util::convert<traccc::spacepoint>(spacePoints, fn);
}*/

}  // namespace Acts::TracccPlugin
