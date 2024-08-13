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

// Custom hash functions
// as some are not defined by std::hash and std::equal_to

struct TracccSpacePointHash{
    std::size_t operator()(const traccc::spacepoint& s) const noexcept {
        return s.meas.measurement_id;
    }
};

struct ActsSpacePointHash{
    std::size_t operator()(const SimSpacePoint& s) const noexcept {
        assert(s.sourceLinks().size() >= 1);
        return s.sourceLinks()[0].get<ActsExamples::IndexSourceLink>().index();
    }
};

/// @brief Converts a traccc space point to an acts space point.
/// @param spacePoint the traccc space point.
/// @param measurementConv the measurement ConversionData (traccc measurement -> acts bound variant measurement).
/// @returns An acts space point.
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

/// @brief Converts a collection of traccc space points and appends the result to the given outputContainer.
/// @param spacepoints the traccc spacepoints.
/// @param meausrementConv the measurement ConversionData (traccc measurement  -> acts bound varaint measurement).
/// @param outputContainer the container to put the converted space points into (an empty container is expected).
/// @returns The spacepoint ConversionData (traccc space point -> acts space point)
template <typename T, typename allocator_t, typename output_container_t>
auto convertSpacePoints(std::vector<traccc::spacepoint, allocator_t>& spacePoints, T& measurementConv, output_container_t& outputContainer){
    auto fn = [&measurementConv](traccc::spacepoint& spacePoint){
        return convertSpacePoint(spacePoint, measurementConv);
    };
    return Util::convert<TracccSpacePointHash, std::equal_to<traccc::spacepoint>>(spacePoints, fn, outputContainer);
}

/// @brief Converts a space point to a traccc space point.
/// @param spacePoint the space point.
/// @param measurementConv the measurement ConversionData ( bound variant measurement -> traccc measurement).
/// @returns A traccc space point.
template <typename T>
traccc::spacepoint convertSpacePoint(const SimSpacePoint spacePoint, T& measurementConv){
    using Scalar = typename traccc::point3::value_type;
    auto idx = spacePoint.sourceLinks()[0].get<ActsExamples::IndexSourceLink>().index();

    traccc::point3 global{
        static_cast<Scalar>(spacePoint.x()), 
        static_cast<Scalar>(spacePoint.y()), 
        static_cast<Scalar>(spacePoint.z())
        };
    return traccc::spacepoint{
        global,
        measurementConv.indexToValue(idx)
    };
}


/// @brief Converts a collection of space points to traccc space points and appends the result to the given outputContainer.
/// @param spacepoints the traccc space points.
/// @param meausrementConv the measurement ConversionData (acts bound varaint measurement -> traccc measurement).
/// @param outputContainer the container to put the converted space points into (an empty container is expected).
/// @returns The spacepoint ConversionData (acts space point -> traccc space point)
template <typename T, typename output_container_t>
auto convertSpacePoints(const SimSpacePointContainer& spacePoints, T& measurementConv, output_container_t& outputContainer){
    auto fn = [&measurementConv](const SimSpacePoint& spacePoint){
        return convertSpacePoint(spacePoint, measurementConv);
    };
    return Util::convert<ActsSpacePointHash, std::equal_to<SimSpacePoint>>(spacePoints, fn, outputContainer);
}

}  // namespace Acts::TracccPlugin
