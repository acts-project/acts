/ This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/CompositeSpacePoint.hpp"


namespace Acts::Experimental{
    class CompositeSpacePointLineSeeder{
        public:
            /// @brief Helper struct describing the line parameters that are tangential
            ///        to a pair of straw circles
            struct TwoCircleTangentPars{
                /// @brief Estimated angle
                double theta{0.};
                /// @brief Estimated intercept
                double y0{0.};
                /// @brief Uncertainty on the angle
                double dTheta{0.};
                /// @brief Uncertainty on the intercept
                double dY0{0.};
            };
            /// @brief Enumeration
            enum class LineAmbi{

            };

            /// @brief Construct the line that is tangential to a pair of two straw circle measurements
            /// @param topHit: First straw hit
            /// @param bottomHit: Second straw hit
            /// @param signTop
            template <CompositeSpacePoint SpacePoint_t>
            static TwoCircleTangentPars constructTangentLine(const SpacePoint_t& topHit, const SpacePoint_t& bottomHit,
            const int signTop, const int signBottom);
        private:
                static constexpr std::array<std::array<int, 2>, 4>{ std::array{1,1}, std::array{1,-1},
                                                                    std::array{-1,-1}, std::array{-1,-1}};

    };
}
#include "Acts/Seeding/CompositeSpacePointLineSeeder.ipp"