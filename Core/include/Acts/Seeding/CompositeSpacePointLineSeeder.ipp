/ This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/CompositeSpacePointLineSeeder.hpp"

namespace Acts::Experimental{


        template <Acts::Experimental::CompositeSpacePoint SpacePoint_t>
        MdtSegmentSeedGenerator::SeedSolution 
            MdtSegmentSeedGenerator::estimateTangentLine(const SpacePoint_t& topHit, const SpacePoint_t& bottomHit,
                                                         const SignComboType& signs) const {

            using namespace Acts::UnitLiterals;
            const auto&[signTop, signBot] = signs;

            SeedSolution solution{};
            const Amg::Vector3D& bottomPos{bottomHit.localPosition()};
            const Amg::Vector3D& topPos{topHit.localPosition()};

            const Amg::Vector3D D = topPos - bottomPos;
            ATH_MSG_VERBOSE(__func__<<"() - "<<__LINE__<<": Bottom position: "<<Amg::toString(bottomPos)
                    <<", driftRadius: "<<(signBot*bottomHit.driftRadius())<<" - top position "<<Amg::toString(topPos)
                    <<", driftRadius: "<<(signTop * topHit.driftRadius()) <<", tube distance: "<<Amg::toString(D));

            const double thetaTubes = std::atan2(D.y(), D.z()); 
            const double distTubes =  Acts::fastHypot(D.y(), D.z());
            constexpr auto covIdx = Acts::toUnderlying(AxisDefs::etaCov);
            const double combDriftUncert{topHit.covariance()[covIdx] + bottomHit.covariance()[covIdx]};
            const double R = signBot * bottomHit.driftRadius() - signTop * topHit.driftRadius();
            solution.theta = thetaTubes - std::asin(std::clamp(R / distTubes, -1., 1.));
            /// Outside the used allowed parameter range
            if (solution.theta < m_cfg.thetaRange[0] || solution.theta > m_cfg.thetaRange[1]) {
                solution.isValid = false;
                return solution;
            }
            const double cosTheta = std::cos(solution.theta);
            const double sinTheta = std::sin(solution.theta);
            solution.y0 = bottomPos.y()*cosTheta - bottomPos.z()*sinTheta + signBot*bottomHit.driftRadius();
            ATH_MSG_VERBOSE("Solution is theta: "<< (solution.theta / 1._degree)<<", y0: "<<solution.y0);
            assert(Acts::abs(topPos.y()*cosTheta - topPos.z()*sinTheta + signTop*topHit.driftRadius() - solution.y0) < 
                    std::numeric_limits<float>::epsilon());
            solution.y0 /= cosTheta;
            if (solution.y0 < m_cfg.interceptRange[0] || solution.y0 > m_cfg.interceptRange[1]) {
                solution.isValid = false;
                return solution;
            }
            const double denomSquare =  1. - Acts::pow(R / distTubes, 2); 
            if (denomSquare < std::numeric_limits<double>::epsilon()){
                ATH_MSG_VERBOSE("Invalid seed, rejecting"); 
                solution.isValid = false;
                 return solution; 
            }
            solution.dTheta =  combDriftUncert / std::sqrt(denomSquare) / distTubes;
            solution.dY0 =  std::hypot(bottomPos.y()*sinTheta + bottomPos.z()*cosTheta, 1.) * solution.dTheta;
            return solution;
        }
}