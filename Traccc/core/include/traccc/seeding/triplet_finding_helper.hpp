/** TRACCC library, part of the ACTS project (R&D line)
 *
 * (c) 2021-2025 CERN for the benefit of the ACTS project
 *
 * Mozilla Public License Version 2.0
 */

#pragma once

// Local include(s).
#include "traccc/definitions/math.hpp"
#include "traccc/edm/spacepoint_collection.hpp"
#include "traccc/seeding/detail/doublet.hpp"
#include "traccc/seeding/detail/lin_circle.hpp"
#include "traccc/seeding/detail/triplet.hpp"

namespace traccc {

// helper function used for both cpu and gpu
struct triplet_finding_helper {
    /// Check if two doublets with common middle spacepoint can form a triplet
    ///
    /// @param spM is middle spacepoint
    /// @param lb is transformed coordinate of middle-bottom doublet
    /// @param lt is transformed coordinate of middle-top doublet
    /// @param config is configuration parameter
    /// @param iSinTheta2 is the square of sin of pitch angle
    /// @param scatteringInRegion2 is the threshold for scattering angle for the
    /// lower pT cut
    /// @param curvature is curvature of triplet
    /// @param impact_parameter is impact parameter of triplet
    ///
    /// @return boolean value for compatibility
    template <typename T>
    static inline TRACCC_HOST_DEVICE bool isCompatible(
        const edm::spacepoint<T>& spM, const lin_circle& lb,
        const lin_circle& lt, const seedfinder_config& config,
        const scalar& iSinTheta2, const scalar& scatteringInRegion2,
        scalar& curvature, scalar& impact_parameter);
};

template <typename T>
bool TRACCC_HOST_DEVICE triplet_finding_helper::isCompatible(
    const edm::spacepoint<T>& spM, const lin_circle& lb, const lin_circle& lt,
    const seedfinder_config& config, const scalar& iSinTheta2,
    const scalar& scatteringInRegion2, scalar& curvature,
    scalar& impact_parameter) {

    // add errors of spB-spM and spM-spT pairs and add the correlation term
    // for errors on spM
    scalar error2 = lt.Er() + lb.Er() +
                    static_cast<scalar>(2.f) *
                        (lb.cotTheta() * lt.cotTheta() * spM.radius_variance() +
                         spM.z_variance()) *
                        lb.iDeltaR() * lt.iDeltaR();

    scalar deltaCotTheta = lb.cotTheta() - lt.cotTheta();
    scalar deltaCotTheta2 = deltaCotTheta * deltaCotTheta;
    scalar error{0.f};
    scalar dCotThetaMinusError2{0.f};

    // if the error is larger than the difference in theta, no need to
    // compare with scattering
    if (deltaCotTheta2 - error2 > 0) {
        deltaCotTheta = math::fabs(deltaCotTheta);
        // if deltaTheta larger than the scattering for the lower pT cut, skip
        error = std::sqrt(error2);
        dCotThetaMinusError2 = deltaCotTheta2 + error2 -
                               static_cast<scalar>(2.) * deltaCotTheta * error;
        // avoid taking root of scatteringInRegion
        // if left side of ">" is positive, both sides of unequality can be
        // squared
        // (scattering is always positive)
        if (dCotThetaMinusError2 > scatteringInRegion2) {
            return false;
        }
    }

    // protects against division by 0
    scalar dU = lt.U() - lb.U();
    if (dU == static_cast<scalar>(0.)) {
        return false;
    }

    // A and B are evaluated as a function of the circumference parameters
    // x_0 and y_0
    scalar A = (lt.V() - lb.V()) / dU;
    scalar S2 = static_cast<scalar>(1.) + A * A;
    scalar B = lb.V() - A * lb.U();
    scalar B2 = B * B;
    // sqrt(S2)/B = 2 * helixradius
    // calculated radius must not be smaller than minimum radius
    if (S2 < B2 * config.minHelixDiameter2) {
        return false;
    }

    // 1/helixradius: (B/sqrt(S2))*2 (we leave everything squared)
    scalar iHelixDiameter2 = B2 / S2;
    // calculate scattering for p(T) calculated from seed curvature
    scalar pT2scatter =
        static_cast<scalar>(4.) * iHelixDiameter2 * config.pT2perRadius;
    // if pT > maxPtScattering, calculate allowed scattering angle using
    // maxPtScattering instead of pt.
    scalar pT =
        config.pTPerHelixRadius * std::sqrt(S2 / B2) / static_cast<scalar>(2.);
    if (pT > config.maxPtScattering) {
        scalar pTscatter = config.highland / config.maxPtScattering;
        pT2scatter = pTscatter * pTscatter;
    }
    // convert p(T) to p scaling by sin^2(theta) AND scale by 1/sin^4(theta)
    // from rad to deltaCotTheta
    scalar p2scatter = pT2scatter * iSinTheta2;
    // if deltaTheta larger than allowed scattering for calculated pT, skip
    if ((deltaCotTheta2 - error2 > static_cast<scalar>(0.)) &&
        (dCotThetaMinusError2 >
         p2scatter * config.sigmaScattering * config.sigmaScattering)) {
        return false;
    }

    // calculate curvature
    curvature = B / std::sqrt(S2);

    // A and B allow calculation of impact params in U/V plane with linear
    // function
    // (in contrast to having to solve a quadratic function in x/y plane)
    impact_parameter = math::fabs((A - B * spM.radius()) * spM.radius());

    if (impact_parameter > config.impactMax) {
        return false;
    }

    return true;
}

}  // namespace traccc
