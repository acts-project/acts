// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Seeding/detail/GlobalPatternFinderAuxiliaries.hpp"

#include "Acts/Seeding/detail/CompSpacePointAuxiliaries.hpp"

#include "Acts/Utilities/StringHelpers.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/VectorHelpers.hpp"
#include "Acts/Surfaces/detail/PlanarHelper.hpp"
#include "Acts/Definitions/Tolerance.hpp"

#include <format>

namespace {
    using namespace Acts::UnitLiterals;
    constexpr double inDeg(double angle) {
        return angle / 1._degree;
    }
}

namespace Acts::Experimental::detail {

template <GlobPatFinderHit Hit_t>
double computePhiVariance(const Hit_t& hit,
                          const Acts::GeometryContext& gctx) {
    const auto* sp = hit.spacePoint();
    using CovIdx = Acts::Experimental::detail::CompSpacePointAuxiliaries::ResidualIdx;
    const Vector3& pos {hit.globalPosition(gctx)};
    const Vector3 phiGrad {Acts::Vector3{-pos.y(), pos.x(), 0.} / 
            Acts::hypotSquare(pos.x(), pos.y())};

    if (sp->isStraw()) {
        const double discVar {Acts::square(sp->driftRadius()) +
            sp->covariance()[Acts::toUnderlying(CovIdx::bending)]};

        return discVar / Acts::hypotSquare(pos.x(), pos.y()) +
            Acts::square(hit.globalSensorDirection(gctx).dot(phiGrad)) * 
                (sp->covariance()[Acts::toUnderlying(CovIdx::nonBending)] - discVar);   
    }
    /// Strip measuremets      
    StripMeasurementDirections stripDirs {hit.stripMeasurementDirections(gctx)};
    /** Helper method to compute the contribution of a 1D measurement to the residual variance */
    auto oneDimContribution = [&](CovIdx idx, const Vector3& measDir) -> double {
        return hit.spacePoint()->covariance()[Acts::toUnderlying(idx)] * 
            Acts::square(measDir.dot(phiGrad));
    };
    if (!sp->measuresLoc1()) {
        /// Only phi strip measurement
        assert(!stripDirs.second);
        return oneDimContribution(CovIdx::nonBending, stripDirs.first) + 
               oneDimContribution(CovIdx::bending, hit.globalSensorDirection(gctx));    
    }
    return oneDimContribution(CovIdx::bending, stripDirs.first) + 
        (stripDirs.second ? oneDimContribution(CovIdx::nonBending, *stripDirs.second) 
                          : oneDimContribution(CovIdx::nonBending, hit.globalSensorDirection(gctx)));            
}

template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
void 
PatternStateAux<Hit_t, Sector_t, Topology_t>::moveLineAnchorHit(
    const OrderedHit& refHit) 
{
    // Treat first the special case where we have only one group
    if (countGroups() < 2) {
        // If we call this method with only one group, it means that we inverted the hit search direction without
        // finding any hit in other groups beside the initial one. So the anchor is the last added hit.
        lineAnchorHit = lastInsertedHit;
        return;
    }
    // Find first the closest group to the reference group among the pattern groups
    const auto& closestGroupIt = std::ranges::min_element(hitsPerGroup, std::ranges::less{},
        [&refHit](const auto& hits){
            if (hits.empty() || 
                    Topology_t::groupIndex(*hits.front()) == Topology_t::groupIndex(*refHit)) {
                return std::numeric_limits<int>::max();
            }
            return std::abs(hits.front().globLayer - refHit.globLayer);
        });

    // Then find the closest hit in that group to the reference hit
    const auto& hits {*closestGroupIt};
    lineAnchorHit = *std::ranges::min_element(hits, std::ranges::less{},
        [&refHit](const OrderedHit& hit){
            return std::abs(hit.globLayer - refHit.globLayer); });
}

template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
void 
PatternStateAux<Hit_t, Sector_t, Topology_t>::updateLineParameters(
    const GeometryContext& gctx,
    const BeamspotInfo& beamSpot) 
{
    if (!needLineUpdate) {
        return;
    }
    Vector3 pos1 {projToPhiPlane(gctx, *lineAnchorHit)};
    Vector3 pos2 {projToPhiPlane(gctx, *lastInsertedHit)};
    Vector3 d {pos2 - pos1};
    leverArm = d.norm();
    
    /** Check whether we have to use the beamspot instead of the anchor hit to draw the line. */
    useBeamspot = leverArm < cfg->minHitDistance4Line;
    if (useBeamspot) {
        linePos = beamSpot.position;
        d = pos2 - beamSpot.position;
        leverArm = d.norm();
    } else {
        linePos = pos1;
    }
    lineDir = d / leverArm;
    needLineUpdate = false;

    ACTS_VERBOSE(__func__<<"() Updated --> linePos: "<<toString(linePos)
        <<", lineDir: "<<toString(lineDir)<<", LeverArm: "
        <<leverArm<<", Use beamspot: "<<useBeamspot);
}
template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
double 
PatternStateAux<Hit_t, Sector_t, Topology_t>::intrinsicVariance(
    const GeometryContext& gctx,
    const Hit_t& hit,
    const Vector3& contractionVector,
    bool isProjected) 
{   
    const auto* sp = hit.spacePoint();
    using CovIdx = Acts::Experimental::detail::CompSpacePointAuxiliaries::ResidualIdx;
    if(sp->isStraw()) {
        /** For straw hits, the math simplifies according to whether the hit is projected or not. */
        const double discVar {Acts::square(sp->driftRadius()) +
            sp->covariance()[Acts::toUnderlying(CovIdx::bending)]};
        if (isProjected) {
            /** If the hit is projected, the contraction vector is J^T * residualDirection,
             *  where J is the Jacobian of the projection. The phi measurement if available is
             *  cancelled by the projection. */
            return discVar * contractionVector.squaredNorm();
        }
        const double vDotRsq {Acts::square(hit.globalSensorDirection(gctx).dot(contractionVector))};
        return discVar * (contractionVector.squaredNorm() - vDotRsq) + 
               vDotRsq * sp->covariance()[Acts::toUnderlying(CovIdx::nonBending)];
    }
    /// Strip measuremets
    StripMeasurementDirections stripDirs {hit.stripMeasurementDirections(gctx)};
    /** Helper method to compute the contribution of a 1D measurement to the residual variance */
    auto oneDimContribution = [&](CovIdx idx, const Vector3& measDir) -> double {
        return hit.spacePoint()->covariance()[Acts::toUnderlying(idx)] * 
            Acts::square(measDir.dot(contractionVector));
    };
    const double bendTerm = oneDimContribution(CovIdx::bending, stripDirs.first);
    /// The contribution of the second measured direction is cancelled in case of projected 
    /// hits and the second measured direction is same as the sensor direction.
    if (isProjected && !stripDirs.second) {
        return bendTerm;
    }
    return bendTerm + (stripDirs.second 
        ? oneDimContribution(CovIdx::nonBending, *stripDirs.second)
        : oneDimContribution(CovIdx::nonBending, hit.globalSensorDirection(gctx)));
}

template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
typename PatternStateAux<Hit_t, Sector_t, Topology_t>::LineTestRes 
PatternStateAux<Hit_t, Sector_t, Topology_t>::computeLineResidual(
    const GeometryContext& gctx, 
    const OrderedHit& testHit,
    const BeamspotInfo& beamSpot) const 
{
    LineTestRes res{};

    /** We project the test hit onto the phi plane only when the test hit does 
        *  not measure phi or when we have no phi layers, otherwise we do not project
        *  so the residual include the error in the phi direction. */
    const bool projectTestHit {!testHit->spacePoint()->measuresLoc0() || 
                               nPhiLayers == 0u};
    const Vector3 testPos {projectTestHit ? projToPhiPlane(gctx,*testHit) 
                                          : testHit->globalPosition(gctx)};
    const Vector3 K {testPos - linePos};
    const double KdotD {K.dot(lineDir)};
    
    const Vector3 R {K - KdotD * lineDir}; // Residual vector
    res.residual = R.norm();
    if (res.residual < Acts::s_epsilon) {
        /** If the residual is very small, it's likely due to bad topology, reject it. */
        res.residual = std::numeric_limits<double>::max();
        return res;
    }
    const Vector3 resDir {R / res.residual}; // Residual direction
    /** Alpha represents the extrapolation distance along the pattern line */
    const double alpha {KdotD / leverArm};
    
    /** Accumulate the derivatives of the scalar residual with respect to the common phi-plane angle. */
    double phiPlaneDerivativeAcc {0.};
    /** Accumulate the derivatives of the scalar residual with respect to the common phi-plane displacement. */
    double phiDisplDerivativeAcc {0.};
    /** Accumulate the covariance contributions of the hits to the residual */
    double residualCovAcc {0.};

    /** @brief Accumulate the covariance contributions of one hit to the residual.
        *         Each hit contributes with its intrinsic covariance and, if projected,  
        *         with the uncertainty in the common phi-plane angle of the pattern.
        *   
        *  Intrinsic covariance: for a projected hit, the residual direction is  
        *  transformed with the projection jacobian J^T * dir.
        *  dir^T * J * cov * J^T * dir = (J^T * dir)^T * cov * (J^T * dir)
        *
        *  Phi plane uncertainty: for a projected hit, the derivative of the scalar  
        *  residual with respect to the common phi-plane angle for the projected hits.
        *
        *  @param hit: the hit for which to compute the covariance
        *  @param isProjected: flag indicating if the hit is projected
        *  @param pos: projected position (needed for the phi-plane derivative of projected hits)
        *  @param preFactor: factor, function of alpha, to scale the covariance contributions */
    auto covarianceTerm = [&](const Hit_t& hit,
                              const Vector3& pos,
                              const double preFactor,
                              bool isProjected) -> void { 
        if (!isProjected) {
            residualCovAcc += Acts::square(preFactor) * 
                intrinsicVariance(gctx, hit, resDir, /*isProjected=*/false);
            return;
        }
        const Vector3 sensorDir {hit.globalSensorDirection(gctx)};
        const double projFactor {sensorDir.dot(resDir) / 
                                 sensorDir.dot(bendPlaneNorm)};
        const Vector3 trfDir {resDir - projFactor * bendPlaneNorm};
        
        residualCovAcc += Acts::square(preFactor) 
            * intrinsicVariance(gctx, hit, trfDir, /*isProjected=*/true);
        phiPlaneDerivativeAcc += preFactor * projFactor * Acts::fastHypot(pos.x(), pos.y());
        phiDisplDerivativeAcc += preFactor * projFactor;
    };

    /** Compute the covariance contributions of the first line point */
    if (useBeamspot) {
        const double covS1 = Acts::square(beamSpot.length) * Acts::square(resDir.z()) + 
                             Acts::square(beamSpot.radius) * (1 - Acts::square(resDir.z()));
        residualCovAcc += Acts::square(alpha - 1) * covS1;
    } else {
        covarianceTerm(*lineAnchorHit, linePos, alpha - 1, /*isProjected=*/true);
    }

    /** Compute the covariance contributions of the second line point */
    const Vector3 pos2 {linePos + leverArm * lineDir};
    covarianceTerm(*lastInsertedHit, pos2, -alpha, /*isProjected=*/true);

    /** Compute the covariance contributions of the test hit */
    covarianceTerm(*testHit, testPos, 1., projectTestHit);

    /** Contribution to the residual given by the projection */
    const double projSigmaTerm {
        Acts::square(phiPlaneDerivativeAcc) * bendPlaneCov[Acts::toUnderlying(BendPlaneCov::ePhiPhi)] + 
        Acts::square(phiDisplDerivativeAcc) * bendPlaneCov[Acts::toUnderlying(BendPlaneCov::eSS)] +
        2 * phiPlaneDerivativeAcc * phiDisplDerivativeAcc * bendPlaneCov[Acts::toUnderlying(BendPlaneCov::ePhiS)]};
    
    res.sigma = std::sqrt(residualCovAcc + projSigmaTerm);
    res.nDof = projectTestHit ? 1u : 2u;

    ACTS_VERBOSE(__func__<<"() "<< brief(*this)<<"\nUse beamspot: "<<useBeamspot
        <<", alpha: "<<alpha<<", Residual: "<<res.residual<<" +- "<<res.sigma
        <<", linePos: "<<toString(linePos)<<", lineDir: "<<toString(lineDir)
        <<", testPos: "<<toString(testPos)<<", resDir: "<<toString(resDir)
        <<", hit pos sigma: "<<std::sqrt(residualCovAcc)
        <<", phi plane sigma: "<<std::sqrt(projSigmaTerm));
    return res;
}
template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
Vector3
PatternStateAux<Hit_t, Sector_t, Topology_t>::projToPhiPlane(
    const Acts::GeometryContext& gctx, 
    const Hit_t& hit) const 
{
    return Acts::PlanarHelper::intersectPlane(hit.globalPosition(gctx), 
                                              hit.globalSensorDirection(gctx),
                                              bendPlaneNorm, 
                                              patPhiOffset * bendPlaneNorm).position();        
}

template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
void 
PatternStateAux<Hit_t, Sector_t, Topology_t>::updatePatternPhi(
    const GeometryContext& gctx,
    const BeamspotInfo& beamSpot)
{
    /** Helper method to align the final solution with the previous one. */
    auto alignPhiAndStore = [&](double newPatphi) -> void {
        const double refPhi {sector.phi()};
        while (newPatphi - refPhi > 90._degree) {
            newPatphi -= 180._degree;
        }
        while (newPatphi - refPhi < -90._degree) {
            newPatphi += 180._degree;
        }
        patPhi = newPatphi;
        bendPlaneNorm = Vector3{-std::sin(patPhi), std::cos(patPhi), 0.};
    };
    if (nPhiLayers == 0) {
        /** If there are no phi hits, we just use the central phi of the sector,
         *  with a standard deviation based on the expanded sector size. The plane passes
         *  through the origin by construction, so s = 0 exactly. */
        alignPhiAndStore(sector.phi());
        patPhiOffset = 0.;
        bendPlaneCov[Acts::toUnderlying(BendPlaneCov::ePhiPhi)] = Acts::square(sector.sectorSize()) / 3.;
        bendPlaneCov[Acts::toUnderlying(BendPlaneCov::eSS)] = 0.;
        bendPlaneCov[Acts::toUnderlying(BendPlaneCov::ePhiS)] = 0.;
        ACTS_VERBOSE(__func__<<"() No phi hits in the pattern, set pattern phi "
            <<inDeg(patPhi)<<", offset "<<patPhiOffset<<", BendPlaneCov "<<print(bendPlaneCov));
        return;
    }
    /** Visits the phi-measuring hits without copying them into a temporary container. */
    auto forEachPhiHit = [&](std::invocable<const Hit_t&> auto&& fn) -> void {
        std::ranges::for_each(hitsPerGroup, [&](const std::vector<OrderedHit>& hits){
            std::ranges::for_each(hits, [&](const OrderedHit& hit){
                if (hit->spacePoint()->measuresLoc0()) {
                    fn(*hit);
                }
            });
        });
        std::ranges::for_each(phiOnlyHits, [&](const Hit_t& hit){
            if (hit.spacePoint()->measuresLoc0()) {
                fn(hit);
            }
        });
    };

    if (nPhiLayers == 1) {
        // When we ahve one hit, we have an exact solution: draw a line through the hit and the beamspot.
        const Hit_t* phiHit {nullptr};
        forEachPhiHit([&](const Hit_t& hit) { phiHit = &hit; });
        assert(phiHit != nullptr);

        const Vector2 bsPos {beamSpot.position.template head<2>()};
        const Vector2 hitPos {phiHit->globalPosition(gctx).template head<2>()};
        const Vector2 d {hitPos - bsPos};
        alignPhiAndStore(std::atan2(d.y(), d.x()));
        patPhiOffset = bendPlaneNorm.dot(beamSpot.position);

        // Coordinates of the hit and the beamspot along the line, e_r = (n_y, -n_x)
        const Vector2 eR {bendPlaneNorm.y(), -bendPlaneNorm.x()};
        // Signed beamspot -> hit distance along e_r (flips sign if alignPhi rotated by 180 deg)
        const double lever {d.dot(eR)};
        // Beamspot & hit coordinates along e_R
        const double rhoBS {bsPos.dot(eR)};
        const double rhoHit {rhoBS + lever};

        const double varBS {Acts::square(beamSpot.radius)};
        const double varHit {intrinsicVariance(gctx, *phiHit, bendPlaneNorm, /*isProjected=*/false)};
        const double invL2 {1. / Acts::square(lever)};
        
        bendPlaneCov[Acts::toUnderlying(BendPlaneCov::ePhiPhi)] = (varHit + varBS) * invL2;
        bendPlaneCov[Acts::toUnderlying(BendPlaneCov::ePhiS)] = - (rhoBS * varHit + rhoHit * varBS) * invL2;
        bendPlaneCov[Acts::toUnderlying(BendPlaneCov::eSS)] = 
            (Acts::square(rhoBS) * varHit + Acts::square(rhoHit) * varBS) * invL2;;
        ACTS_VERBOSE(__func__<<"() One phi hit in the pattern, set pattern phi "
            <<inDeg(patPhi)<<", offset "<<patPhiOffset<<", BendPlaneCov "<<print(bendPlaneCov));
        return;
    }

    // Use the previous phi value as the first estimate to freeze the weights of the linear regression
    const Vector3 n0 {-std::sin(patPhi), std::cos(patPhi), 0.};

    /** Positions are taken relative to the beam spot: this avoids the cancellation in the
     *  centred moments and puts the beam spot at the origin, so it only adds its weight. */
    const Vector2 bsPos {beamSpot.position.template head<2>()};
    double weightSum {0.}, Sxx {0.}, Sxy {0.}, Syy {0.};
    Vector2 centroid {Vector2::Zero()}; //Sx, Sy
    forEachPhiHit([&](const Hit_t& hit){
        const double variance {
            intrinsicVariance(gctx, hit, n0, /*isProjected=*/false)
        };
        if (variance < Acts::s_epsilon) {
            ACTS_WARNING(__func__<<"() Skipping phi hit with zero variance: "<<*hit.spacePoint());
            return;
        }
        double weight {1. / variance};
        const Vector2 d {hit.globalPosition(gctx).template head<2>() - bsPos};
        const Vector2 wd {weight * d};

        weightSum += weight;
        centroid += wd;
        Sxx += wd.x() * d.x();
        Sxy += wd.x() * d.y();
        Syy += wd.y() * d.y();
    });
    // Beam spot as an additional point measurement, sitting at the origin of the shifted frame
    weightSum += 1. / Acts::square(beamSpot.radius);

    if (weightSum < Acts::s_epsilon) {
        throw std::runtime_error(std::format(
            "Unexpected to have zero total weight when updating pattern phi with {} phi hits and beamspot", 
            nPhiLayers));
    }
    const double weightSumInv {1. / weightSum};
    centroid *= weightSumInv;

    // Weighted centered scatter matrix (translation invariant).
    const double Mxx {Sxx - weightSum * Acts::square(centroid.x())};
    const double Mxy {Sxy - weightSum * centroid.x() * centroid.y()};
    const double Myy {Syy - weightSum * Acts::square(centroid.y())};

    const double eigenvalSeparation {Acts::fastHypot(Mxx - Myy, 2. * Mxy)};
    if (eigenvalSeparation < Acts::s_epsilon) {
        throw std::runtime_error(std::format(
            "Unexpected to have zero eigenvalue separation when updating pattern phi with {} phi hits and beamspot", 
            nPhiLayers));
    }
    const double eigenValSeparationInv {1. / eigenvalSeparation};

    alignPhiAndStore(0.5 * std::atan2(2. * Mxy, Mxx - Myy));

    // Back to the original frame
    centroid += bsPos;
    patPhiOffset = bendPlaneNorm.template head<2>().dot(centroid);
    // Centroid coordinate along the line, e_r = (n_y, -n_x)
    const Vector2 eR {bendPlaneNorm.y(), -bendPlaneNorm.x()};
    const double rhoCentroid {centroid.dot(eR)};

    bendPlaneCov[Acts::toUnderlying(BendPlaneCov::ePhiPhi)] = eigenValSeparationInv;
    bendPlaneCov[Acts::toUnderlying(BendPlaneCov::ePhiS)] = - rhoCentroid * eigenValSeparationInv;
    bendPlaneCov[Acts::toUnderlying(BendPlaneCov::eSS)] = 
        Acts::square(rhoCentroid) * eigenValSeparationInv + weightSumInv;
    
    ACTS_VERBOSE(__func__<<"() Updated pattern phi "
        <<inDeg(patPhi)<<", offset "<<patPhiOffset<<", BendPlaneCov "<<print(bendPlaneCov)
        <<", with "<<static_cast<int>(nPhiLayers)<<" phi hits and beamspot");
}

template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
void 
PatternStateAux<Hit_t, Sector_t, Topology_t>::addHit(
    const GeometryContext& gctx,
    const OrderedHit& hit,
    const double residual,
    const double resSigma,
    const BeamspotInfo& beamSpot) 
{
    /** Add the new hit */
    hitsPerGroup[Topology_t::groupIndex(*hit)].push_back(hit);

    /** Update the hit counts */
    if (hit->isPrecision()) {
        nPrecisionLayers++;
    } else {
        nTriggerLayers++;
    }
    if (hit->spacePoint()->measuresLoc0()) {
        nPhiLayers++;
        updatePatternPhi(gctx, beamSpot);
    }

    /** Update the pointers to previous layer hit */
    const bool isNewSGroup {Topology_t::groupIndex(*hit) != 
                            Topology_t::groupIndex(*lastInsertedHit)};
    prevLayerHit = lastInsertedHit;
    lastInsertedHit = hit;

    meanNormResidual2 += Acts::square(residual / resSigma);
    lastResSigma = resSigma;
    lastResidual = residual;

    /** If the new compatible hit is in a different group, update the line anchor */
    if (isNewSGroup) {
        moveLineAnchorHit(hit);
    }
    needLineUpdate = true;
}

template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
void 
PatternStateAux<Hit_t, Sector_t, Topology_t>::overWriteHit(
    const GeometryContext& gctx,
    const OrderedHit& newHit,
    const double newResidual,
    const double newResSigma,
    const BeamspotInfo& beamSpot) 
{
    const auto group {Topology_t::groupIndex(*newHit)};
    if (group != Topology_t::groupIndex(*lastInsertedHit) || 
            lastInsertedHit.globLayer != newHit.globLayer) {
        throw std::runtime_error(std::format(
            "Trying to overwrite a hit in group/layer {}/{} with another one from group/layer {}/{}", 
                static_cast<int>(Topology_t::groupIndex(*lastInsertedHit)), 
                lastInsertedHit.globLayer, static_cast<int>(group), newHit.globLayer));
    }
    /* We expect to overwrite hits of the same type (precision/trigger), since we only branch when we have 
        * compatible hits in the same layer, except for sTGC hits, where we have pad and strips in the same layer */
    if (lastInsertedHit->isPrecision() != newHit->isPrecision()) {
        if (newHit->isPrecision()) {
            nPrecisionLayers++;
            nTriggerLayers--;
        } else {
            nPrecisionLayers--;
            nTriggerLayers++;
        }
    }
    /** Update the phi counts */
    bool updatePhi {false};
    if (lastInsertedHit->spacePoint()->measuresLoc0()) {
        nPhiLayers--;
        updatePhi = true;
    }
    if (newHit->spacePoint()->measuresLoc0()) {
        nPhiLayers++;
        updatePhi = true;
    }
    /** Update the residual */
    meanNormResidual2 += Acts::square(newResidual / newResSigma) - 
                         Acts::square(lastResidual / lastResSigma);
    lastResSigma = newResSigma;
    lastResidual = newResidual;

    auto& stHits {hitsPerGroup[group]};
    if (stHits.back() != lastInsertedHit) {
        std::stringstream ss {};
        ss << "Trying to overwrite a hit that is not the last inserted hit in group/layer " 
        << static_cast<int>(group) << "/" << lastInsertedHit.globLayer << "\n";
        ss << "Last inserted hit: " << *lastInsertedHit->spacePoint() << "\n";
        ss << "Last hit in group: " << *stHits.back()->spacePoint();
        throw std::runtime_error(ss.str());
    }
    stHits.pop_back();

    /** Add the new hit */
    stHits.push_back(newHit);
    lastInsertedHit = newHit;

    if (updatePhi) {
        updatePatternPhi(gctx, beamSpot);
    }
    needLineUpdate = true;
}

template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
bool 
PatternStateAux<Hit_t, Sector_t, Topology_t>::isInPattern(
    const Hit_t& hit) const 
{
    const auto& hits {hitsPerGroup[Topology_t::groupIndex(hit)]};
    return std::ranges::find_if(hits, 
        [&hit](const OrderedHit& c){ return *c == hit; }) != hits.end();                                           
}

template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
double 
PatternStateAux<Hit_t, Sector_t, Topology_t>::getMeanResidual2() const {
    if (isFinalized) {
        return meanNormResidual2;
    }
    return meanNormResidual2 / nBendingLayers();
}

template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
typename Topology_t::LayerIdx 
PatternStateAux<Hit_t, Sector_t, Topology_t>::nBendingLayers() const {
    return nPrecisionLayers + nTriggerLayers;
}

template<GlobPatFinderHit Hit_t,
         SectorType Sector_t,
         PatternTopology<Hit_t> Topology_t>
typename Topology_t::GroupIdx 
PatternStateAux<Hit_t, Sector_t, Topology_t>::countGroups() const 
{
    return std::ranges::count_if(hitsPerGroup, 
        [](const auto& hits){ return !hits.empty(); });
}

template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
typename PatternStateAux<Hit_t, Sector_t, Topology_t>::PatternPrintView 
PatternStateAux<Hit_t, Sector_t, Topology_t>::brief(const PatternStateAux& p) {
    return {p, /*detailed=*/false};
}

template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
typename PatternStateAux<Hit_t, Sector_t, Topology_t>::PatternPrintView 
PatternStateAux<Hit_t, Sector_t, Topology_t>::detailed(const PatternStateAux& p) {
    return {p, /*detailed=*/true};
}


template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
PatternStateAux<Hit_t, Sector_t, Topology_t>::PatternStateAux(
    const GeometryContext& gctx,
    const OrderedHit& seed,
    const typename Sector_t::Index_t sector,    
    const BeamspotInfo& beamSpot,      
    const Config* cfg,
    const Acts::Logger* logger)
        : cfg{cfg},
          m_logger{logger},
          lastInsertedHit{seed},
          prevLayerHit{seed},
          lineAnchorHit{seed},
          patTheta{VectorHelpers::theta(seed->globalPosition(gctx))},
          sector{sector} {
                
        /** Add the new hit */
        hitsPerGroup[Topology_t::groupIndex(*seed)].push_back(seed);

        /** Update the hit counts */
        if (seed->isPrecision()) {
            nPrecisionLayers++;
        } else {
            nTriggerLayers++;
        }
        if (seed->spacePoint()->measuresLoc0()) {
            nPhiLayers++;
        }
        updatePatternPhi(gctx, beamSpot);
        needLineUpdate = true;
    }

template<GlobPatFinderHit Hit_t, 
         SectorType Sector_t, 
         PatternTopology<Hit_t> Topology_t>
void PatternStateAux<Hit_t, Sector_t, Topology_t>::print(
    std::ostream& ostr, 
    bool detailed) const 
{
    using namespace Acts::UnitLiterals;
    auto inDeg = [](double angle) { return angle / 1._degree; };
    ostr<<"PatternStateAux Exp Sector: "<<static_cast<int>(sector.sector())
        <<", Theta: "<<inDeg(patTheta) << ", BendPlane Phi "<<inDeg(patPhi)
        <<", displ "<<patPhiOffset<<", cov: "<<print(bendPlaneCov);
    ostr<<", nPrec: "<<static_cast<int>(nPrecisionLayers)<<", nEtaNonPrec: "
        <<static_cast<int>(nTriggerLayers)<<", nPhi: "<<static_cast<int>(nPhiLayers);
    ostr<<", mean norma res sq: "<<getMeanResidual2()<<", dir: "<<toString(lineDir);
    ostr<<", Hit per group: \n";
    for (typename Topology_t::GroupIdx g = 0u; g < Topology_t::nGroups; ++g) {
        const auto& hits {hitsPerGroup[g]};
        if (hits.empty()) {
            continue;
        }

        ostr<<"  Group "<<static_cast<int>(g)<<" has "<<hits.size()<<" hits ";
        if (detailed) {
            ostr<<"\n";
            for (const auto& hit : hits) {
                ostr<<"    "<<*hit->spacePoint()<<", globLay: "<<static_cast<int>(hit.globLayer)<<"\n";
            }
        }
    }
    if (!detailed) {
        ostr<<"\n";
        ostr<<"  Last inserted hit: "<<lastInsertedHit<<"\n";
        ostr<<"  Previous layer hit: "<<prevLayerHit<<"\n";
        ostr<<"  Line anchor hit: "<<lineAnchorHit<<"\n";
    }
}

template <GlobPatFinderHit Hit_t, 
          SectorType Sector_t, 
          PatternTopology<Hit_t> Topology_t>
void
PatternStateAux<Hit_t, Sector_t, Topology_t>::OrderedHit::print(std::ostream& ostr) const {
    ostr<<*hitPtr->spacePoint()<< ", group: " << static_cast<int>(Topology_t::groupIndex(*hitPtr)) 
        <<", globLay: "<<static_cast<int>(globLayer);
}
}