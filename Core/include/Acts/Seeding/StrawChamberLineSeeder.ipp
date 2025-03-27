// This file is part of the ACTS project.
//
// Copyright (C) 2025 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "Acts/Utilities/StringHelpers.hpp"
#include "Acts/Utilities/UnitVectors.hpp"
#include "Acts/Utilities/MathHelpers.hpp"
#include "Acts/Utilities/Enumerate.hpp"

#include "Acts/Surfaces/detail/LineHelper.hpp"
#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/Definitions/Common.hpp"

#include <format>


namespace Acts{
    namespace detail {
        /** @brief Simple helper utility checking whether the value is in the closed interval
         *  @param cutRange: Interval in which the value needs to reside
         *  @param value: Value to check */
        constexpr bool inRange(const std::array<double, 2>& cutRange, const double value) {
            return cutRange[0] <= value && cutRange[1] >= value;
        }
        /** @brief Calculate the unsigned residual of a segment line parametrized by position + direction
         *         with a straw measurement.
         *  @param linePos: Arbitrary point on the segment line
         *  @param lineDir: Direction of the segment line
         *  @param strawMeas: Referecne to the straw tube measurement. */
        template<StationSpacePoint UncalibSp_t>
        constexpr double calcStrawResidual(const Vector3& linePos, 
                                           const Vector3& lineDir,
                                           const UncalibSp_t& strawMeas) {
            double dist = LineHelper::signedDistance(linePos, lineDir, strawMeas.localPosition(), strawMeas.sensorDirection());
            Vector2 res{linePos.x() - strawMeas.localPosition().x(), std::abs(dist) - strawMeas.driftRadius()};
            const Eigen::Map<const ActsSquareMatrix<2>> cov(strawMeas.covariance().template block<2,2>(0,0).data());
            return std::sqrt(res.dot(cov.inverse()*res));
        }
        template<StationSpacePoint UncalibSp_t>
        constexpr int calcStrawSign(const Vector3& linePos,
                                    const Vector3& lineDir,
                                    const UncalibSp_t& strawMeas) {
            return LineHelper::signedDistance(linePos, lineDir, strawMeas.localPosition(), strawMeas.sensorDirection()) > 0 ? 1 : -1;
        }
        template<StationSpacePointContainer CalibSpCont_t>
            constexpr std::vector<int> calcStrawSigns(const Vector3& linePos,
                                                      const Vector3& lineDir,
                                                      const CalibSpCont_t& measVec) {
                std::vector<int> signs{};
                signs.reserve(measVec.size());
                std::ranges::transform(measVec, std::back_inserter(signs),
                                        [&linePos, &lineDir](const CalibSpCont_t::value_type& straw ){
                                            return calcStrawSign(linePos, lineDir, *straw);
                                        });
                return signs;
            }
    }
    template <StationSpacePointContainer UnCalibCont_t,
              StationSpacePointSorter<UnCalibCont_t> Sorter_t,
              StationSpacePointContainer CalibSp_t,
              StationSpacePointCalibrator<UnCalibCont_t, CalibSp_t> Calibrator_t>
    StrawChamberLineSeeder<UnCalibCont_t, Sorter_t, CalibSp_t,Calibrator_t>::StrawChamberLineSeeder(const UnCalibCont_t& seedHits,
                                                                                                  Config&& cfg,
                                                                                                  std::unique_ptr<const Logger> logObj):
        m_hitLayers{seedHits},
        m_cfg{std::move(cfg)},
        m_logger{std::move(logObj)} {
            
        /** The StrawChamberLine needs to have at least 2 straw layers */
        if (m_hitLayers.strawHits().size() < 2) {
            ACTS_VERBOSE(__func__<<"() "<<__LINE__<<" - Not enough straw layers have been parsed");
            return;
        }
        if (std::ranges::find_if(m_hitLayers.strawHits(), [this](const UnCalibCont_t& vec){
                return vec.size() > m_cfg.busyLayerLimit;
            }) != m_hitLayers.strawHits().end()) {
            ACTS_VERBOSE(__func__<<"() "<<__LINE__<<" - Detected at least one busy layer with more than "<<m_cfg.busyLayerLimit
                <<". Will not attempt the external seed as first seed");
            m_cfg.startWithPattern = false;
        }
        // Set the start for the upper layer
        m_upperLayer = m_hitLayers.strawHits().size()-1; 

        /** Check whether the first layer is too busy */
        while (m_lowerLayer < m_upperLayer && m_hitLayers.strawHits()[m_lowerLayer].size() >m_cfg.busyLayerLimit){
            ACTS_VERBOSE(__func__<<"() "<<__LINE__<<" - Lower layer "<<m_lowerLayer<<" has too many hits. Don't use it as seeding layer");
            ++m_lowerLayer;
        }
        /** Check whether the lower layer is too busy */
        while (m_lowerLayer < m_upperLayer && m_hitLayers.strawHits()[m_upperLayer].size() > m_cfg.busyLayerLimit) {
            ACTS_VERBOSE(__func__<<"() "<<__LINE__<<" - Upper layer "<<m_upperLayer<<" has too many hits. Don't use it as seeding layer");
            --m_upperLayer;
        }
    }
       template <StationSpacePointContainer UnCalibCont_t,
              StationSpacePointSorter<UnCalibCont_t> Sorter_t,
              StationSpacePointContainer CalibSp_t,
              StationSpacePointCalibrator<UnCalibCont_t, CalibSp_t> Calibrator_t>
    std::ostream& StrawChamberLineSeeder<UnCalibCont_t, Sorter_t, CalibSp_t,Calibrator_t>::SeedSolution::print(std::ostream& ostr) const {
        ostr<<"theta: "<<theta<<" pm "<<dTheta<<", ";
        ostr<<"Y0: "<<Y0<<" pm "<<dY0<<", ";
        return ostr;
    }
       template <StationSpacePointContainer UnCalibCont_t,
              StationSpacePointSorter<UnCalibCont_t> Sorter_t,
              StationSpacePointContainer CalibSp_t,
              StationSpacePointCalibrator<UnCalibCont_t, CalibSp_t> Calibrator_t>
    std::optional<typename StrawChamberLineSeeder<UnCalibCont_t, Sorter_t, CalibSp_t,Calibrator_t>::DriftCircleSeed>
        StrawChamberLineSeeder<UnCalibCont_t, Sorter_t, CalibSp_t,Calibrator_t>::generateSeed(const CalibrationContext& ctx) {
            std::optional<DriftCircleSeed> found{std::nullopt};
            while (m_lowerLayer < m_upperLayer) {
                const UnCalibCont_t& lower = m_hitLayers.strawHits().at(m_lowerLayer);
                const UnCalibCont_t& upper = m_hitLayers.strawHits().at(m_upperLayer);
                ACTS_VERBOSE(__func__<<"() "<<__LINE__<<" - Layers with hits: "<<m_hitLayers.strawHits().size()
                                <<" -- next bottom hit: "<<m_lowerLayer<<", hit: "<<m_lowerHitIndex
                                <<" ("<<lower.size()<<"), topHit " <<m_upperLayer<<", "<<m_upperHitIndex
                                <<" ("<<upper.size()<<") - ambiguity ("<<s_signCombos[m_signComboIndex][0]<<";"<<s_signCombos[m_signComboIndex][1]<<")");
    
                found = buildSeed(ctx, upper.at(m_upperHitIndex), lower.at(m_lowerHitIndex), s_signCombos.at(m_signComboIndex));
                /// Increment for the next candidate
                moveToNextCandidate();
                /// If a candidate is built return it. Otherwise continue the process
                if (found) {
                    return found;
                }
            }
            return std::nullopt; 
    }
    template <StationSpacePointContainer UnCalibCont_t,
              StationSpacePointSorter<UnCalibCont_t> Sorter_t,
              StationSpacePointContainer CalibSp_t,
              StationSpacePointCalibrator<UnCalibCont_t, CalibSp_t> Calibrator_t>
    void StrawChamberLineSeeder<UnCalibCont_t, Sorter_t, CalibSp_t,Calibrator_t>::moveToNextCandidate() {
        const UnCalibCont_t& lower = m_hitLayers.strawHits()[m_lowerLayer];
        const UnCalibCont_t& upper = m_hitLayers.strawHits()[m_upperLayer];
        /// Vary the left-right solutions 
        if (++m_signComboIndex < s_signCombos.size()) {
            return;
        }
        m_signComboIndex = 0; 
            
        /// Move to the next hit in the lower layer
            if (++m_lowerHitIndex < lower.size()) {
            return;
        }
        m_lowerHitIndex=0;
        /// Move to the next hit in the upper layer
        if (++m_upperHitIndex < upper.size()) {
            return;
        }
        m_upperHitIndex = 0;
        /// All combinations of hits & lines in both layers are processed
        /// Switch to the next lowerLayer. But skip the busy ones according to the configuration
        while (m_lowerLayer < m_upperLayer && m_hitLayers.strawHits()[++m_lowerLayer].size() > m_cfg.busyLayerLimit) {
        } 
            
        if (m_lowerLayer < m_upperLayer) {
            return;
        }
        m_lowerLayer = 0; 
        while (m_lowerLayer < m_upperLayer && m_hitLayers.strawHits()[--m_upperLayer].size() > m_cfg.busyLayerLimit){
        }
    }
     
    template <StationSpacePointContainer UnCalibCont_t,
              StationSpacePointSorter<UnCalibCont_t> Sorter_t,
              StationSpacePointContainer CalibSp_t,
              StationSpacePointCalibrator<UnCalibCont_t, CalibSp_t> Calibrator_t>
    std::optional<typename StrawChamberLineSeeder<UnCalibCont_t, Sorter_t, CalibSp_t,Calibrator_t>::DriftCircleSeed>
        StrawChamberLineSeeder<UnCalibCont_t, Sorter_t, CalibSp_t,Calibrator_t>::buildSeed(const CalibrationContext& ctx,
                                                                                           const UncalibSp_t& topHit, 
                                                                                           const UncalibSp_t& bottomHit, 
                                                                                           const SignCombo_t& signs) {
        
        // Fetch the signs
        const auto&[signTop, signBot] = signs;
        /// Calculate the relative radius
        double R = signBot *bottomHit->driftRadius() - signTop * topHit->driftRadius(); 
        const Vector3& bottomPos{bottomHit->localPosition()};
        const Vector3& topPos{topHit->localPosition()};
        /// Calculate the distance between the two and their relative angle
        const Vector3 D = topPos - bottomPos;
        const double thetaTubes = std::atan2(D.y(), D.z()); 
        const double distTubes =  fastHypot(D.y(), D.z());
        ACTS_VERBOSE(__func__<<"() "<<__LINE__<<" - Try to build new 2 circle seed from bottom Hit: "<<toString(bottomPos)<<", r: "<<bottomHit->driftRadius()
                    <<", top hit: "<<toString(topPos)<<", r: "<<topHit->driftRadius()
                    <<" --> tube distance: "<<toString(D)<<", mag: "<<distTubes<<", theta: "<<thetaTubes);
        
        /// Calculate the seed theta.
        double theta{thetaTubes - std::asin(std::clamp(R / distTubes, -1., 1.))};
                        
        Vector3 seedDir = makeDirectionFromPhiTheta(90.*UnitConstants::degree, theta);
        double y0 = bottomPos.y() * seedDir.z() - bottomPos.z() * seedDir.y() + signBot * bottomHit->driftRadius();
        
        DriftCircleSeed candidateSeed{};

        double combDriftUncert{std::sqrt(bottomHit->covariance()(eY, eY) + 
                                         topHit->covariance()(eY, eY))};


        candidateSeed.parameters[eBoundLoc0] = y0 / seedDir.z();
        candidateSeed.parameters[eBoundTheta] = theta;
        
        /// Check that the initial estimate of the seed is in range
        if (!detail::inRange(m_cfg.thetaRange, candidateSeed.parameters[eBoundTheta]) ||
            !detail::inRange(m_cfg.interceptRange, candidateSeed.parameters[eBoundLoc0])) {
            ACTS_VERBOSE(__func__<<"() "<<__LINE__<<" Seed parameters are out of range");
            return std::nullopt;
        }
        
        const Vector3 seedPos = candidateSeed.parameters[eBoundLoc0]  * Vector3::UnitY();

        assert(std::abs(topPos.y()*seedDir.z() - topPos.z() * seedDir.y() + signTop*topHit->driftRadius() - Y0) < std::numeric_limits<float>::epsilon() );
        ACTS_VERBOSE(__func__<<"() "<<__LINE__<<": Candidate seed theta: "<<theta<<", tanTheta: "<<(seedDir.y() / seedDir.z())
                    <<", y0: "<<candidateSeed.parameters[eBoundLoc0]);

        SeedSolution solCandidate{};
        solCandidate.Y0 = candidateSeed.parameters[eBoundLoc0];
        solCandidate.theta = theta;
        /// d/dx asin(x) = 1 / sqrt(1- x*x)
        const double denomSquare =  1. - std::pow(R / distTubes, 2); 
        if (denomSquare < std::numeric_limits<double>::epsilon()){
            ACTS_VERBOSE("Invalid seed, rejecting"); 
            return std::nullopt; 
        }
        solCandidate.dTheta =  combDriftUncert / std::sqrt(denomSquare) / distTubes;
        solCandidate.dY0 =  fastHypot(-bottomPos.y()*seedDir.y() + bottomPos.z()*seedDir.z(), 1.) * solCandidate.dTheta;
        ACTS_VERBOSE(__func__<<"() "<<__LINE__<<": Test new "<<solCandidate<<". "<<m_seenSolutions.size());
        if (std::ranges::find_if(m_seenSolutions,
            [&solCandidate, this] (const SeedSolution& seen) {
                const double deltaY = std::abs(seen.Y0 - solCandidate.Y0);
                const double limitY = fastHypot(seen.dY0, solCandidate.dY0);
                const double dTheta = std::abs(seen.theta - solCandidate.theta);
                const double limitTh = fastHypot(seen.dTheta, solCandidate.dTheta);
                ACTS_VERBOSE(__func__<<"() "<<__LINE__<<": "<<seen
                        <<std::format(" delta Y: {:.2f} {:} {:.2f}", deltaY, deltaY < limitY ? '<' : '>', limitY)
                        <<std::format(" delta theta: {:.2f} {:} {:.2f}", dTheta, dTheta < limitTh ? '<' : '>', limitTh) );
                    return deltaY < limitY && dTheta < limitTh;;
            }) != m_seenSolutions.end()){
            ACTS_VERBOSE(__func__<<"() "<<__LINE__<<": Reject due to similarity");
            return std::nullopt;
        }
        /** Collect all hits close to the seed line */
        for (const auto& [layerNr,  hitsInLayer] : enumerate(m_hitLayers.strawHits())) {
            ACTS_VERBOSE( __func__<<"() "<<__LINE__<<": "<<hitsInLayer.size()<<" hits in layer "<<(layerNr +1));
            bool hadGoodHit{false};
            for (const UncalibSp_t& testMe : hitsInLayer) {
                const double pull = detail::calcStrawResidual(seedPos, seedDir, *testMe);            
                ACTS_VERBOSE(__func__<<"() "<<__LINE__<<": Test hit: "<<toString(testMe->localPosition())
                            <<"radius: "<<testMe->driftRadius()<<", pull: "<<pull);
                if (pull < m_cfg.hitPullCut) {
                    hadGoodHit = true;
                    solCandidate.seedHits.emplace_back(testMe);
                    ++candidateSeed.nStrawHits;
                }/// what ever comes after is not matching onto the segment 
                else if (hadGoodHit) {
                    break;
                } 
            }
            /** Reject seeds with too little straw hit association */
            const unsigned hitCut = std::max(1.*m_cfg.nStrawHitCut, m_cfg.nStrawLayHitCut * m_hitLayers.strawHits().size()); 

            if (1.*candidateSeed.nStrawHits < hitCut) {
                ACTS_VERBOSE(__func__<<"() "<<__LINE__<<": Too few hits associated "<<candidateSeed.nStrawHits
                          <<", expect at least "<<hitCut<<" hits.");
                return std::nullopt;
            }
             /* Calculate the left-right signs of the used hits */
            if (m_cfg.overlapCorridor) {
                solCandidate.solutionSigns = detail::calcStrawSigns(seedPos, seedDir, solCandidate.seedHits);
                ACTS_VERBOSE(__func__<<"() "<<__LINE__<<": Circle solutions for seed - "<<solCandidate);
                /** Last check whether another seed with the same left-right combination hasn't already been found */
                for (const SeedSolution& accepted : m_seenSolutions) {
                    unsigned int nOverlap{0};
                    std::vector<int> corridor = detail::calcStrawSigns(seedPos, seedDir, accepted.seedHits);
                    for (unsigned int l =0; l < accepted.seedHits.size(); ++l) {
                        nOverlap += (corridor[l] == accepted.solutionSigns[l]);
                    }
                    /// The seed basically generates a new line that's in he same left-right corridor compared 
                    /// to a previously found solution. There's no need to return that seed again
                    if (nOverlap == corridor.size() && accepted.seedHits.size() >= solCandidate.seedHits.size()) {
                        ACTS_VERBOSE(__func__<<"() "<<__LINE__<<": Same set of hits collected within the same corridor.");
                        return std::nullopt;
                    }
                }
            }
            /** If we found a long straw hit seed, then ensure that all
               *  subsequent seeds have at least the same amount of straw hits hits. */
            if (m_cfg.tightenHitCut) {
                m_cfg.nStrawHitCut = std::max(m_cfg.nStrawHitCut, candidateSeed.nStrawHits);
            }
            ++m_nGenSeeds;
            
            if (m_cfg.fastSeedFit) {
                if (!m_cfg.fastSegFitWithT0) {
                    fitDriftCircles(candidateSeed);
                } else {
                    fitDriftCirclesWithT0(ctx, candidateSeed);
                }
            }
        }
        return candidateSeed;
    }
    
    template <StationSpacePointContainer UnCalibCont_t,
    StationSpacePointSorter<UnCalibCont_t> Sorter_t,
    StationSpacePointContainer CalibSp_t,
    StationSpacePointCalibrator<UnCalibCont_t, CalibSp_t> Calibrator_t>
    typename StrawChamberLineSeeder<UnCalibCont_t, Sorter_t, CalibSp_t,Calibrator_t>::SeedFitAuxilliaries
    
    StrawChamberLineSeeder<UnCalibCont_t, Sorter_t, CalibSp_t,Calibrator_t>::estimateAuxillaries(const DriftCircleSeed& seed) const {

        SeedFitAuxilliaries aux{};
        /// Seed direction vector
        const Vector3 seedDir= makeDirectionFromPhiTheta(90.*UnitConstants::degree,
                                                         seed.parameters[eBoundTheta]);
        /// y0Prime = y0 * cos(theta)
        const double y0 = seed.parameters[eBoundLoc0] * seedDir.z();

        aux.invCovs.reserve(seed.measurements.size());
        aux.driftSigns.reserve(seed.measurements.size());
        double norm{0.};
        /// Calculate the centre of gravity of the segment seed which is 
        /// the sum over the measurement's positions weighted by their inverse uncertainty
        for (const auto& hit : seed.measurements) {
            const double invCov = 1./ hit->covariance()(eY, eY);
            const Vector3& pos{hit->localPosition()};
            
            const int sign = y0  - pos.y() * seedDir.z() + pos.z()* seedDir.y() > 0 ? 1 : -1;

            aux.centerOfGrav+= invCov * pos;
            aux.invCovs.push_back(invCov);
            aux.driftSigns.push_back(sign);

            norm += invCov;
        }
   
        aux.covNorm = 1./ norm;
        aux.centerOfGrav *= aux.covNorm;
        /// Calculate the fit constants
        for (const auto [covIdx, hit] : Acts::enumerate(seed.measurements)) {
            const double& invCov = aux.invCovs[covIdx];
            const int& sign = aux.driftSigns[covIdx];
            const Vector3 pos = hit->localPosition() - aux.centerOfGrav;
            const double signedCov = invCov * sign;
            aux.T_zzyy += invCov * (std::pow(pos.z(), 2) - std::pow(pos.y(), 2));
            aux.T_yz   += invCov * pos.y()*pos.z();
            aux.T_rz   += signedCov * pos.z() * hit->driftRadius();
            aux.T_ry   += signedCov * pos.y() * hit->driftRadius();
            aux.fitY0  += signedCov * aux.covNorm * hit->driftRadius();
        }
        ACTS_VERBOSE(__func__<<"() "<<__LINE__<<": Estimated T_zzyy: "<<aux.T_zzyy
                 <<", T_yz: "<<aux.T_yz<<", T_rz: "<<aux.T_rz
                 <<", T_ry: "<<aux.T_ry<<", centre "<<toString(aux.centerOfGrav)<<", y0: "<<aux.fitY0
                 <<", norm: "<<aux.covNorm<<"/"<<norm);
        return aux;
    }
    
       template <StationSpacePointContainer UnCalibCont_t,
              StationSpacePointSorter<UnCalibCont_t> Sorter_t,
              StationSpacePointContainer CalibSp_t,
              StationSpacePointCalibrator<UnCalibCont_t, CalibSp_t> Calibrator_t>
    void StrawChamberLineSeeder<UnCalibCont_t, Sorter_t, CalibSp_t,Calibrator_t>::fitDriftCircles(DriftCircleSeed& inSeed) const {

        const SeedFitAuxilliaries auxVars = estimateAuxillaries(inSeed);
     
        double theta = inSeed.parameters[eBoundTheta];
        /// Now it's time to use the guestimate
        const double thetaGuess = std::atan2( 2.*(auxVars.T_yz - auxVars.T_rz), auxVars.T_zzyy) / 2.;

        ACTS_VERBOSE(__func__<<"() "<<__LINE__<<": Start fast fit seed: "<<theta
                    <<", guess: "<<thetaGuess<<", y0: "<<inSeed.parameters[eBoundLoc0]
                    <<", fitY0: "<<auxVars.fitY0<<", centre: "<<toString(auxVars.centerOfGrav));
        
        /*** Dummy helper struct to conviniently wrap cos & sin */
        struct sincos{
            sincos(const double alpha):
                cs{std::cos(alpha)},
                sn{std::sin(alpha)}{}
            double cs{0.};
            double sn{0.};
        };
        //// 
        theta = thetaGuess;
        sincos thetaCS{theta};
        bool converged{false};
        while (!converged && inSeed.nIter++ <= m_cfg.nMaxIter) {
            const sincos twoTheta{2.*theta};
            const double thetaPrime = 0.5*auxVars.T_zzyy *twoTheta.sn - auxVars.T_yz * twoTheta.cs 
                                    - auxVars.T_rz * thetaCS.cs - auxVars.T_ry * thetaCS.sn;
            if (std::abs(thetaPrime) < m_cfg.precCutOff){
                converged = true;
                break;
            }

            const double thetaTwoPrime =  auxVars.T_zzyy * twoTheta.cs + 2.* auxVars.T_yz * twoTheta.sn 
                                       + auxVars.T_rz * thetaCS.sn - auxVars.T_ry * thetaCS.cs;
            const double update = thetaPrime / thetaTwoPrime;
            ACTS_VERBOSE(__func__<<"() "<<__LINE__<<": Fit iteration # "<<inSeed.nIter
                        <<" -- theta: "<<theta<<", thetaPrime: "<<thetaPrime
                        <<", thetaTwoPrime: "<<thetaTwoPrime<<" -- "<<std::format("{:.8f}", update)
                        <<" --> next theta "<<(theta - thetaPrime / thetaTwoPrime));

            if (std::abs(update) < m_cfg.precCutOff) {
                converged = true;
                break;
            }
            theta -= update;
            thetaCS = sincos{theta};
        }
        if (!converged) {
           return;
        }
        double fitY0 = (auxVars.centerOfGrav.y() *thetaCS.cs - auxVars.centerOfGrav.z() * thetaCS.sn + auxVars.fitY0) / thetaCS.cs;
        ACTS_VERBOSE(__func__<<"() "<<__LINE__<<": Drift circle fit converged within "
                   <<inSeed.nIter<<" iterations giving "<<toString(inSeed.parameters)<<", chi2: "
                   <<inSeed.chi2<<" - theta: "<<(theta / UnitConstants::degree)<<", y0: "<<fitY0);
        inSeed.parameters[eBoundTheta] = theta;
        inSeed.parameters[eBoundLoc0] = fitY0;       
    }
       template <StationSpacePointContainer UnCalibCont_t,
              StationSpacePointSorter<UnCalibCont_t> Sorter_t,
              StationSpacePointContainer CalibSp_t,
              StationSpacePointCalibrator<UnCalibCont_t, CalibSp_t> Calibrator_t>
    void StrawChamberLineSeeder<UnCalibCont_t, Sorter_t, CalibSp_t,Calibrator_t>:: fitDriftCirclesWithT0(const CalibrationContext& ctx, 
                                                                   DriftCircleSeed& candidateSeed) const {

    }
}