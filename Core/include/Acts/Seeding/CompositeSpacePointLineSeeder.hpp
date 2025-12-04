// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/CompositeSpacePoint.hpp"
#include "Acts/EventData/CompositeSpacePointCalibrator.hpp"
#include "Acts/EventData/CompositeSpacePointSorter.hpp"
#include "Acts/Seeding/CompositeSpacePointLineFitter.hpp"
#include "Acts/Seeding/detail/CompSpacePointAuxiliaries.hpp"
#include "Acts/Utilities/Delegate.hpp"

namespace Acts::Experimental {

using Line_t = detail::CompSpacePointAuxiliaries::Line_t;
using ParIdx = detail::CompSpacePointAuxiliaries::FitParIndex;

namespace detail {
/// @brief Abrivation to dereference a pointer to a const object
template <typename T>
using PointerConstDref_t = std::add_const_t<RemovePointer_t<T>>&;
}  // namespace detail

class CompositeSpacePointLineSeeder {
 public:
  using Vector = detail::CompSpacePointAuxiliaries::Vector;
  using SeedParam_t = CompositeSpacePointLineFitter::ParamVec_t;
  /// @brief typedef to the underlying space point object
  template <CompositeSpacePoint SpacePoint_t>
  using Selector_t = CompositeSpacePointLineFitter::Selector_t<SpacePoint_t>;
  using AbortSelector_t = Delegate<bool(std::size_t layerIdx)>;
  template <CompositeSpacePointContainer Cont_t>
  using SpacePoint_t = RemovePointer_t<typename Cont_t::value_type>;

  /// @brief The SeedSolution bookkeeps the straw space points used.
  ///         In principle, it's nothing else than a vector of pointers.
  ///         However, if the client decides to provide a vector of unique_ptrs,
  ///         we need to wrap the objects around a std::reference_wrapper
  template <typename Trait_t>
  struct UncalibSpBook_t {
    using type = Trait_t;
  };
  template <SmartPointerConcept Trait_t>
  struct UncalibSpBook_t<Trait_t> {
    using type = std::reference_wrapper<const Trait_t>;
  };

  template <CompositeSpacePoint SpacePoint_t>
  using SpInSeedBookCont_t =
      std::vector<typename UncalibSpBook_t<SpacePoint_t>::type>;

  struct Config {
    /// @brief Cut on the theta angle
    std::array<double, 2> thetaRange{0, 180. * UnitConstants::degree};
    /// @brief Cut on the intercept range
    std::array<double, 2> interceptRange{-20. * UnitConstants::m,
                                         20. * UnitConstants::m};

    /// @brief do not apply any cuts on the seed parameters
    bool noCutsOnSeedParams{false};
    /// @brief Upper cut on the hit chi2 w.r.t. seed in order to be associated to the seed
    double hitPullCut{5.};
    /// @brief How many drift circles may be on a layer to be used for seeding
    std::size_t busyLayerLimit{2};
    /// @brief How many drift circle hits needs the seed to contain in order to be valid
    std::size_t nStrawHitCut{3};
    /// @brief Hit cut based on the fraction of collected tube layers.
    ///        The seed must pass the tighter of the two requirements.
    double nStrawLayHitCut{2. / 3.};
    /// @brief Once a seed with even more than initially required hits is found,
    ///        reject all following seeds with less hits
    bool tightenHitCut{true};
    /// @brief Check whether a new seed candidate shares the same left-right solution with already accepted ones
    ///          Reject the seed if it has the same amount of hits
    bool overlapCorridor{true};
    /// @brief Recalibrate the seed drift circles from the initial estimate
    bool recalibSeedCircles{false};
  };

  struct SeedParameters {
    /// @brief: Theta of the line
    double theta{0.};
    /// @brief Estimated intercept
    double y0{0.};
    /// @brief Uncertainty on the angle
    double dTheta{0.};
    /// @brief Uncertainty on the intercept
    double dY0{0.};
    /// @brief Line parameters
    SeedParam_t lineParams{};
    /// @brief Vector of radial signs of the valid hits
    std::vector<int> solutionSigns{};
    /// @brief Seed chi2
    double chi2{0.};
    /// @brief Number of straw hits on the seed
    std::size_t nStrawHits{0};

    /// @brief Stringstream output operator
    friend std::ostream& operator<<(std::ostream& ostr,
                                    const SeedParameters& sol) {
      return sol.print(ostr);
    }
    /// @brief Prints the seed solution to the screen
    std::ostream& print(std::ostream& ostr) const;
  };

  /// @brief Seed object returned by the seeder. The seed contains the initial parameter estimate w.r.t
  ///         to the central plane surface inside the chamber.///Note* the
  /// parameter q/p is set to zero and does not serve any purpose here.
  ///
  ///         Further the seed contains the list of associated measurements. If
  /// the seeder is run in the fast 2D fit mode, the number of iterations w.r.t.
  /// this fit procedure is appended as well.

  /// @brief Cache of all solutions seen thus far
  template <CompositeSpacePointContainer Cont_t>
  struct SeedSolution : public SeedParameters {
    explicit SeedSolution(const SeedParameters& pars) : SeedParameters(pars) {};
    SeedSolution() = default;
    /// @brief Used hits in the seed
    Cont_t seedHits{};
  };

  template <CompositeSpacePointContainer Cont_t,
            CompositeSpacePointSorter<Cont_t> Splitter_t,
            CompositeSpacePointContainer CalibCont_t,
            CompositeSpacePointCalibrator<Cont_t, CalibCont_t> Calibrator_t>
  struct SeedOptions {
    // @brief Splitter holding the straw and strip hits
    std::unique_ptr<Splitter_t> splitter{};
    using Sp_t = SpacePoint_t<Cont_t>;
    // @brief Good hit selector
    Selector_t<Sp_t> selector{};
    AbortSelector_t abortSelector{};
    // @brief Calibrator
    const Calibrator_t* calibrator{nullptr};
    // @brief Experiment specific calibration context
    const CalibrationContext* calibContext{};
    // @brief radius of the straw tubes used to reject hits outside the tube
    double strawRadius{0.};

    /// @brief  @brief Index of the upper layer under consideration for the seeding
    std::size_t upperLayer{0};
    /// @brief Index of the lower layer under consideration for the seeding
    std::size_t lowerLayer{0};
    /// @brief  Index of the hit in the lower layer under consideration for the seeding
    std::size_t lowerHitIndex{0};
    /// @brief  Index of the hit in the upper layer under consideration for the seeding
    std::size_t upperHitIndex{0};
    /// @brief  Index of the sign combination under consideration for the seeding
    std::size_t signComboIndex{0};

    /// @brief Try at the first time the external seed parameters as candidate
    bool startWithPattern{false};

    /// @brief Estimated parameters from pattern
    SeedParam_t patternParams{};

    std::vector<SeedSolution<Cont_t>> seenSolutions{};
    std::size_t nGenSeeds{0};
    std::size_t nStrawCut{0};

    /// @brief Stringstream output operator
    friend std::ostream& operator<<(std::ostream& ostr,
                                    const SeedOptions& opts) {
      return opts.print(ostr);
    }
    /// @brief Prints the seed solution to the screen
    std::ostream& print(std::ostream& ostr) const;
  };

  /// @brief Enumeration to pick one of the four tangent lines to
  ///       the straw circle pair.
  enum class TangentAmbi : std::uint8_t {
    RR = 0,  //< Both circles are on the right side
    RL = 1,  //< The top circle is on the right and the bottom circle on the
             // left < side
    LR = 2,  //< The top circle is  on the left and the bottom circle on the
             //< right side
    LL = 3,  //< Both circles are on the left side
  };

  /// @brief Converts the line tangent ambiguity into a string
  static std::string toString(const TangentAmbi ambi);
  /// @brief Translate the combination of two drift signs into the proper
  ///        tangent ambiguity enum value
  /// @param signTop: Left/right sign of the top straw tube
  /// @param signBottom: Left/right sign of the bottom straw tube
  static constexpr TangentAmbi encodeAmbiguity(const int signTop,
                                               const int signBottom);
  /// @brief Construct the line that is tangential to a pair of two straw circle measurements
  /// @param topHit: First straw hit
  /// @param bottomHit: Second straw hit
  /// @param ambi: Left right ambiguity of the bottom & top hit
  template <CompositeSpacePoint Spt_t>
  static SeedParameters constructTangentLine(const Spt_t& topHit,
                                             const Spt_t& bottomHit,
                                             const TangentAmbi ambi);

  template <CompositeSpacePointContainer Cont_t,
            CompositeSpacePointSorter<Cont_t> Splitter_t,
            CompositeSpacePointContainer CalibCont_t,
            CompositeSpacePointCalibrator<Cont_t, CalibCont_t> Calibrator_t>
  bool prepareSeedOptions(SeedOptions<Cont_t, Splitter_t, CalibCont_t,
                                      Calibrator_t>& options) const;

  template <CompositeSpacePointContainer Cont_t,
            CompositeSpacePointSorter<Cont_t> Splitter_t,
            CompositeSpacePointContainer CalibCont_t,
            CompositeSpacePointCalibrator<Cont_t, CalibCont_t> Calibrator_t>
  std::optional<SeedSolution<CalibCont_t>> nextSeed(
      SeedOptions<Cont_t, Splitter_t, CalibCont_t, Calibrator_t>& options)
      const;

  /// @brief Class constructor
  /// @param cfg Reference to the seeder configuration object
  /// @param logger Logger object used for debug print out
  explicit CompositeSpacePointLineSeeder(
      const Config& cfg,
      std::unique_ptr<const Logger> logger = getDefaultLogger(
          "CompositeSpacePointLineSeeder", Logging::Level::INFO));

  /// @brief Reference to the logger object
  const Logger& logger() const { return *m_logger; }

  /// @brief Creates the direction vector from the reference hit used to
  ///        construct the tangent seed and the result on theta
  /// @param refHit: Reference hit to define the local axes (Bottom hit)
  /// @param tanAngle: Theta value from the TwoCircleTangentPars
  template <CompositeSpacePoint SpacePoint_t>
  static Vector makeDirection(const SpacePoint_t& refHit,
                              const double tanAngle);

 private:
  static constexpr std::array<std::array<int, 2>, 4> s_signCombo{
      std::array{1, 1}, std::array{1, -1}, std::array{-1, 1},
      std::array{-1, -1}};

  template <CompositeSpacePointContainer UnCalibCont_t>
  bool moveToNextHit(const UnCalibCont_t& hitVec,
                     const Selector_t<SpacePoint_t<UnCalibCont_t>>& selector,
                     std::size_t& hitIdx) const;

  template <CompositeSpacePointContainer UnCalibCont_t>
  bool firstGoodHit(const UnCalibCont_t& hitVec,
                    const Selector_t<SpacePoint_t<UnCalibCont_t>>& selector,
                    std::size_t& hitIdx) const;

  template <CompositeSpacePointContainer Cont_t,
            CompositeSpacePointSorter<Cont_t> Splitter_t,
            CompositeSpacePointContainer CalibCont_t,
            CompositeSpacePointCalibrator<Cont_t, CalibCont_t> Calibrator_t>
  void moveToNextCandidate(SeedOptions<Cont_t, Splitter_t, CalibCont_t,
                                       Calibrator_t>& options) const;

  template <CompositeSpacePointContainer Cont_t,
            CompositeSpacePointSorter<Cont_t> Splitter_t,
            CompositeSpacePointContainer CalibCont_t,
            CompositeSpacePointCalibrator<Cont_t, CalibCont_t> Calibrator_t>
  std::optional<SeedSolution<CalibCont_t>> buildSeed(
      SeedOptions<Cont_t, Splitter_t, CalibCont_t, Calibrator_t>& options)
      const;

  SeedParam_t constructLine(const double theta, const double y0,
                            SeedParam_t patternParams) const;

  /// @brief check if the seed line is valid within the configured cuts
  bool isValidLine(SeedParameters seedSol) const;

  Config m_cfg{};
  /// @brief Logger instance
  std::unique_ptr<const Logger> m_logger{};
};
}  // namespace Acts::Experimental
#include "Acts/Seeding/CompositeSpacePointLineSeeder.ipp"
