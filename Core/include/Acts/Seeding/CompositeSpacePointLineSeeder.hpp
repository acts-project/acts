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
#include "Acts/Seeding/detail/CompSpacePointAuxiliaries.hpp"
#include "Acts/Utilities/Delegate.hpp"

namespace Acts::Experimental {

using Line_t = detail::CompSpacePointAuxiliaries::Line_t;
using Line_t = detail::CompSpacePointAuxiliaries::Line_t;

namespace detail {
/// @brief Abrivation to dereference a pointer to a const object
template <typename T>
using PointerConstDref_t = std::add_const_t<RemovePointer_t<T>>&;
}  // namespace detail

/// @brief Define the concept of the space point measurement sorter. The sorter shall take a collection
///         of station space points and split them into straw && strip hits.
///         Hits
/// from each category are then subdivided further into the particular detector
/// layers.
///
///      A possible implementation of the CompositeSpacePointSorter needs to
///      have
/// the following attributes
///
///      using SpVec_t =  Standard container satisfiyng the
/// CompositeSpacePointContainer concept
///
///
///      const std::vector<SpVec_t>& strawHits();
///      const std::vector<SpVec_t>& stripHits();
///  Each SpVec_t contains all measurements from a particular detector layer
template <typename Splitter_t, typename SpacePointCont_t>
concept CompositeSpacePointSorter =
    CompositeSpacePointContainer<SpacePointCont_t> &&
    requires(Splitter_t sorter) {
      /// @brief Return the straw-hit space point sorted by straw layer
      {
        sorter.strawHits()
      } -> std::same_as<const std::vector<SpacePointCont_t>&>;
      /// @brief Return the strip-hit  space points sorted by detector layer
      {
        sorter.stripHits()
      } -> std::same_as<const std::vector<SpacePointCont_t>&>;
    };

// template <CompositeSpacePointContainer UnCalibCont_t,
//           CompositeSpacePointSorter<UnCalibCont_t> Splitter_t,
//           CompositeSpacePointContainer CalibCont_t,
//           CompositeSpacePointCalibrator<UnCalibCont_t, CalibCont_t>
//           Calibrator_t>
class CompositeSpacePointLineSeeder {
 public:
  using Vector = detail::CompSpacePointAuxiliaries::Vector;
  using SeedParam_t = ActsVector<6>;
  /// @brief typedef to the underlying space point object
  template <CompositeSpacePoint SpacePoint_t>
  using Selector_t = Delegate<bool(const SpacePoint_t&)>;
  using AbortSelector_t = Delegate<bool(uint layerIdx)>;
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
    std::array<double, 2> thetaRange{0, 179. * UnitConstants::degree};
    /// @brief Cut on the intercept range
    std::array<double, 2> interceptRange{-20. * UnitConstants::m,
                                         20. * UnitConstants::m};
    /// @brief Upper cut on the hit chi2 w.r.t. seed in order to be associated to the seed
    double hitPullCut{5.};
    /// @brief Try at the first time the external seed parameters as candidate
    bool startWithPattern{false};
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
    /// @brief Pointer to the space point calibrator
    // const ISpacePointCalibrator* calibrator{nullptr};
    /// @brief Toggle whether the seed is rapidly refitted
    bool fastSeedFit{true};
    /// @brief Toggle whether an initial t0 fit shall be executed
    bool fastSegFitWithT0{false};
    /// @brief Maximum number of iterations in the fast segment fit
    std::size_t nMaxIter{100};
    /// @brief Precision cut off in the fast segment fit
    double precCutOff{1.e-6};
    /// @brief Pointer to the calibrator object
    // const Calibrator_t* calibrator{nullptr};
    /// @brief Good hit selector
    // Selector_t<UncalibSp_t> selector{};
    /* data */
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
    Line_t line{};
    /// @brief Vector of radial signs of the valid hits
    std::vector<int> solutionSigns{};
    /// @brief Iterations to obtain the seed with fast fit
    std::size_t nIter{0};
    /// @brief Seed chi2
    double chi2{0.};
    /// @brief Number of straw hits on the seed
    std::size_t nStrawHits{0};
    /// @brief Estimated parameters from pattern
    Line_t::ParamVector patternParams{};

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
    SeedSolution(const SeedParameters& pars) : SeedParameters(pars) {};
    SeedSolution() {};
    /// @brief Used hits in the seed
    // SpInSeedBookCont_t<Spt_t> seedHits{};
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

    std::size_t upperLayer{0};
    std::size_t lowerLayer{0};
    std::size_t lowerHitIndex{0};
    std::size_t upperHitIndex{0};

    uint m_signComboIndex{0};

    /// @brief Try at the first time the external seed parameters as candidate
    bool startWithPattern{false};

    /// @brief Estimated parameters from pattern
    Line_t::ParamVector patternParams{};
    double t0Estimate{0.};

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
  /// @brief Translate the array of two drift signs into the proper
  ///        tangent ambiguity enum value
  /// @param signTop: Left/right sign of the top straw tube
  /// @param signBottom: Left/right sign of the bottom straw tube
  static constexpr TangentAmbi encodeAmbiguity(const std::array<int, 2> ambi);
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
          "CompositeSpacePointLineSeeder", Logging::Level::DEBUG));

  /// @brief Reference to the logger object
  const Logger& logger() const { return *m_logger; }

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

  Line_t constructLine(const double theta, const double y0,
                       Line_t::ParamVector patternPars) const;

  /// @brief check if the seed line is valid within the configured cuts
  bool isValidLine(SeedParameters seedSol) const;

  Config m_cfg{};
  /// @brief Logger instance
  std::unique_ptr<const Logger> m_logger{};
};
}  // namespace Acts::Experimental
#include "Acts/Seeding/CompositeSpacePointLineSeeder.ipp"
