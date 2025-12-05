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
#include "Acts/Seeding/detail/CompSpacePointAuxiliaries.hpp"
#include "Acts/Utilities/Delegate.hpp"

namespace Acts::Experimental {

template <typename SeedFiller_t, typename UnCalibSp_t, typename CalibCont_t>
concept SeedFiller =
    CompositeSpacePoint<UnCalibSp_t> &&
    CompositeSpacePointContainer<CalibCont_t> &&
    requires(const SeedFiller_t& filler, const CalibrationContext& cctx,
             const Vector3& pos, const Vector3& dir, const double t0,
             const UnCalibSp_t& testSp, CalibCont_t& seedContainer) {
      { filler.goodCandidate(cctx, testSp) } -> std::same_as<bool>;
      {
        filler.candidatePull(cctx, pos, dir, t0, testSp)
      } -> std::same_as<double>;
      {
        filler.newContainer(cctx)
      } -> std::same_as<std::unique_ptr<CalibCont_t>>;
      {
        filler.append(cctx, pos, dir, t0, testSp, seedContainer)
      } -> std::same_as<void>;
    };

class CompositeSpacePointLineSeeder {
 public:
  /// @brief Configuration of the cuts to sort out generated
  ///        seeds with poor quality.
  struct Config {
    /// @brief Cut on the theta angle
    std::array<double, 2> thetaRange{0, 180. * UnitConstants::degree};
    /// @brief Cut on the intercept range
    std::array<double, 2> interceptRange{-20. * UnitConstants::m,
                                         20. * UnitConstants::m};

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
  };
  /// @brief Class constructor
  /// @param cfg Reference to the seeder configuration object
  /// @param logger Logger object used for debug print out
  explicit CompositeSpacePointLineSeeder(
      const Config& cfg,
      std::unique_ptr<const Logger> logger = getDefaultLogger(
          "CompositeSpacePointLineSeeder", Logging::Level::INFO));

  /// @brief Use the assignment of the parameter indices from the CompSpacePointAuxiliaries
  using ParIdx = detail::CompSpacePointAuxiliaries::FitParIndex;
  /// @brief Use the vector from the CompSpacePointAuxiliaires
  using Vector = detail::CompSpacePointAuxiliaries::Vector;
  /// @brief Vector containing the 5 straight segment line parameters
  using SeedParam_t = std::array<double, toUnderlying(ParIdx::nPars)>;

  /// @brief Helper struct describing the line parameters that are
  ///        tangential to a pair of straw measurements
  struct TwoCircleTangentPars {
    /// @brief Estimated angle
    double theta{0.};
    /// @brief Estimated intercept
    double y0{0.};
    /// @brief Uncertainty on the angle
    double dTheta{0.};
    /// @brief Uncertainty on the intercept
    double dY0{0.};
    /// @brief Definition of the print operator
    /// @param ostr: Mutable reference to the stream to print to
    /// @param pars: The parameters to be printed
    friend std::ostream& operator<<(std::ostream& ostr,
                                    const TwoCircleTangentPars& pars) {
      pars.print(ostr);
      return ostr;
    }

   private:
    /// @brief Actual implementation of the printing
    /// @param ostr: Mutable reference to the stream to print to
    void print(std::ostream& ostr) const;
  };

 private:
  /// @brief Cache object of a constructed & valid seed solution.
  ///        It basically consists out of the generated parameters.
  ///        the straw hits contributing to the seed & the left/right
  ///        ambiguities given the parameters of the solutions.
  ///        To avoid the copy of the memory, the hits are encoded as a
  ///        pair of indices representing the straw layer & hit number.
  template <CompositeSpacePointSorter Splitter_t>
  struct SeedSolution : public TwoCircleTangentPars {
    /// @brief Constructor taking the constructed tangential parameters &
    ///        the pointer to the splitter to associate the hits to the seed
    /// @param pars: Theta & intercept describing the tangential line
    /// @param layerSorter: Pointer to the sorter object carrying a sorted
    ///                     collection of hits that are split per logical layer
    explicit SeedSolution(const TwoCircleTangentPars& pars,
                          const Splitter_t* layerSorter)
        : TwoCircleTangentPars{pars}, spSplitter{layerSorter} {}

    /// @brief Helper function to calculate the straw signs
    std::vector<int> leftRightSigns(const TwoCircleTangentPars& pars) const;

    /// @brief Vector of radial signs of the valid hits
    std::vector<int> solutionSigns{};
    /// @brief Seed chi2
    double chi2{0.};
    /// @brief Number of straw hits on the seed
    std::size_t nStrawHits{0};
    /// @brief Definition of the print operator
    /// @param ostr: Mutable reference to the stream to print to
    /// @param sol: The seed solution to be printed
    friend std::ostream& operator<<(std::ostream& ostr,
                                    const SeedSolution& sol) {
      sol.print(ostr);
      return ostr;
    }

   private:
    /// @brief Pointer to the space point per layer splitter to gain access to the
    ///        container
    const Splitter_t* m_spSplitter{nullptr};
    /// @brief Set of hits collected onto the seed. For each element
    ///        the first index represents the layer &
    ///        the second one the particular hit in that layer
    std::vector<std::pair<std::size_t, std::size_t>> m_seedHits{};

    /// @brief Prints the seed solution to the screen
    void print(std::ostream& ostr) const;
  };

  template <CompositeSpacePoint Sp_t>
  using Selector_t = Delegate<bool(const Sp_t&)>;
  using AbortSelector_t = Delegate<bool(std::size_t layerIdx)>;
  template <CompositeSpacePointContainer Cont_t>
  using SpacePoint_t = RemovePointer_t<typename Cont_t::value_type>;

 public:
  template <CompositeSpacePointContainer Cont_t,
            CompositeSpacePointSorter<Cont_t> Splitter_t,
            CompositeSpacePointContainer CalibCont_t,
            CompositeSpacePointCalibrator<Cont_t, CalibCont_t> Calibrator_t>

  struct SeedOptions {
    /// @brief Splitter holding the straw and strip hits
    std::unique_ptr<Splitter_t> splitter{};
    using Sp_t = SpacePoint_t<Cont_t>;
    /// @brief Good hit selector
    Selector_t<Sp_t> selector{};
    AbortSelector_t abortSelector{};
    // @brief Experiment specific calibration context
    CalibrationContext calibContext{};
    // @brief radius of the straw tubes used to reject hits outside the tube
    double strawRadius{15. * UnitConstants::mm};

    /// @brief Try at the first time the external seed parameters as candidate
    bool startWithPattern{false};

    /// @brief Estimated parameters from pattern
    SeedParam_t patternParams{};

    std::size_t nGenSeeds{0};
    std::size_t nStrawCut{0};

    /// @brief Stringstream output operator
    friend std::ostream& operator<<(std::ostream& ostr,
                                    const SeedOptions& opts) {
      return opts.print(ostr);
    }
    /// @brief Prints the seed solution to the screen
    std::ostream& print(std::ostream& ostr) const;

   private:
    std::vector<SeedSolution<Cont_t>> seenSolutions{};
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
  static TangentAmbi encodeAmbiguity(const int signTop, const int signBottom);
  /// @brief Construct the line that is tangential to a pair of two straw circle measurements
  /// @param topHit: First straw hit
  /// @param bottomHit: Second straw hit
  /// @param ambi: Left right ambiguity of the bottom & top hit
  template <CompositeSpacePoint Sp_t>
  static TwoCirclePars constructTangentLine(const Sp_t& topHit,
                                            const Sp_t& bottomHit,
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

  /// @brief Creates the direction vector from the reference hit used to
  ///        construct the tangent seed and the result on theta
  /// @param refHit: Reference hit to define the local axes (Bottom hit)
  /// @param tanAngle: Theta value from the TwoCircleTangentPars
  template <CompositeSpacePoint Spt_t>
  static Vector makeDirection(const Spt_t& refHit, const double tanAngle);

 private:
  /// @brief Reference to the logger object
  const Logger& logger() const { return *m_logger; }

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

  /// @brief  Construct the final seed parameters by combining the initial
  ///         pattern parameters with the parameter from two circle tangent
  /// @param parTheta: Theta angle of the two circle tangent solution
  /// @param parY0: Intercept of the two circle tangent solution
  /// @param patternParams: Parameter estimate from the hit pattern
  SeedParam_t constructLine(const double parTheta, const double parY0,
                            const SeedParam_t& patternParams) const;
  /// @brief Check whether the generated seed parameters are within the ranges defined by the used
  /// @param tangentPars: Reference to the seed parameters to check
  bool isValidLine(const TwoCircleTangentPars& tangentPars) const;
  /// @brief Configuration object
  Config m_cfg{};
  /// @brief Logger instance
  std::unique_ptr<const Logger> m_logger{};
  /// @brief Array encoding the four possible left right solutions.
  ///        The first (second) index encodes the ambiguity of the
  ///        bottom (top) straw measurement. The order of the index
  ///        pairs is coupled to the TangentAmbi index
  static constexpr std::array<std::array<int, 2>, 4> s_signCombo{
      std::array{1, 1}, std::array{1, -1}, std::array{-1, 1},
      std::array{-1, -1}};
};
}  // namespace Acts::Experimental
#include "Acts/Seeding/CompositeSpacePointLineSeeder.ipp"
