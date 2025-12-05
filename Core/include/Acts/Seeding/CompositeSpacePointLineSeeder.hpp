// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/EventData/CompositeSpacePoint.hpp"
#include "Acts/Seeding/detail/CompSpacePointAuxiliaries.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"
#include "Acts/Utilities/Delegate.hpp"

namespace Acts::Experimental {
namespace detail {
/// @brief Concept for a utility class to fill the space points from
///        one container to another. The filler is templated over the
///        SpacePoint type from the uncalibrated source container over
///        the target space point container type, which may differ
///        The filler is responsible for the conversion between the two
///        space point formats and to apply a re-calibation if needed
template <typename SeedFiller_t, typename UnCalibSp_t, typename CalibCont_t>
concept CompositeSpacePointSeedFiller =
    CompositeSpacePoint<UnCalibSp_t> &&
    CompositeSpacePointContainer<CalibCont_t> &&
    requires(const SeedFiller_t& filler, const CalibrationContext& cctx,
             const Vector3& pos, const Vector3& dir, const double t0,
             const UnCalibSp_t& testSp, CalibCont_t& seedContainer) {
      /// @brief Utility function to choose the good straw space points
      ///        for seeding
      /// @param testSp: Reference to the straw-type measurement to test
      { filler.goodCandidate(testSp) } -> std::same_as<bool>;
      /// @brief Calculates the pull of the space point w.r.t. to the
      ///        candidate seed line. To improve the pull's precision
      ///        the function may call the calibrator in the backend
      /// @param cctx: Reference to the calibration context to pipe
      ///              the hook for conditions access to the caller
      /// @param pos: Position of the cancidate seed line
      /// @param dir: Direction of the candidate seed line
      /// @param t0: Offse in the time of arrival of the particle
      /// @param testSp: Reference to the straw space point to test
      {
        filler.candidatePull(cctx, pos, dir, t0, testSp)
      } -> std::same_as<double>;
      /// @brief Creates a new empty container to construct a new segment seed
      /// @param cctx: Calibration context in case that the container shall be
      ///              part of a predefined memory block
      { filler.newContainer(cctx) } -> std::same_as<CalibCont_t&&>;
      /// @brief Appends the candidate candidate space point to the segment
      ///        seed container & optionally calibrates the parameters
      /// @param cctx: Reference to the calibration context to pipe
      ///              the hook for conditions access to the caller
      /// @param pos: Position of the cancidate seed line
      /// @param dir: Direction of the candidate seed line
      /// @param t0: Offse in the time of arrival of the particle
      /// @param testSp: Reference to the straw space point to test
      /// @param seedContainer: Reference to the target container to
      ///                       which the space point is appended to
      {
        filler.append(cctx, pos, dir, t0, testSp, seedContainer)
      } -> std::same_as<void>;
    };
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
    requires(const Splitter_t& sorter) {
      /// @brief Return the straw-hit space point sorted by straw layer
      {
        sorter.strawHits()
      } -> std::same_as<const std::vector<SpacePointCont_t>&>;
      /// @brief Return the strip-hit  space points sorted by detector layer
      {
        sorter.stripHits()
      } -> std::same_as<const std::vector<SpacePointCont_t>&>;
    };

template <typename SeedAuxiliary_t, typename UnCalibCont_t,
          typename CalibCont_t>
concept CompSpacePointSeederDelegate =
    CompositeSpacePointSeedFiller<
        SeedAuxiliary_t, RemovePointer_t<typename UnCalibCont_t::value_type>,
        CalibCont_t> &&
    CompositeSpacePointSorter<SeedAuxiliary_t, UnCalibCont_t>;

}  // namespace detail

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
    /// @brief Flag indicating which solution is constructed
    TangentAmbi ambi{TangentAmbi::LL};
    /// @brief Definition of the print operator
    /// @param ostr: Mutable reference to the stream to print to
    /// @param pars: The parameters to be printed
    friend std::ostream& operator<<(std::ostream& ostr,
                                    const TwoCircleTangentPars& pars) {
      pars.print(ostr);
      return ostr;
    }

   protected:
    /// @brief Actual implementation of the printing
    /// @param ostr: Mutable reference to the stream to print to
    virtual void print(std::ostream& ostr) const;
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
  static TwoCircleTangentPars constructTangentLine(const Sp_t& topHit,
                                                   const Sp_t& bottomHit,
                                                   const TangentAmbi ambi);
  /// @brief Creates the direction vector from the reference hit used to
  ///        construct the tangent seed and the result on theta
  /// @param refHit: Reference hit to define the local axes (Bottom hit)
  /// @param tanAngle: Theta value from the TwoCircleTangentPars
  template <CompositeSpacePoint Spt_t>
  static Vector makeDirection(const Spt_t& refHit, const double tanAngle);

 private:
  /// @brief Cache object of a constructed & valid seed solution.
  ///        It basically consists out of the generated parameters.
  ///        the straw hits contributing to the seed & the left/right
  ///        ambiguities given the parameters of the solutions.
  ///        To avoid the copy of the memory, the hits are encoded as a
  ///        pair of indices representing the straw layer & hit number.
  template <CompositeSpacePointContainer UnCalibCont_t,
            detail::CompositeSpacePointSorter<UnCalibCont_t> Splitter_t>
  struct SeedSolution : public TwoCircleTangentPars {
    /// @brief Constructor taking the constructed tangential parameters &
    ///        the pointer to the splitter to associate the hits to the seed
    /// @param pars: Theta & intercept describing the tangential line
    /// @param layerSorter: Pointer to the sorter object carrying a sorted
    ///                     collection of hits that are split per logical layer
    explicit SeedSolution(const TwoCircleTangentPars& pars,
                          const Splitter_t* layerSorter)
        : TwoCircleTangentPars{pars}, m_splitter{layerSorter} {}

    /// @brief Abrivation of the underlying space point reference
    using SpacePoint_t = ConstDeRef_t<typename UnCalibCont_t::value_type>;
    /// @brief Helper function to calculate the straw signs of the seed hits
    ///        cached by this solution w.r.t. an external line
    /// @param seedPos: Reference point of the segment line
    /// @param seedDir: Direction of the segment line
    std::vector<int> leftRightAmbiguity(const Vector& seedPos,
                                        const Vector3& seedDir) const;

    /// @brief Returns the number of hits cached in the seed
    std::size_t size() const { return m_seedHits.size(); }
    /// @brief Returns the i-th seed hit
    /// @param idx: Index of the hit to to return
    SpacePoint_t getHit(const std::size_t idx) const;
    /// @brief Appends a new seed hit to the solution
    void append(const std::size_t layIdx, const std::size_t hitIdx);
    /// @brief Vector of the associate left-rignt ambiguities
    std::vector<int> solutionSigns{};

   private:
    /// @brief Pointer to the space point per layer splitter to gain access to the
    ///        input space point container
    const Splitter_t* m_splitter{nullptr};
    /// @brief Set of hits collected onto the seed. For each element
    ///        the first index represents the layer &
    ///        the second one the particular hit in that layer
    std::vector<std::pair<std::size_t, std::size_t>> m_seedHits{};
    /// @brief Prints the seed solution to the screen
    /// @param ostr: Mutable reference to the stream to print to
    void print(std::ostream& ostr) const override final;
  };

 public:
  /// @brief Struct that represents a segment seed
  template <CompositeSpacePointContainer contType_t>
  struct SegmentSeed {
    /// @brief Constructor taking the seed parameters &&
    ///        a new hit container
    /// @param _pars:
    /// @param _hits
    explicit SegmentSeed(SeedParam_t _pars, contType_t&& _hits) noexcept
        : parameters{std::move(_pars)}, hits{std::move(_hits)} {}
    /// @brief Seed line parameters
    SeedParam_t parameters;
    /// @brief Collection of hits
    contType_t hits;
  };

  template <CompositeSpacePointContainer UncalibCont_t,
            CompositeSpacePointContainer CalibCont_t,
            detail::CompSpacePointSeederDelegate<UncalibCont_t, CalibCont_t>
                Delegate_t>
  struct SeedOptions {
    /// @brief Splitter holding the straw and strip hits
    std::unique_ptr<Delegate_t> delegate{};
    /// @brief radius of the straw tubes used to reject hits outside the tube
    double strawRadius{15. * UnitConstants::mm};
    /// @brief Try at the first time the external seed parameters as candidate
    bool startWithPattern{false};
    /// @brief Estimated parameters from pattern
    SeedParam_t patternParams{};
    /// @brief Number of generated seeds
    std::size_t nStrawCut{0ul};
    /// @brief Stringstream output operator
    friend std::ostream& operator<<(std::ostream& ostr,
                                    const SeedOptions& opts) {
      opts.print(ostr);
      return ostr;
    }

    std::size_t nGenSeeds() const;

    friend CompositeSpacePointLineSeeder;

   private:
    /// @brief Prints the seed solution to the screen
    void print(std::ostream& ostr) const;
    /// @brief List of straw measurement already constructed straw measurement seeds
    std::vector<SeedSolution<UncalibCont_t, Delegate_t>> m_seenSolutions{};
    /// @brief  @brief Index of the upper layer under consideration for the seeding
    std::optional<std::size_t> m_upperLayer{std::nullopt};
    /// @brief Index of the lower layer under consideration for the seeding
    std::optional<std::size_t> m_lowerLayer{std::nullopt};
    /// @brief  Index of the hit in the lower layer under consideration for the seeding
    std::size_t m_lowerHitIndex{0ul};
    /// @brief  Index of the hit in the upper layer under consideration for the seeding
    std::size_t m_upperHitIndex{0ul};
    /// @brief  Index of the sign combination under consideration for the seeding
    std::size_t m_signComboIndex{0ul};
  };

  template <CompositeSpacePointContainer UncalibCont_t,
            CompositeSpacePointContainer CalibCont_t,
            detail::CompSpacePointSeederDelegate<UncalibCont_t, CalibCont_t>
                Delegate_t>
  std::optional<SegmentSeed<CalibCont_t>> nextSeed(
      const CalibrationContext& cctx,
      SeedOptions<UncalibCont_t, CalibCont_t, Delegate_t>& options) const;

 private:
  /// @brief Reference to the logger object
  const Logger& logger() const { return *m_logger; }
  /// @brief Abrivation of the selector delegate to skip invalid straw hits in the seed
  template <CompositeSpacePointContainer Cont_t>
  using Selector_t = Delegate<bool(ConstDeRef_t<typename Cont_t::value_type>)>;
  /// @brief
  template <CompositeSpacePointContainer Cont_t>
  using StrawLayers_t = std::vector<Cont_t>;
  /// @brief Moves to the hit index to the next good hit inside the layer.
  ///        The index is incremented until the underlying hit is accepted
  ///        by the selector or all hits in the container were tried
  /// @param hitVec: Reference to  the straw hits inside the layer
  /// @param selector: Delegate method to skip bad bad hits
  /// @param hitIdx: Mutable reference to the index that incremented
  template <CompositeSpacePointContainer UnCalibCont_t>
  bool moveToNextHit(const UnCalibCont_t& hitVec,
                     const Selector_t<UnCalibCont_t>& selector,
                     std::size_t& hitIdx) const;
  /// @brief Sets the parsed index to the first good hit inside the straw layer.
  /// @param hitVec: Reference to  the straw hits inside the layer
  /// @param selector: Delegate method to skip bad bad hits
  /// @param hitIdx: Mutable reference to the index that incremented
  template <CompositeSpacePointContainer UnCalibCont_t>
  bool firstGoodHit(const UnCalibCont_t& hitVec,
                    const Selector_t<UnCalibCont_t>& selector,
                    std::size_t& hitIdx) const;
  /// @brief Move the layer index towards the possible value or,if the
  ///        layer index is not yet initializes the lyaer index to
  //         the next possible value.
  /// @param strawLayers: List of all straw hits split into the particular layers
  /// @param selector: Delegate method to skip bad hits
  /// @param boundary: Boundary value that the layer index must not cross
  /// @param layerIndex: Mutable reference to the layer index that needs to be moved
  /// @param hitIdx: Mutable reference to the associated hit index inside the layer
  /// @param moveForward: Flag toggling whether the layer index shall be incremented or
  ///                     decremented.
  template <CompositeSpacePointContainer UnCalibCont_t>
  bool nextLayer(const StrawLayers_t<UnCalibCont_t>& strawLayers,
                 const Selector_t<UnCalibCont_t>& selector,
                 const std::size_t boundary,
                 std::optional<std::size_t>& layerIndex, std::size_t& hitIdx,
                 bool moveForward) const;
#ifdef STONJEK

  template <CompositeSpacePointContainer Cont_t,
            CompositeSpacePointSorter<Cont_t> Splitter_t,
            CompositeSpacePointContainer CalibCont_t,
            CompositeSpacePointCalibrator<Cont_t, CalibCont_t> Calibrator_t>
  std::optional<SeedSolution<CalibCont_t>> buildSeed(
      SeedOptions<Cont_t, Splitter_t, CalibCont_t, Calibrator_t>& options)
      const;

  template <CompositeSpacePointContainer Cont_t,
            CompositeSpacePointSorter<Cont_t> Splitter_t,
            CompositeSpacePointContainer CalibCont_t,
            CompositeSpacePointCalibrator<Cont_t, CalibCont_t> Calibrator_t>
  void moveToNextCandidate(SeedOptions<Cont_t, Splitter_t, CalibCont_t,
                                       Calibrator_t>& options) const;
#endif
  /// @brief Construct the final seed parameters by combining the initial
  ///        pattern parameters with the parameter from two circle tangent
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
