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
             const UnCalibSp_t& testSp, CalibCont_t& seedContainer,
             const std::size_t lowerLayer, const std::size_t upperLayer) {
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
        filler.candidateChi2(cctx, pos, dir, t0, testSp)
      } -> std::same_as<double>;
      /// @brief Returns the radius of the straw tube in which the space point
      ///        is recorded
      /// @param testSp: Reference to the straw space point
      { filler.strawRadius(testSp) } -> std::same_as<double>;
      /// @brief Creates a new empty container to construct a new segment seed
      /// @param cctx: Calibration context in case that the container shall be
      ///              part of a predefined memory block
      { filler.newContainer(cctx) } -> std::same_as<CalibCont_t>;
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
      /// @brief Helper method to send a stop signal to the line seeder, if for instance,
      ///        the two layers are too close to each other. The method is
      ///        called after every layer update. If true is returned no seeds
      ///        are produced further
      /// @param lowerLayer: Index of the current lower straw layer
      /// @param upperLayer: Index of the current upper straw layer
      { filler.stopSeeding(lowerLayer, upperLayer) } -> std::same_as<bool>;
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
///                       CompositeSpacePointContainer concept
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

/// @brief Concept of the interface for the auxiliary class such that the
///        CompositeSpacePointLineSeeder can construct segment seeds
template <typename SeedAuxiliary_t, typename UnCalibCont_t,
          typename CalibCont_t>
concept CompSpacePointSeederDelegate =
    CompositeSpacePointSeedFiller<
        SeedAuxiliary_t, RemovePointer_t<typename UnCalibCont_t::value_type>,
        CalibCont_t> &&
    CompositeSpacePointSorter<SeedAuxiliary_t, UnCalibCont_t>;

}  // namespace detail

/// @brief Initial line parameters from a pattern recognition like
///        the Hough transform are often not suitable for a line fit
///        as the resolution of the hough bins usually exceeds the size
///        of the straws.
///        The CompositeSpacePointLineSeeder refines the parameters
///        and the selected measurements such that both become
///        candidates for a stright line fit. The user needs to
///        split the straw measurements per logical straw layer.
///        Further, the interface needs to provide some auxiliary
///        methods to interact with an empty space point container
///        & to calculate a calbrated candidate pull.
///        From these ingredients, the `CompositeSpacePointLineSeeder`
///        iterates from the outermost layers at both ends and tries
///        to construct new candidates. Tangent lines are constructed
///        to a pair of circles from each seeding layer and then straw
///        measurements from the other layers are tried to be added
///        onto the line. If the number of straws exceed the threshold,
///        compatible strip measurements from each strip layer are added.
class CompositeSpacePointLineSeeder {
 public:
  /// @brief Use the assignment of the parameter indices from the CompSpacePointAuxiliaries
  using ParIdx = detail::CompSpacePointAuxiliaries::FitParIndex;
  /// @brief Use the assignment of the parameter indices from the CompSpacePointAuxiliaries
  using CovIdx = detail::CompSpacePointAuxiliaries::ResidualIdx;
  /// @brief Use the vector from the CompSpacePointAuxiliaires
  using Vector = detail::CompSpacePointAuxiliaries::Vector;
  /// @brief Vector containing the 5 straight segment line parameters
  using SeedParam_t = std::array<double, toUnderlying(ParIdx::nPars)>;
  /// @brief Abrivation of the straight line. The first element is the
  ///        reference position and the second element is the direction
  using Line_t = std::pair<Vector, Vector>;

  /// @brief Configuration of the cuts to sort out generated
  ///        seeds with poor quality.
  struct Config {
    /// @brief Cut on the theta angle
    std::array<double, 2> thetaRange{0., 0.};
    /// @brief Cut on the intercept range
    std::array<double, 2> interceptRange{0., 0.};
    /// @brief Upper cut on the hit chi2 w.r.t. seed in order to be associated to the seed
    double hitPullCut{5.};
    /// @brief How many drift circles may be on a layer to be used for seeding
    std::size_t busyLayerLimit{2};
    /// @brief Layers may contain measurements with bad hits and hence the
    bool busyLimitCountGood{true};
    /// @brief Try at the first time the external seed parameters as candidate
    bool startWithPattern{false};
    /// @brief Use explicitly the line distance and the driftRadius to calculate
    ///        the pull from the seed line to the space point.
    bool useSimpleStrawPull{true};
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
  /// @brief Return the configuration object of the seeder
  const Config& config() const { return m_cfg; }

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
    /// @brief Default destructor
    virtual ~TwoCircleTangentPars() = default;

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
                          const Splitter_t& layerSorter)
        : TwoCircleTangentPars{pars}, m_splitter{layerSorter} {}

    /// @brief Abrivation of the underlying space point reference
    using SpacePoint_t = ConstDeRef_t<typename UnCalibCont_t::value_type>;
    /// @brief Helper function to calculate the straw signs of the seed hits
    ///        cached by this solution w.r.t. an external line
    /// @param seedPos: Reference point of the segment line
    /// @param seedDir: Direction of the segment line
    std::vector<int> leftRightAmbiguity(const Vector& seedPos,
                                        const Vector3& seedDir) const;

    /// @brief Returns the number of hits cached by the seed
    std::size_t size() const { return m_seedHits.size(); }
    /// @brief Returns the i-th seed hit
    /// @param idx: Index of the hit to to return
    SpacePoint_t getHit(const std::size_t idx) const;
    /// @brief Appends a new seed hit to the solution
    void append(const std::size_t layIdx, const std::size_t hitIdx);
    /// @brief Vector of the associate left-rignt ambiguities
    std::vector<int> solutionSigns{};
    ///@brief Number of good straw measurements
    std::size_t nStrawHits{0ul};

   private:
    /// @brief Pointer to the space point per layer splitter to gain access to the
    ///        input space point container
    const Splitter_t& m_splitter;
    /// @brief Set of hits collected onto the seed. For each element
    ///        the first index represents the layer &
    ///        the second one the particular hit in that layer
    std::vector<std::pair<std::size_t, std::size_t>> m_seedHits{};
    /// @brief Prints the seed solution to the screen
    /// @param ostr: Mutable reference to the stream to print to
    void print(std::ostream& ostr) const final;
  };

 public:
  /// @brief Helper struct to pack the parameters and the associated
  ///        measurements into a common object. Returned by the
  ///        central nextSeed method (cf. below)
  template <CompositeSpacePointContainer contType_t>
  struct SegmentSeed {
    /// @brief Constructor taking the seed parameters &&
    ///        a new hit container
    /// @param _pars: The seed line parameter
    /// @param _hits  A new empty container to be filled
    explicit SegmentSeed(SeedParam_t _pars, contType_t&& _hits) noexcept
        : parameters{_pars}, hits{std::move(_hits)} {}
    /// @brief Seed line parameters
    SeedParam_t parameters;
    /// @brief Collection of hits
    contType_t hits;
  };

  /// @brief Central auxiliary struct to steer the seeding process.
  ///        First, the user needs to implement the experiment specific
  ///        CompSpacePointSeederDelegate over which the template of the
  ///        SeedingState is then specified. The base class provides the
  ///        straw hits split per logical layer and the SeedingState holds
  ///        the indices to select iteratively a straw measurement each
  ///        from the first and last layer. Also it keeps track of the
  ///        previously constructed seeds.
  template <CompositeSpacePointContainer UncalibCont_t,
            CompositeSpacePointContainer CalibCont_t,
            detail::CompSpacePointSeederDelegate<UncalibCont_t, CalibCont_t>
                Delegate_t>
  struct SeedingState : public Delegate_t {
    /// @brief Declare the public constructor by explicitly forwarding the
    ///        constructor arguments to the (protected) base class constructor
    /// @param initialPars: Initial parameters from an external pattern seed
    ///                     (Needed to combine the precision parameters with a
    ///                     non-precision estimate)
    /// @param args: Arguments to be forwarded to the base class constructor
    template <typename... args_t>
    explicit SeedingState(const SeedParam_t& initialPars, args_t&&... args)
        : Delegate_t{std::forward<args_t>(args)...},
          m_initialPars{initialPars} {}
    /// @brief Stringstream output operator
    friend std::ostream& operator<<(std::ostream& ostr,
                                    const SeedingState& opts) {
      opts.print(ostr);
      return ostr;
    }
    /// @brief Return the number of generated seeds
    std::size_t nGenSeeds() const { return m_seenSolutions.size(); }
    /// @brief Returns the pattern parameters
    const SeedParam_t& initialParameters() const { return m_initialPars; }
    /// @brief Grant the embedding class access to the private members
    friend CompositeSpacePointLineSeeder;

   private:
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
    /// @brief Number of minimum straw hits a seed must have
    std::size_t m_nStrawCut{0ul};
    /// @brief Flag toggling whether the upper of the lower layer shall be moved
    bool m_moveUpLayer{true};
    /// @brief Flag toggling whether the pattern parameters shall be returned as
    ///        first seed
    bool m_patternSeedProduced{false};
    /// @brief Prints the seed solution to the screen
    void print(std::ostream& ostr) const;
    /// @brief Estimated parameters from pattern
    SeedParam_t m_initialPars{};
  };
  /// @brief Main interface method provided by the SeederClass. The user instantiates
  ///        a SeedingState object containing all the straw hit candidates from
  ///        which the seed shall be constructed. Then, the nextSeed() returns
  ///        the next best seed candidate which can then be fitted. The user
  ///        continues to call the method until a nullopt is returned.
  /// @param cctx: Experiment specific calibration context to be piped back to the
  ///              caller such that the space points may be calibrated during
  ///              the seeding process.
  /// @param state: Mutable reference to the SeedingState object from which all the
  ///               segment seeds are constructed.
  template <CompositeSpacePointContainer UncalibCont_t,
            CompositeSpacePointContainer CalibCont_t,
            detail::CompSpacePointSeederDelegate<UncalibCont_t, CalibCont_t>
                Delegate_t>
  std::optional<SegmentSeed<CalibCont_t>> nextSeed(
      const CalibrationContext& cctx,
      SeedingState<UncalibCont_t, CalibCont_t, Delegate_t>& state) const;

 private:
  /// @brief Reference to the logger object
  const Logger& logger() const { return *m_logger; }
  /// @brief Abrivation of the selector delegate to skip invalid straw hits in the seed
  template <CompositeSpacePointContainer Cont_t>
  using Selector_t = Delegate<bool(ConstDeRef_t<typename Cont_t::value_type>)>;
  /// @brief Abrivation of the split hit containers
  template <CompositeSpacePointContainer Cont_t>
  using StrawLayers_t = std::vector<Cont_t>;
  /// @brief Counts the number of hits inside the container. Depending on whether
  ///        the busyLimitCountGood flag is true, bad hits are not considered
  /// @param container: Reference to the container which size is to be evaluated
  /// @param selector: Delegate method to skip bad bad hits
  template <CompositeSpacePointContainer Cont_t>
  std::size_t countHits(const Cont_t& container,
                        const Selector_t<Cont_t>& selector) const;
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
  /// @brief Move the layer and hit indices inside the state towards the next candidate
  ///        pair. First, the L-R ambiguities are incremented, then it is
  ///        searched for the next pair inside the lower && upper layer pair.
  ///        Finally, the indices are moved towards the next layer
  /// @param selector: Delegate method to skip bad hits
  /// @param state: Mutable reference to the SeedingState object carring the state indices
  template <CompositeSpacePointContainer UncalibCont_t,
            CompositeSpacePointContainer CalibCont_t,
            detail::CompSpacePointSeederDelegate<UncalibCont_t, CalibCont_t>
                Delegate_t>
  void moveToNextCandidate(
      const Selector_t<UncalibCont_t>& selector,
      SeedingState<UncalibCont_t, CalibCont_t, Delegate_t>& state) const;
  /// @brief Attempts to construct the next seed from the given configuration of
  ///        seed circles. The seed needs to contain a minimum number of other
  ///        straw hits and there must be no other previously constructed seed
  ///        with the same Left-Right solution
  /// @param cctx: Calibration context to be piped to the experiment's implementation
  ///              such that conditions data access becomes possible
  /// @param selector: Delegate method to skip bad hits
  /// @param state: Mutable reference to the SeedingState object carring the state indices
  template <CompositeSpacePointContainer UncalibCont_t,
            CompositeSpacePointContainer CalibCont_t,
            detail::CompSpacePointSeederDelegate<UncalibCont_t, CalibCont_t>
                Delegate_t>
  std::optional<SegmentSeed<CalibCont_t>> buildSeed(
      const CalibrationContext& cctx, const Selector_t<UncalibCont_t>& selector,
      SeedingState<UncalibCont_t, CalibCont_t, Delegate_t>& state) const;
  /// @brief Checks whether the new seed candidate passes the quality cuts on
  ///        the number of good straw hits and whether it is not within the
  ///        same overlap corridor as previously produced seeds
  /// @param tangentSeed: Pair of reference position & direction constructed
  ///                     from the two line tangent seed
  /// @param newSolution: The new seed solution that's to be tested
  /// @param state: The cache carrying the already produced solutions
  template <CompositeSpacePointContainer UncalibCont_t,
            CompositeSpacePointContainer CalibCont_t,
            detail::CompSpacePointSeederDelegate<UncalibCont_t, CalibCont_t>
                Delegate_t>
  bool passSeedCuts(
      const Line_t& tangentSeed,
      SeedSolution<UncalibCont_t, Delegate_t>& newSolution,
      SeedingState<UncalibCont_t, CalibCont_t, Delegate_t>& state) const;
  /// @brief Converts the accepted seed solution to the segment seed returned by
  ///        nextSeed and adds the strip measurements to the seed. The solution
  ///        is then appended to the state
  /// @param cctx: Calibration context to be piped to the experiment's implementation
  ///              such that conditions data access becomes possible
  /// @param tangentSeed: Position and direction constructed from the current tangent seed
  /// @param state: Mutable reference to the state from which the strip measurements are drawn
  ///               and to which the newSolution is then appended
  /// @param newSolution: Current tangent seed solution object holding the straw measurements
  ///                     to be put onto the seed.
  template <CompositeSpacePointContainer UncalibCont_t,
            CompositeSpacePointContainer CalibCont_t,
            detail::CompSpacePointSeederDelegate<UncalibCont_t, CalibCont_t>
                Delegate_t>
  SegmentSeed<CalibCont_t> consructSegmentSeed(
      const CalibrationContext& cctx, const Line_t& tangentSeed,
      SeedingState<UncalibCont_t, CalibCont_t, Delegate_t>& state,
      SeedSolution<UncalibCont_t, Delegate_t>&& newSolution) const;
  /// @brief Construct the final seed parameters by combining the initial
  ///        pattern parameters with the parameter from two circle tangent
  /// @param tangentSeed: Pair of reference position & direction constructed
  ///                     from the two line tangent seed
  /// @param patternParams: Parameter estimate from the hit pattern
  SeedParam_t combineWithPattern(const Line_t& tangentSeed,
                                 const SeedParam_t& patternParams) const;
  /// @brief Constructs a line from the parsed seed parameters. The
  ///        first element is the reference point && the second one
  ///        is the direction
  /// @param pars: Reference to the line parameters from which the line is created
  Line_t makeLine(const SeedParam_t& pars) const;
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
