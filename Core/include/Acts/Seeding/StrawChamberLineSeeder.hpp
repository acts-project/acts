// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.
#pragma once

#include "Acts/Definitions/TrackParametrization.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/StationSpacePoint.hpp"
#include "Acts/Surfaces/detail/LineHelper.hpp"
#include "Acts/Utilities/CalibrationContext.hpp"

#include <array>
#include <functional>
#include <optional>

namespace Acts {
namespace detail {
/** @brief Abrivation to derference a pointer to a const object */
template <typename T>
using pointer_const_dref_t = std::add_const_t<std::remove_pointer_t<T>>&;

}  // namespace detail
/** @brief Define the concept of the space point measurement sorter. The sorter shall take a collection
 *         of station space points and split them into straw && strip hits. Hits
 * from each category are then subdivided further into the particular detector
 * layers.
 *
 *      A possible implementation of the StationSpacePointSorter needs to have
 * the following attribtues
 *
 *      using SpVec_t =  Standard container satisfiyng the
 * StationSpacePointContainer concept
 *
 *
 *      const std::vector<SpVec_t>& strawHits();
 *      const std::vector<SpVec_t>& stripHits();
 *  Each SpVec_t contains all measurements from a particular detector layer   */
template <typename Splitter_t, typename SpacePointCont_t>
concept StationSpacePointSorter =
    StationSpacePointContainer<SpacePointCont_t> &&
    std::constructible_from<Splitter_t, const SpacePointCont_t&> &&
    requires(Splitter_t sorter) {
      /** @brief Return the straw-hit space point sorted by straw layer */
      {
        sorter.strawHits()
      } -> std::same_as<const std::vector<SpacePointCont_t>&>;
      /** @brief Return the strip-hit  space points sorted by detector layer */
      {
        sorter.stripHits()
      } -> std::same_as<const std::vector<SpacePointCont_t>&>;
    };

/** @brief The StationSpacePointCalibrator refines the StationsSpacePoints by using the points of closest\
 *         approach to the external track. In case, of the straw measurements,
 * the calibrator further provides the derivatives of the electron-drift  <->
 * straw radius relations */

template <typename Calibrator_t, typename UnCalibCont_t, typename CalibCont_t>
concept StationSpacePointCalibrator =
    StationSpacePointContainer<UnCalibCont_t> &&
    StationSpacePointContainer<CalibCont_t> &&

    requires(const Calibrator_t calibrator, const UnCalibCont_t& uncalibCont,
             CalibCont_t& calibCont, const Vector3& trackPos,
             const Vector3& trackDir, const double trackT0,
             const CalibrationContext& ctx,
             detail::pointer_const_dref_t<typename UnCalibCont_t::value_type>
                 singleUnCalibSp,
             detail::pointer_const_dref_t<typename CalibCont_t::value_type>
                 singleCalibSp) {
      /** @brief Calibrate the entire input space point container using the external track parameters
       *  @param ctx: Calibration context to access the calibration constants (Experiment specific)
       *  @param trackPos: Position of the track / segment
       *  @param trackDir: Direction of the track / segment
       *  @param trackT0: Time offset w.r.t. the nominal time of flight calculated by (globPos) / c
       *  @param uncalibCont: Const reference to the calibrated input container to calibrate */
      {
        calibrator.calibrate(ctx, trackPos, trackDir, trackT0, uncalibCont)
      } -> std::same_as<CalibCont_t>;
      /** @brief Append a single measurement to the calibration container using the external track parameters
       *  @param ctx: Calibration context to access the calibration constants (Experiment specific)
       *  @param trackPos: Position of the track / segment
       *  @param trackDir: Direction of the track / segment
       *  @param trackT0: Time offset w.r.t. the nominal time of flight calculated by (globPos) / c
       *  @param singleUnCalibSp: Reference to the measurement to calibrate
       *  @param calibCont: Reference to the mutable calibration container */
      {
        calibrator.calibrate(ctx, trackPos, trackDir, trackT0, singleUnCalibSp,
                             calibCont)
      };
      /** @brief Returns the drift velocity of the straw measurement's radius - time relation
       *         which is defined as the first derivative of the relation.
       *  @param ctx: Calibration context to access the calibration constants (Experiment specific)
       *  @param singleCalibSp: Reference to the calibrated space point */
      { calibrator.driftVelocity(ctx, singleCalibSp) } -> std::same_as<double>;
      /** @brief Returns the drift acceleration of the straw measurement's radius - time relation
       *         which is defined as the second derivative of the relation
       *  @param ctx: Calibration context to access the calibration constants (Experiment specific)
       *  @param singleCalibSp: Reference to the calibrated space point */
      {
        calibrator.driftAcceleration(ctx, singleCalibSp)
      } -> std::same_as<double>;
    };

template <StationSpacePointContainer UnCalibCont_t,
          StationSpacePointSorter<UnCalibCont_t> Splitter_t,
          StationSpacePointContainer CalibCont_t,
          StationSpacePointCalibrator<UnCalibCont_t, CalibCont_t> Calibrator_t>
class StrawChamberLineSeeder {
 public:
  using SeedParam_t = ActsVector<6>;
  struct Config {
    /** @brief Cut on the theta angle */
    std::array<double, 2> thetaRange{0, 179. * UnitConstants::degree};
    /** @brief Cut on the intercept range */
    std::array<double, 2> interceptRange{-20. * UnitConstants::m,
                                         20. * UnitConstants::m};
    /** @brief Upper cut on the hit chi2 w.r.t. seed in order to be associated to the seed*/
    double hitPullCut{5.};
    /** @brief Try at the first time the external seed parameters as candidate */
    bool startWithPattern{false};
    /** @brief How many drift circles may be on a layer to be used for seeding */
    unsigned int busyLayerLimit{2};
    /** @brief How many drift circle hits needs the seed to contain in order to be valid */
    unsigned int nStrawHitCut{3};
    /** @brief Hit cut based on the fraction of collected tube layers.
     *         The seed must pass the tighter of the two requirements.  */
    double nStrawLayHitCut{2. / 3.};
    /** @brief Once a seed with even more than initially required hits is found,
     *         reject all following seeds with less hits */
    bool tightenHitCut{true};
    /** @brief Check whether a new seed candidate shares the same left-right solution with already accepted ones
     *         Reject the seed if it has the same amount of hits */
    bool overlapCorridor{true};
    /** @brief Recalibrate the seed drift circles from the initial estimate  */
    bool recalibSeedCircles{false};
    /** @brief Pointer to the space point calibrator */
    // const ISpacePointCalibrator* calibrator{nullptr};
    /** @brief Toggle whether the seed is rapidly refitted */
    bool fastSeedFit{true};
    /** @brief Toggle whether an initial t0 fit shall be executed */
    bool fastSegFitWithT0{false};
    /** @brief Maximum number of iterations in the fast segment fit */
    unsigned int nMaxIter{100};
    /** @brief Precision cut off in the fast segment fit */
    double precCutOff{1.e-6};
    /** @brief Pointer to the calibrator object */
    const Calibrator_t* calibrator{nullptr};
  };
  /** @brief Seed object returned by the seeder. The seed contains the initial parameter estimate w.r.t
   *         to the central plane surface inside the chamber. *Note* the
   * parameter q/p is set to zero and does not serve any purpose here.
   *
   *         Further the seed contains the list of associated measurements. If
   * the seeder is run in the fast 2D fit mode, the number of iterations w.r.t.
   * this fit procedure is appended as well. */
  struct DriftCircleSeed {
    /** @param Initial straight line seed parameters  */
    SeedParam_t parameters{SeedParam_t::Zero()};
    /** @brief List of calibrated hits on the seed  */
    CalibCont_t measurements{};
    /** @brief Iterations to obtain the seed with fast fit */
    unsigned int nIter{0};
    /** @brief Seed chi2 */
    double chi2{0.};
    /** @brief Number of straw hits on the seed */
    unsigned int nStrawHits{0};
  };
  /** @brief Standard constructor taking grainly externally parameters,
   *         the list of hits to run the seeding on and a generic configuration
   *         object to steer the seeder.
   *  @param grainSeedPars: Initial estimate of the segment line parameters coming from
   *                         e.g. a Hough transform. The parameters may include
   * both directions and intercept. The components, in the non-precision plane
   * are later used to derive the calibrations of the space points using the
   * calibrator
   * @param seedHits: List of straw & strip hits to run the seeding on. Seeds are always generated from
   *                  the straw hits and compatible strip hits are added later
   * once the straw seed passes the configured selection criteria
   * @param cfg: Configuration object to steer the seeder
   * @param logObj: Standard acts logger for later debugging */
  StrawChamberLineSeeder(const SeedParam_t& grainSeedPars,
                         const UnCalibCont_t& seedHits, const Config& cfg,
                         std::unique_ptr<const Logger> logObj);

  /** @brief Remove the copy constructor */
  StrawChamberLineSeeder(const StrawChamberLineSeeder& other) = delete;
  /** @brief Remove the copy assignment */
  StrawChamberLineSeeder& operator=(const StrawChamberLineSeeder& other) =
      delete;

  /** @brief Returns the number of already generated seeds */
  unsigned int numGenerated() const { return m_nGenSeeds; }
  /** @brief Main interface function to fetch a new seed from the seeder. The seeder iterates
   *         through the two most outer straw layers and tries to construct
   * tangent lines to two circles and then to collect other close-by
   * measurements. If no more seeds may be found, a std::nullopt is returned
   *  @param ctx: Calibration context to access the calibration constants (Experiment specific) */
  std::optional<DriftCircleSeed> generateSeed(const CalibrationContext& ctx);

 private:
  /** @brief typedef to the underlying space point object */
  using UncalibSp_t = typename UnCalibCont_t::value_type;

  /** @brief The SeedSolution bookkeeps the straw space points used.
   *         In principle, it's nothing else than a vector of pointers.
   *         However, if the client decides to provide a vector of unique_ptrs,
   *         we need to wrap the objects around a std::reference_wrapper */
  template <typename Trait_t>
  struct UncalibSpBook_t {
    using type = Trait_t;
  };
  template <SmartPointerConcept Trait_t>
  struct UncalibSpBook_t<Trait_t> {
    using type = std::reference_wrapper<const Trait_t>;
  };

  using SpInSeedBookCont_t =
      std::vector<typename UncalibSpBook_t<UncalibSp_t>::type>;
  /** @brief Parameters from an external seed  */
  SeedParam_t m_coarsePars{SeedParam_t::Zero()};
  /** @brief Instance of the layer hit sorter */
  const Splitter_t m_hitLayers;
  Config m_cfg{};
  std::unique_ptr<const Logger> m_logger{};
  /** @brief Sign combination to draw the 4 lines tangent to 2 drift circles */
  using SignCombo_t = std::array<int, 2>;
  constexpr static std::array<SignCombo_t, 4> s_signCombos{
      std::array{1, 1}, std::array{1, -1}, std::array{-1, -1},
      std::array{-1, 1}};

  /** @brief Cache of all solutions seen thus far */
  struct SeedSolution {
    /** @brief: Theta of the line */
    double theta{0.};
    /** @brief Intersecpt of the line */
    double Y0{0.};
    /** @brief: Uncertainty on the slope*/
    double dTheta{0.};
    /** @brief: Uncertainty on the intercept */
    double dY0{0.};
    /** @brief Used hits in the seed */
    SpInSeedBookCont_t seedHits{};
    /** @brief Vector of radial signs of the valid hits */
    std::vector<int> solutionSigns{};
    /** @brief Stringstream output operator */
    friend std::ostream& operator<<(std::ostream& ostr,
                                    const SeedSolution& sol) {
      return sol.print(ostr);
    }
    /** @brief Prints the seed solution to the screen */
    std::ostream& print(std::ostream& ostr) const;
  };

  /** @brief Return the reference to the logger */
  const Logger& logger() const { return *m_logger; }

  /** @brief Prepares the generator to generate the seed from the next pair of drift circles */
  void moveToNextCandidate();

  /** @brief Tries to build the seed from the two hits. Fails if the solution is invalid
   *         or if the seed has already been built before
   *  @param ctx: Calibration context to access the calibration constants (Experiment specific)
   *  @param topHit: Hit candidate from the upper layer
   *  @param bottomHit: Hit candidate from the lower layer
   *  @param sign: Object encoding whether the tangent is left / right  */
  std::optional<DriftCircleSeed> buildSeed(const CalibrationContext& ctx,
                                           const UncalibSp_t& topHit,
                                           const UncalibSp_t& bottomHit,
                                           const SignCombo_t& signs);

  /** @brief Auxillary struct to calculate the fast-fit constants */
  struct SeedFitAuxilliaries {
    /** @brief Tube position center weigthed with inverse covariances */
    Vector3 centerOfGrav{Vector3::Zero()};
    /** @brief Vector of inverse covariances */
    std::vector<double> invCovs{};
    /** @brief Vector of drfit signs */
    std::vector<int> driftSigns{};
    /** @brief Covariance norm */
    double covNorm{0.};
    /** @brief Expectation value of T_{z}^{2} - T_{y}^{2} */
    double T_zzyy{0.};
    /** @brief Expectation value of T_{y} * T_{z} */
    double T_yz{0.};
    /** @brief Expectation value of T_{z} * r  */
    double T_rz{0.};
    /** @brief Expectation value of T_{y} * r  */
    double T_ry{0.};
    /** @brief Prediced y0 given as the expection value of the radii
     *         divided by the inverse covariance sum. */
    double fitY0{0.};
  };
  /** @brief Extension of the auxillary fit constants needed for the
   *         seed refinement when T0 is floating s*/
  struct SeedFitAuxWithT0 : public SeedFitAuxilliaries {
    /** @brief Constructor */
    SeedFitAuxWithT0(SeedFitAuxilliaries&& parent)
        : SeedFitAuxilliaries{std::move(parent)} {}
    /** @brief Expectation value of T_{y} * v */
    double T_vy{0.};
    /** @brief Expectation value of T_{z} * v */
    double T_vz{0.};
    /** @brief Expectation value of T_{y} * a */
    double T_ay{0.};
    /** @brief Expectation value of T_{z} * a */
    double T_az{0.};
    /** @brief Expectation value of r * v */
    double R_vr{0.};
    /** @brief Expectation value of v * v */
    double R_vv{0.};
    /** @brief Expectation value of r * a */
    double R_va{0.};
    /** @brief First derivative of the fitted Y0 */
    double fitY0Prime{0.};
    /** @brief Second derivative of the ftted Y0 */
    double fitY0TwoPrime{0.};
  };

  /** @brief Executes a fast fit to refine the seed parameters
   *         from the parsed drift circles. The angle & the intercept
   *         are estimated by a 1D fit without recalibrating the drift circless
   *  @param candidateSeed: Candidate seed to be refined */
  void fitDriftCircles(DriftCircleSeed& candidateSeed) const;

  /** @brief Executes a fast fit to refine the seed parameters from the
   *         parsed drift circles. The three parameters theta, intercept & t0
   *         are estimated by a 2D fit.
   *  @param ctx: Calibration context to access the calibration constants (Experiment specific)
   *  @param candidateSeed: Candidate seed to be refined*/
  void fitDriftCirclesWithT0(const CalibrationContext& ctx,
                             DriftCircleSeed& candidateSeed) const;

  /** @brief Estimates the auxillary parameters needed to fastly fit the
   *         seed parameters to the drift circle measurements ignoring the time
   * component
   *  @param candidateSeed: Reference to the drift circle candidate seed */
  SeedFitAuxilliaries estimateAuxillaries(
      const DriftCircleSeed& candidateSeed) const;
  /** @brief Estimates the auxillary parameters needed to fastly fit the
   *         seed parameters to the drift circle measurements including the time component
   *  @param ctx: Calibration context to access the calibration constants (Experiment specific)
   *  @param candidateSeed: Reference to the drift circle candidate seed */
  SeedFitAuxWithT0 estimateAuxillaries(
      const CalibrationContext& ctx,
      const DriftCircleSeed& candidateSeed) const;

  /** @brief Combines the seed parameters in the precision plane with the
   *         two parametes in the non bending direction. The precision angle
   *         and offset is taken from the estimate but taking a potential phi
   *         convolution into account
   *  @param y0: Offset in eBoundLoc1
   *  @param theta: Precision angle
   *  @param outPars: Reference to the parameters where the result is written to*/
  void combineWithExternalSeed(const double y0, const double theta,
                               SeedParam_t& outPars) const;
  /** @brief Considered layer to pick the top drift circle from*/
  std::size_t m_upperLayer{0};
  /** @brief Considered layer to pick the bottom drift circle from*/
  std::size_t m_lowerLayer{0};
  /** @brief Explicit hit to pick in the selected bottom layer */
  std::size_t m_lowerHitIndex{0};
  /** @brief Explicit hit to pick in the selected top layer */
  std::size_t m_upperHitIndex{0};
  /** @brief Index of the left-right ambiguity between the circles */
  std::size_t m_signComboIndex{0};
  /** @brief Vector caching equivalent solutions to avoid double seeding */
  std::vector<SeedSolution> m_seenSolutions{};
  /** Counter on how many seeds have been generated */
  unsigned int m_nGenSeeds{0};
};

}  // namespace Acts
#include "Acts/Seeding/StrawChamberLineSeeder.ipp"