// This file is part of the Acts project.
//
// Copyright (C) 2021 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

// Workaround for building on clang+libstdc++
#include "Acts/Utilities/detail/ReferenceWrapperAnyCompat.hpp"

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Definitions/Units.hpp"
#include "Acts/EventData/MultiComponentBoundTrackParameters.hpp"
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"
#include "Acts/Propagator/DefaultExtension.hpp"
#include "Acts/Propagator/DenseEnvironmentExtension.hpp"
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/EigenStepperError.hpp"
#include "Acts/Propagator/StepperExtensionList.hpp"
#include "Acts/Propagator/detail/Auctioneer.hpp"
#include "Acts/Propagator/detail/SteppingHelper.hpp"
#include "Acts/Utilities/Intersection.hpp"
#include "Acts/Utilities/Result.hpp"

#include <algorithm>
#include <cmath>
#include <functional>
#include <limits>
#include <numeric>
#include <sstream>
#include <vector>

#include "MultiStepperError.hpp"

// #define PRINT_STEPSIZE_CHANGE

namespace Acts {

using namespace Acts::UnitLiterals;

/// @brief Reducer struct for the Loop MultiEigenStepper which reduces the
/// multicomponent state to simply by summing the weighted values
struct WeightedComponentReducerLoop {
  template <typename component_t>
  static Vector3 toVector3(const std::vector<component_t>& comps,
                           const FreeIndices i) {
    return std::accumulate(
        begin(comps), end(comps), Vector3{Vector3::Zero()},
        [i](const auto& sum, const auto& cmp) -> Vector3 {
          return sum + cmp.weight * cmp.state.pars.template segment<3>(i);
        });
  }

  template <typename stepper_state_t>
  static Vector3 position(const stepper_state_t& s) {
    return toVector3(s.components, eFreePos0);
  }

  template <typename stepper_state_t>
  static Vector3 direction(const stepper_state_t& s) {
    return toVector3(s.components, eFreeDir0).normalized();
  }

  template <typename stepper_state_t>
  static ActsScalar momentum(const stepper_state_t& s) {
    return std::accumulate(
        begin(s.components), end(s.components), ActsScalar{0.},
        [](const auto& sum, const auto& cmp) -> ActsScalar {
          return sum +
                 cmp.weight * (1 / (cmp.state.pars[eFreeQOverP] / cmp.state.q));
        });
  }

  template <typename stepper_state_t>
  static ActsScalar charge(const stepper_state_t& s) {
    return std::accumulate(begin(s.components), end(s.components),
                           ActsScalar{0.},
                           [](const auto& sum, const auto& cmp) -> ActsScalar {
                             return sum + cmp.weight * cmp.state.q;
                           });
  }

  template <typename stepper_state_t>
  static ActsScalar time(const stepper_state_t& s) {
    return std::accumulate(
        begin(s.components), end(s.components), ActsScalar{0.},
        [](const auto& sum, const auto& cmp) -> ActsScalar {
          return sum + cmp.weight * cmp.state.pars[eFreeTime];
        });
  }

  template <typename stepper_state_t>
  static FreeVector pars(const stepper_state_t& s) {
    return std::accumulate(begin(s.components), end(s.components),
                           FreeVector{FreeVector::Zero()},
                           [](const auto& sum, const auto& cmp) -> FreeVector {
                             return sum + cmp.weight * cmp.state.pars;
                           });
  }

  template <typename stepper_state_t>
  static FreeVector cov(const stepper_state_t& s) {
    return std::accumulate(begin(s.components), end(s.components),
                           FreeMatrix{FreeMatrix::Zero()},
                           [](const auto& sum, const auto& cmp) -> FreeMatrix {
                             return sum + cmp.weight * cmp.state.cov;
                           });
  }
};

struct MaxMomentumReducerLoop {
  template <typename component_t>
  static const auto& maxMomenutmIt(const std::vector<component_t>& cmps) {
    return *std::max_element(cmps.begin(), cmps.end(),
                             [&](const auto& a, const auto& b) {
                               return std::abs(a.state.pars[eFreeQOverP]) >
                                      std::abs(b.state.pars[eFreeQOverP]);
                             });
  }

  template <typename stepper_state_t>
  static Vector3 position(const stepper_state_t& s) {
    return maxMomenutmIt(s.components)
        .state.pars.template segment<3>(eFreePos0);
  }

  template <typename stepper_state_t>
  static Vector3 direction(const stepper_state_t& s) {
    return maxMomenutmIt(s.components)
        .state.pars.template segment<3>(eFreeDir0);
  }

  template <typename stepper_state_t>
  static ActsScalar momentum(const stepper_state_t& s) {
    const auto& cmp = maxMomenutmIt(s.components);
    return 1.0 / (cmp.state.pars[eFreeQOverP] / cmp.state.q);
  }

  template <typename stepper_state_t>
  static ActsScalar charge(const stepper_state_t& s) {
    return maxMomenutmIt(s.components).state.q;
  }

  template <typename stepper_state_t>
  static ActsScalar time(const stepper_state_t& s) {
    return maxMomenutmIt(s.components).state.pars[eFreeTime];
  }

  template <typename stepper_state_t>
  static FreeVector pars(const stepper_state_t& s) {
    return maxMomenutmIt(s.components).state.pars;
  }

  template <typename stepper_state_t>
  static FreeVector cov(const stepper_state_t& s) {
    return maxMomenutmIt(s.components).state.cov;
  }
};

/// @brief Stepper based on the EigenStepper, but handles Multi-Component Tracks
/// (e.g., for the GSF)
template <typename extensionlist_t,
          typename component_reducer_t = WeightedComponentReducerLoop,
          typename auctioneer_t = detail::VoidAuctioneer>
class MultiEigenStepperLoop
    : public EigenStepper<extensionlist_t, auctioneer_t> {
  const LoggerWrapper logger;

  /// @brief Limits the number of steps after at least one component reached the surface
  std::size_t m_stepLimitAfterFirstComponentOnSurface = 50;

 public:
  /// @brief Typedef to the Single-Component Eigen Stepper
  using SingleStepper = EigenStepper<extensionlist_t, auctioneer_t>;

  /// @brief Typedef to the State of the single component Stepper
  using SingleState = typename SingleStepper::State;

  /// @brief Use the definitions from the Single-stepper
  using typename SingleStepper::BoundState;
  using typename SingleStepper::Covariance;
  using typename SingleStepper::CurvilinearState;
  using typename SingleStepper::Jacobian;

  /// @brief The reducer type
  using Reducer = component_reducer_t;

  /// @brief How many components can this stepper manage?
  static constexpr int maxComponents = std::numeric_limits<int>::max();

  struct State {
    State() = delete;

    /// Constructor from the initial bound track parameters
    ///
    /// @tparam charge_t Type of the bound parameter charge
    ///
    /// @param [in] gctx is the context object for the geometry
    /// @param [in] mctx is the context object for the magnetic field
    /// @param [in] par The track parameters at start
    /// @param [in] ndir The navigation direciton w.r.t momentum
    /// @param [in] ssize is the maximum step size
    /// @param [in] stolerance is the stepping tolerance
    ///
    /// @note the covariance matrix is copied when needed
    template <typename charge_t>
    explicit State(
        const GeometryContext& gctx, const MagneticFieldContext& mctx,
        const std::shared_ptr<const MagneticFieldProvider>& bField,
        const MultiComponentBoundTrackParameters<charge_t>& multipars,
        NavigationDirection ndir = forward,
        double ssize = std::numeric_limits<double>::max(),
        double stolerance = s_onSurfaceTolerance)
        : navDir(ndir), geoContext(gctx), magContext(mctx) {
      throw_assert(!multipars.components().empty(),
                   "Empty MultiComponent state");

      for (const auto& [weight, single_component] : multipars.components()) {
        components.push_back(
            {SingleState(gctx, bField->makeCache(mctx), single_component, ndir,
                         ssize, stolerance),
             weight, Intersection3D::Status::reachable});
      }

      if (components.front().state.covTransport) {
        covTransport = true;
      }
    }

    /// The struct that stores the individual components
    struct Component {
      SingleState state;
      ActsScalar weight;
      Intersection3D::Status status;
    };

    /// The components of which the state consists
    std::vector<Component> components;

    /// Required through stepper concept
    /// TODO how can they be interpreted for a Multistepper
    bool covTransport = false;
    Covariance cov = Covariance::Zero();
    NavigationDirection navDir;
    double pathAccumulated = 0.;
    int steps = 0;

    /// geoContext
    std::reference_wrapper<const GeometryContext> geoContext;

    /// MagneticFieldContext
    std::reference_wrapper<const MagneticFieldContext> magContext;

    /// Step-limit counter which limits the number of steps when one component
    /// reached a surface
    std::optional<std::size_t> stepCounterAfterFirstComponentOnSurface;
  };

  /// A state type which can be used for a function which accepts only
  /// single-component states and stepper
  template <typename navigation_t, typename options_t>
  struct SinglePropState {
    SinglePropState(SingleState& s, navigation_t& n, options_t& o,
                    GeometryContext g)
        : stepping(s), navigation(n), options(o), geoContext(g) {}
    SingleState& stepping;
    navigation_t& navigation;
    options_t& options;
    GeometryContext geoContext;
  };

  template <typename navigation_t, typename options_t>
  SinglePropState(SingleState, navigation_t, options_t, GeometryContext)
      -> SinglePropState<navigation_t, options_t>;

  /// A proxy struct which allows access to a single component of the
  /// multi-component state
  struct ComponentProxy {
    const State& m_state;
    typename State::Component& m_cmp;

    auto& status() { return m_cmp.status; }
    auto status() const { return m_cmp.status; }
    auto& weight() { return m_cmp.weight; }
    auto weight() const { return m_cmp.weight; }
    auto& charge() { return m_cmp.state.charge; }
    auto charge() const { return m_cmp.state.charge; }
    auto& pathLength() { return m_cmp.state.pathAccumulated; }
    auto pathLength() const { return m_cmp.state.pathAccumulated; }
    auto& pars() { return m_cmp.state.pars; }
    const auto& pars() const { return m_cmp.state.pars; }
    auto& derivative() { return m_cmp.state.derivative; }
    const auto& derivative() const { return m_cmp.state.derivative; }
    auto& jacTransport() { return m_cmp.state.jacTransport; }
    const auto& jacTransport() const { return m_cmp.state.jacTransport; }
    auto& cov() { return m_cmp.state.cov; }
    const auto& cov() const { return m_cmp.state.cov; }
    auto& jacobian() { return m_cmp.state.jacobian; }
    const auto& jacobian() const { return m_cmp.state.jacobian; }
    auto& jacToGlobal() { return m_cmp.state.jacToGlobal; }
    const auto& jacToGlobal() const { return m_cmp.state.jacToGlobal; }

    template <typename propagator_state_t>
    const auto& singleState(const propagator_state_t& state) const {
      static_assert(
          std::is_same_v<SingleState,
                         decltype(state.stepping.components.front().state)>);
      return SinglePropState{m_cmp.state, state.navigation, state.options,
                             state.geoContext};
    }

    template <typename propagator_state_t>
    auto singleState(propagator_state_t& state) {
      static_assert(
          std::is_same_v<SingleState,
                         decltype(state.stepping.components.front().state)>);

      return SinglePropState(m_cmp.state, state.navigation, state.options,
                             state.geoContext);
    }

    const auto& singleStepper(const MultiEigenStepperLoop& stepper) const {
      return static_cast<const SingleStepper&>(stepper);
    }

    Result<BoundState> boundState(const Surface& surface, bool transportCov) {
      if (m_cmp.status == Intersection3D::Status::onSurface) {
        return detail::boundState(m_state.geoContext, cov(), jacobian(),
                                  jacTransport(), derivative(), jacToGlobal(),
                                  pars(), m_state.covTransport && transportCov,
                                  m_cmp.state.pathAccumulated, surface);
      } else {
        return MultiStepperError::ComponentNotOnSurface;
      }
    }

    void update(const FreeVector& freeParams, const BoundVector& boundParams,
                const Covariance& covariance, const Surface& surface) {
      m_cmp.state.pars = freeParams;
      m_cmp.state.cov = covariance;
      m_cmp.state.jacToGlobal =
          surface.boundToFreeJacobian(m_state.geoContext, boundParams);
    }
  };

  struct ConstComponentProxy {
    const State& m_state;
    const typename State::Component& m_cmp;

    auto status() const { return m_cmp.status; }
    auto weight() const { return m_cmp.weight; }
    auto charge() const { return m_cmp.state.charge; }
    auto pathLength() const { return m_cmp.state.pathAccumulated; }
    const auto& pars() const { return m_cmp.state.pars; }
    const auto& derivative() const { return m_cmp.state.derivative; }
    const auto& jacTransport() const { return m_cmp.state.jacTransport; }
    const auto& cov() const { return m_cmp.state.cov; }
    const auto& jacobian() const { return m_cmp.state.jacobian; }
    const auto& jacToGlobal() const { return m_cmp.state.jacToGlobal; }

    template <typename propagator_state_t>
    const auto& singleState(const propagator_state_t& state) const {
      static_assert(
          std::is_same_v<SingleState,
                         decltype(state.stepping.components.front().state)>);
      return SinglePropState{m_cmp.state, state.navigation, state.options,
                             state.geoContext};
    }

    const auto& singleStepper(const MultiEigenStepperLoop& stepper) const {
      return static_cast<const SingleStepper&>(stepper);
    }
  };

  /// Creates an iterable which can be plugged into a range-based for-loop to
  /// iterate over components
  /// @note Use a for-loop with by-value semantics, since the Iterable returns a
  /// proxy internally holding a reference
  auto componentIterable(State& state) const {
    struct Iterator {
      typename decltype(state.components)::iterator it;
      const State& state;

      // clang-format off
      auto& operator++() { ++it; return *this; }
      auto operator!=(const Iterator& other) const { return it != other.it; }
      auto operator*() const { return ComponentProxy{state, *it}; }
      // clang-format on
    };

    struct Iterable {
      State& state;

      auto begin() { return Iterator{state.components.begin(), state}; }
      auto end() { return Iterator{state.components.end(), state}; }
    };

    return Iterable{state};
  }

  /// Creates an constant iterable which can be plugged into a range-based
  /// for-loop to iterate over components
  /// @note Use a for-loop with by-value semantics, since the Iterable returns a
  /// proxy internally holding a reference
  auto constComponentIterable(const State& state) const {
    struct ConstIterator {
      typename decltype(state.components)::const_iterator it;
      const State& state;

      // clang-format off
      auto& operator++() { ++it; return *this; }
      auto operator!=(const ConstIterator& other) const { return it != other.it; }
      auto operator*() const { return ConstComponentProxy{state, *it}; }
      // clang-format on
    };

    struct Iterable {
      const State& state;

      auto begin() const {
        return ConstIterator{state.components.cbegin(), state};
      }
      auto end() const { return ConstIterator{state.components.cend(), state}; }
    };

    return Iterable{state};
  }

  /// Get the number of components
  std::size_t numberComponents(const State& state) const {
    return state.components.size();
  }

  /// Remove missed components from the component state
  void removeMissedComponents(State& state) const {
    auto new_end = std::remove_if(
        state.components.begin(), state.components.end(), [](const auto& cmp) {
          return cmp.status == Intersection3D::Status::missed;
        });

    state.components.erase(new_end, state.components.end());
  }

  /// Reset the number of components
  void clearComponents(State& state) const { state.components.clear(); }

  /// Add a component to the Multistepper
  /// @note This function makes no garantuees about how new components are
  /// initialized, it is up to the caller to ensure that all components are
  /// valid in the end.
  /// @note The returned component-proxy is only garantueed to be valid until
  /// the component number is again modified
  template <typename charge_t>
  Result<ComponentProxy> addComponent(
      State& state, const SingleBoundTrackParameters<charge_t>& pars,
      double weight) const {
    state.components.push_back(
        {SingleState(state.geoContext,
                     SingleStepper::m_bField->makeCache(state.magContext), pars,
                     state.navDir),
         weight, Intersection3D::Status::onSurface});

    return ComponentProxy{state, state.components.back()};
  }

  /// Construct and initialize a state
  template <typename charge_t>
  State makeState(std::reference_wrapper<const GeometryContext> gctx,
                  std::reference_wrapper<const MagneticFieldContext> mctx,
                  const MultiComponentBoundTrackParameters<charge_t>& par,
                  NavigationDirection ndir = forward,
                  double ssize = std::numeric_limits<double>::max(),
                  double stolerance = s_onSurfaceTolerance) const {
    return State(gctx, mctx, SingleStepper::m_bField, par, ndir, ssize,
                 stolerance);
  }

  /// Constructor
  MultiEigenStepperLoop(std::shared_ptr<const MagneticFieldProvider> bField,
                        LoggerWrapper l = getDummyLogger())
      : EigenStepper<extensionlist_t, auctioneer_t>(bField), logger(l) {}

  void setOverstepLimit(double oLimit) {
    SingleStepper::m_overstepLimit = oLimit;
  }

  /// Get the field for the stepping, it checks first if the access is still
  /// within the Cell, and updates the cell if necessary.
  ///
  /// @param [in,out] state is the propagation state associated with the track
  ///                 the magnetic field cell is used (and potentially updated)
  /// @param [in] pos is the field position
  Result<Vector3> getField(State& state, const Vector3& pos) const {
    // get the field from the cell
    return SingleStepper::getField(state.components.front().state, pos);
  }

  /// Global particle position accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  Vector3 position(const State& state) const {
    return Reducer::position(state);
  }

  auto position(std::size_t i, const State& state) const {
    return SingleStepper::position(state.components.at(i).state);
  }

  /// Momentum direction accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  Vector3 direction(const State& state) const {
    return Reducer::direction(state);
  }

  auto direction(std::size_t i, const State& state) const {
    return SingleStepper::direction(state.components.at(i).state);
  }

  /// Absolute momentum accessor
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double momentum(const State& state) const { return Reducer::momentum(state); }

  auto momentum(std::size_t i, const State& state) const {
    return SingleStepper::momentum(state.components.at(i).state);
  }

  /// Charge access
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double charge(const State& state) const { return Reducer::charge(state); }

  /// Time access
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double time(const State& state) const { return Reducer::time(state); }

  auto time(std::size_t i, const State& state) const {
    return SingleStepper::time(state.components.at(i).state);
  }

  /// Update surface status
  ///
  /// It checks the status to the reference surface & updates
  /// the step size accordingly
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param surface [in] The surface provided
  /// @param bcheck [in] The boundary check for this status update
  Intersection3D::Status updateSurfaceStatus(
      State& state, const Surface& surface, const BoundaryCheck& bcheck,
      LoggerWrapper /*extLogger*/ = getDummyLogger()) const {
    using Status = Intersection3D::Status;

    std::array<int, 4> counts = {0, 0, 0, 0};

#ifdef PRINT_STEPSIZE_CHANGE
    const std::string before = outputStepSize(state);
#endif

    for (auto& component : state.components) {
      component.status = detail::updateSingleSurfaceStatus<SingleStepper>(
          *this, component.state, surface, bcheck, logger);
      ++counts[static_cast<std::size_t>(component.status)];
    }

    ACTS_VERBOSE("Component status wrt "
                 << surface.geometryId() << " at {"
                 << surface.center(state.geoContext).transpose() << "}:\t"
                 << [&]() {
                      std::stringstream ss;
                      for (auto& component : state.components) {
                        ss << component.status << "\t";
                      }
                      return ss.str();
                    }());

#ifdef PRINT_STEPSIZE_CHANGE
    std::cout << "MultiStepperLoop::updateSurfaceStatus(...):\n"
              << "\tBEFORE" << before << "\n"
              << "\tAFTER" << outputStepSize(state) << std::endl;
#endif

    // Switch on stepCounter if one or more components reached a surface, but
    // some are still in progress of reaching the surface
    if (!state.stepCounterAfterFirstComponentOnSurface &&
        counts[static_cast<std::size_t>(Status::onSurface)] > 0 &&
        counts[static_cast<std::size_t>(Status::reachable)] > 0) {
      state.stepCounterAfterFirstComponentOnSurface = 0;
      ACTS_VERBOSE("started stepCounterAfterFirstComponentOnSurface");
    }

    // This is a 'any_of' criterium. As long as any of the components has a
    // certain state, this determines the total state (in the order of a
    // somewhat importance)
    if (counts[static_cast<std::size_t>(Status::reachable)] > 0) {
      return Status::reachable;
    } else if (counts[static_cast<std::size_t>(Status::onSurface)] > 0) {
      state.stepCounterAfterFirstComponentOnSurface.reset();
      return Status::onSurface;
    } else if (counts[static_cast<std::size_t>(Status::unreachable)] > 0) {
      return Status::unreachable;
    } else {
      return Status::missed;
    }
  }

  /// Update step size
  ///
  /// This method intersects the provided surface and update the navigation
  /// step estimation accordingly (hence it changes the state). It also
  /// returns the status of the intersection to trigger onSurface in case
  /// the surface is reached.
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param oIntersection [in] The ObjectIntersection to layer, boundary, etc
  /// @param release [in] boolean to trigger step size release
  template <typename object_intersection_t>
  void updateStepSize(State& state, const object_intersection_t& oIntersection,
                      bool release = true) const {
#ifdef PRINT_STEPSIZE_CHANGE
    const std::string before = outputStepSize(state);
#endif
    const Surface& surface = *oIntersection.representation;

    for (auto& component : state.components) {
      auto intersection = surface.intersect(
          component.state.geoContext, SingleStepper::position(component.state),
          SingleStepper::direction(component.state), true);

      // We don't know whatever was done to manipulate the intersection before
      // (e.g. in Layer.ipp:240), so we trust and just adjust the sign
      if (std::signbit(oIntersection.intersection.pathLength) !=
          std::signbit(intersection.intersection.pathLength)) {
        intersection.intersection.pathLength *= -1;
      }

      if (std::signbit(oIntersection.alternative.pathLength) !=
          std::signbit(intersection.alternative.pathLength)) {
        intersection.alternative.pathLength *= -1;
      }

      SingleStepper::updateStepSize(component.state, intersection, release);
    }

#ifdef PRINT_STEPSIZE_CHANGE
    std::cout << "MultiStepperLoop::updateStepSize(...):\n"
              << "\tBEFORE" << before << "\n"
              << "\tAFTER" << outputStepSize(state) << std::endl;
#endif
  }

  /// Set Step size - explicitely with a double
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  /// @param stepSize [in] The step size value
  /// @param stype [in] The step size type to be set
  void setStepSize(State& state, double stepSize,
                   ConstrainedStep::Type stype = ConstrainedStep::actor,
                   bool release = true) const {
#ifdef PRINT_STEPSIZE_CHANGE
    const std::string before = outputStepSize(state);
#endif

    for (auto& component : state.components) {
      SingleStepper::setStepSize(component.state, stepSize, stype, release);
    }

#ifdef PRINT_STEPSIZE_CHANGE
    std::cout << "MultiStepperLoop::setStepSize(...):\n"
              << "\tBEFORE" << before << "\n"
              << "\tAFTER" << outputStepSize(state) << std::endl;
#endif
  }

  /// Get the step size
  /// TODO add documentation
  double getStepSize(const State& state, ConstrainedStep::Type stype) const {
    return SingleStepper::getStepSize(
        std::min_element(begin(state.components), end(state.components),
                         [this, stype](const auto& a, const auto& b) {
                           return SingleStepper::getStepSize(a.state, stype) <
                                  SingleStepper::getStepSize(b.state, stype);
                         })
            ->state,
        stype);
  }

  /// Release the Step size
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  void releaseStepSize(State& state) const {
#ifdef PRINT_STEPSIZE_CHANGE
    const std::string before = outputStepSize(state);
#endif

    for (auto& component : state.components) {
      SingleStepper::releaseStepSize(component.state);
    }

#ifdef PRINT_STEPSIZE_CHANGE
    std::cout << "MultiStepperLoop::releaseStepSize(...):\n"
              << "\tBEFORE" << before << "\n"
              << "\tAFTER" << outputStepSize(state) << std::endl;
#endif
  }

  /// Output the Step Size - single component
  ///
  /// @param state [in,out] The stepping state (thread-local cache)
  std::string outputStepSize(const State& state) const {
    std::stringstream ss;
    for (const auto& component : state.components)
      ss << component.state.stepSize.toString() << " || ";

    return ss.str();
  }

  /// Overstep limit
  ///
  /// @param state [in] The stepping state (thread-local cache)
  double overstepLimit(const State& state) const {
    // A dynamic overstep limit could sit here
    return SingleStepper::overstepLimit(state.components.front().state);
  }

  /// @brief Resets the state
  ///
  /// @param [in, out] state State of the stepper
  /// @param [in] boundParams Parameters in bound parametrisation
  /// @param [in] cov Covariance matrix
  /// @param [in] surface The reference surface of the bound parameters
  /// @param [in] navDir Navigation direction
  /// @param [in] stepSize Step size
  void resetState(
      State& state, const BoundVector& boundParams, const BoundSymMatrix& cov,
      const Surface& surface, const NavigationDirection navDir = forward,
      const double stepSize = std::numeric_limits<double>::max()) const {
    for (auto& component : state.components) {
      SingleStepper::resetState(component.state, boundParams, cov, surface,
                                navDir, stepSize);
    }
  }

  /// Create and return the bound state at the current position
  ///
  /// @brief This transports (if necessary) the covariance
  /// to the surface and creates a bound state. It does not check
  /// if the transported state is at the surface, this needs to
  /// be guaranteed by the propagator
  ///
  /// @param [in] state State that will be presented as @c BoundState
  /// @param [in] surface The surface to which we bind the state
  /// @param [in] transportCov Flag steering covariance transport
  ///
  /// @return A bound state:
  ///   - the parameters at the surface
  ///   - the stepwise jacobian towards it (from last bound)
  ///   - and the path length (from start - for ordering)
  Result<BoundState> boundState(State& state, const Surface& surface,
                                bool transportCov = true) const {
    if (numberComponents(state) == 1) {
      return SingleStepper::boundState(state.components.front().state, surface,
                                       transportCov);
    } else {
      std::vector<std::pair<double, BoundState>> bs_vec;
      double accumulatedPathLength = 0.0;
      int failedBoundTransforms = 0;

      for (auto i = 0ul; i < numberComponents(state); ++i) {
        auto bs = SingleStepper::boundState(state.components[i].state, surface,
                                            transportCov);

        if (bs.ok()) {
          bs_vec.push_back({state.components[i].weight, *bs});
          accumulatedPathLength += std::get<double>(*bs);
        } else {
          failedBoundTransforms++;
        }
      }

      const auto [params, cov] = detail::combineComponentRange(
          bs_vec.begin(), bs_vec.end(), [&](const auto& wbs) {
            const auto& bp = std::get<BoundTrackParameters>(wbs.second);
            return std::tie(wbs.first, bp.parameters(), bp.covariance());
          });

      if (failedBoundTransforms > 0) {
        ACTS_WARNING("Multi component bound state: "
                     << failedBoundTransforms << " of "
                     << numberComponents(state) << " transforms failed");
      }

      // TODO Jacobian for multi component state not defined really?
      return BoundState{
          BoundTrackParameters(surface.getSharedPtr(), params, cov),
          Jacobian::Zero(), accumulatedPathLength / bs_vec.size()};
    }
  }

  /// Create and return a curvilinear state at the current position
  ///
  /// @brief This transports (if necessary) the covariance
  /// to the current position and creates a curvilinear state.
  ///
  /// @param [in] state State that will be presented as @c CurvilinearState
  /// @param [in] transportCov Flag steering covariance transport
  ///
  /// @return A curvilinear state:
  ///   - the curvilinear parameters at given position
  ///   - the stepweise jacobian towards it (from last bound)
  ///   - and the path length (from start - for ordering)
  CurvilinearState curvilinearState(State& state,
                                    bool transportCov = true) const {
    // TODO this is not correct, but not needed right now somewhere...
    return SingleStepper::curvilinearState(state.components.front().state,
                                           transportCov);
  }

  /// Method to update a stepper state to the some parameters
  ///
  /// @param [in,out] state State object that will be updated
  /// @param [in] pars Parameters that will be written into @p state
  /// TODO is this function useful for a MultiStepper?
  void update(State& state, const FreeVector& parameters,
              const Covariance& covariance) const {
    for (auto& component : state.components) {
      SingleStepper::update(component.state, parameters, covariance);
    }
  }

  /// Method to update momentum, direction and p
  ///
  /// @param [in,out] state State object that will be updated
  /// @param [in] uposition the updated position
  /// @param [in] udirection the updated direction
  /// @param [in] up the updated momentum value
  /// TODO is this function useful for a MultiStepper?
  void update(State& state, const Vector3& uposition, const Vector3& udirection,
              double up, double time) const {
    for (auto& component : state.components) {
      SingleStepper::update(component.state, uposition, udirection, up, time);
    }
  }

  /// Method to update the components individually
  template <typename component_rep_t>
  void updateComponents(State& state, const std::vector<component_rep_t>& cmps,
                        const Surface& surface) const {
    state.components.clear();

    for (const auto& cmp : cmps) {
      using Opt = std::optional<BoundSymMatrix>;
      auto& tp = *cmp.trackStateProxy;
      auto track_pars = BoundTrackParameters(
          surface.getSharedPtr(), tp.filtered(),
          (state.covTransport ? Opt{BoundSymMatrix{tp.filteredCovariance()}}
                              : Opt{}));

      state.components.push_back(
          {SingleStepper::makeState(state.geoContext, state.magContext,
                                    std::move(track_pars), state.navDir),
           cmp.weight, Intersection3D::Status::onSurface});

      state.components.back().state.jacobian = cmp.jacobian;
      state.components.back().state.derivative = cmp.derivative;
      state.components.back().state.jacToGlobal = cmp.jacToGlobal;
      state.components.back().state.jacTransport = cmp.jacTransport;
    }
  }

  /// Method for on-demand transport of the covariance
  /// to a new curvilinear frame at current  position,
  /// or direction of the state
  ///
  /// @param [in,out] state State of the stepper
  ///
  /// @return the full transport jacobian
  void transportCovarianceToCurvilinear(State& state) const {
    for (auto& component : state.components) {
      SingleStepper::transportCovarianceToCurvilinear(component.state);
    }
  }

  /// Method for on-demand transport of the covariance
  /// to a new curvilinear frame at current position,
  /// or direction of the state
  ///
  /// @tparam surface_t the Surface type
  ///
  /// @param [in,out] state State of the stepper
  /// @param [in] surface is the surface to which the covariance is forwarded
  /// to
  /// @note no check is done if the position is actually on the surface
  void transportCovarianceToBound(State& state, const Surface& surface) const {
    for (auto& component : state.components) {
      SingleStepper::transportCovarianceToBound(component.state, surface);
    }
  }

  /// Perform a Runge-Kutta track parameter propagation step
  ///
  /// @param [in,out] state is the propagation state associated with the track
  /// parameters that are being propagated.
  ///
  /// The state contains the desired step size. It can be negative during
  /// backwards track propagation, and since we're using an adaptive
  /// algorithm, it can be modified by the stepper class during propagation.
  template <typename propagator_state_t>
  Result<double> step(propagator_state_t& state) const {
    State& stepping = state.stepping;

    // Emit warning if charge is not the same for all componenents
    {
      std::stringstream ss;
      bool charge_ambigous = false;
      for (const auto& cmp : stepping.components) {
        ss << cmp.state.q << " ";
        if (cmp.state.q != stepping.components.front().state.q) {
          charge_ambigous = true;
        }
      }

      if (charge_ambigous) {
        ACTS_VERBOSE(stepping.steps << "Charge of components is ambigous: "
                                    << ss.str());
      } else {
        ACTS_VERBOSE(stepping.steps << "Charge of components: " << ss.str());
      }
    }

    // Update step count
    stepping.steps++;

    // Check if we abort because of m_stepLimitAfterFirstComponentOnSurface
    if (stepping.stepCounterAfterFirstComponentOnSurface) {
      (*stepping.stepCounterAfterFirstComponentOnSurface)++;

      // If the limit is reached, remove all components which are not on a
      // surface, reweight the components, perform no step and return 0
      if (*stepping.stepCounterAfterFirstComponentOnSurface >=
          m_stepLimitAfterFirstComponentOnSurface) {
        auto& cmps = stepping.components;

        // It is not possible to remove components from the vector, since the
        // GSF actor relies on the fact that the ordering and number of
        // components does not change
        for (auto& cmp : cmps) {
          if (cmp.status != Intersection3D::Status::onSurface) {
            cmp.status = Intersection3D::Status::missed;
            cmp.weight = 0.0;
            cmp.state.pars.template segment<3>(eFreeDir0) = Vector3::Zero();
          }
        }

        // Reweight
        const auto sum_of_weights = std::accumulate(
            begin(cmps), end(cmps), ActsScalar{0},
            [](auto sum, const auto& cmp) { return sum + cmp.weight; });
        for (auto& cmp : cmps) {
          cmp.weight /= sum_of_weights;
        }

        ACTS_VERBOSE(
            "hit m_stepLimitAfterFirstComponentOnSurface, "
            "perform no step");

        stepping.stepCounterAfterFirstComponentOnSurface.reset();

        return 0.0;
      }
    }


    // Loop over all components and collect results in vector, write some
    // summary information to a stringstream
    std::vector<Result<double>> results;
    std::stringstream ss;

    for (auto& component : stepping.components) {
      // We must also propagate missed components for the case that all
      // components miss the target we need to retarget
      if (component.status == Intersection3D::Status::onSurface) {
        ss << "cmp skipped\t";
        continue;
      }

      SinglePropState single_state{component.state, state.navigation,
                                   state.options, state.geoContext};
      results.push_back(SingleStepper::step(single_state));

      if (results.back().ok()) {
        ss << *results.back() << "\t";
      } else {
        ss << "step error: " << results.back().error() << "\t";
      }
    }

    // Return no component was updated
    if (results.empty()) {
      return 0.0;
    }

    // Collect pointers to results which are ok, since Result is not copyable
    std::vector<Result<double>*> ok_results;
    for (auto& res : results) {
      if (res.ok()) {
        ok_results.push_back(&res);
      }
    }

    // Return error if there is no ok result
    if (ok_results.empty()) {
      return GsfError::AllComponentsSteppingError;
    }

    // Print the summary
    if (ok_results.size() == results.size()) {
      ACTS_VERBOSE("Performed steps: " << ss.str());
    } else {
      ACTS_WARNING("Performed steps with errors: " << ss.str());
    }

    // Compute the average stepsize for the return value and the
    // pathAccumulated
    const auto avg_step =
        std::accumulate(begin(ok_results), end(ok_results), 0.,
                        [](auto sum, auto res) { return sum + res->value(); }) /
        static_cast<double>(ok_results.size());
    stepping.pathAccumulated += avg_step;

    return avg_step;
  }
};

}  // namespace Acts
