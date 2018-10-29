// This file is part of the Acts project.
//
// Copyright (C) 2016-2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Propagator/detail/extension_list_implementation.hpp"
#include "Acts/Utilities/detail/Extendable.hpp"
#include "Acts/Utilities/detail/MPL/all_of.hpp"
#include "Acts/Utilities/detail/MPL/has_duplicates.hpp"

namespace Acts {

/// The ExtensionList allows to add an arbitrary number of step evaluation
/// algorithms for the RKN4 evaluation. These can be categorised in two general
/// types:
/// a) Step evaluation that should not be evaluated along other extensions, or
/// at least would overwrite partial results. This means that in the best case
/// unnecessary/redundant calculation would be performed, in the worst case the
/// evaluation would go wrong.
/// b) The step evaluation remains untouched and only further calculations are
/// performed (like additional features or data gathering) that can be treated
/// as independent of the basic step evaluation in type a). These types can be
/// added but do not require special treatment in order to not ruin the step
/// evaluation.
/// The concept of the auctioneers aims in the first place to judge which
/// extension of category a) is the one to go. Although every extension can
/// judge if it is valid based on the data given from the state of stepper,
/// multiple extensions from type a) could fulfill their dependencies. Since an
/// extension does not know about other extensions, the decision for the best
/// extension for the step can only be estimated on a global scope. This is the
/// job of the auctioneers.
/// TODO: An anticipation of an optimal concept of the input (and maybe alos the
/// output) of the call operator of an auctioneer cannot be performed at the
/// current stage but the concept of passing booblean vectors could be extended
/// to vectors of ints or doubles (this is the equivalent to a bid which every
/// extension can make for the upcoming step). At the current stage, a
/// bid-system would be pure guessing.

/// @brief Auctioneer that takes all extensions as valid that state to be valid
struct VoidAuctioneer
{
  /// @brief Default constructor
  VoidAuctioneer() = default;

  /// @brief Call operator that just returns the list of candidates as valids
  ///
  /// @param [in] vCandidates Candidates that are treated as valid extensions
  /// @return The to vCandidates identical list of valid extensions
  std::vector<bool>
  operator()(std::vector<bool> vCandidates)
  {
    return std::move(vCandidates);
  }
};

/// @brief Auctioneer that states only the first valid extensions as indeed
/// valid extension
struct FirstValidAuctioneer
{
  /// @brief Default constructor
  FirstValidAuctioneer() = default;

  /// @brief Call operator that states the first valid extension as the only
  /// valid extension
  ///
  /// @param [in] vCandidates Candidates for a valid extension
  /// @return List with at most one valid extension
  std::vector<bool>
  operator()(std::vector<bool> vCandidates)
  {
    // Indicator if the first valid was already found
    bool firstValidFound = false;
    //~ for (bool& v : vCandidates) {
    for (unsigned int i = 0; i < vCandidates.size(); i++) {
      // If a valid extensions is already found, set all following to false
      if (firstValidFound) vCandidates[i] = false;
      // If the first valid isn't found yet, toggle the flag on the first found
      else if (vCandidates[i])
        firstValidFound = true;
    }
    return std::move(vCandidates);
  }
};

/// @brief Container of extensions used in the stepper of the propagation
/// @tparam extensions Types of the extensions
template <typename... extensions>
struct ExtensionList : private detail::Extendable<extensions...>
{
private:
  // Checkout for duplicates in the extensions
  static_assert(not detail::has_duplicates_v<extensions...>,
                "same extension type specified several times");

  // Access to all extensions
  using detail::Extendable<extensions...>::tuple;

  using impl = detail::extension_list_impl<extensions...>;

  // Vector of valid extensions for a step
  std::vector<bool> validExtensions;

  /// @brief Evaluation function to set valid extensions for an upcoming
  /// integration step
  ///
  /// @tparam stepper_state_t Type of the state of the stepper
  /// @tparam auctioneer_t Type of the auctioneer for post-processing of valid
  /// extension candidates
  /// @param [in] state State of the stepper
  template <typename stepper_state_t, typename auctioneer_t = VoidAuctioneer>
  void
  validExtensionForStep(const stepper_state_t& state)
  {
    auctioneer_t      auctioneer;
    std::vector<bool> validExtensionCandidates;
    // Ask all extensions for a boolean statement of their validity
    impl::validExtensionForStep(tuple(), state, validExtensionCandidates);
    // Post-process the vector in an auctioneer
    validExtensions = auctioneer(std::move(validExtensionCandidates));
  }

public:
  // Access to an extension
  using detail::Extendable<extensions...>::get;

  /// @brief This functions implies the call of the method k(). It collects all
  /// extensions and arguments and passes them forward for evaluation and
  /// returns a boolean.
  template <typename stepper_state_t, typename auctioneer_t = VoidAuctioneer>
  bool
  k(const stepper_state_t& state,
    Vector3D&              knew,
    const Vector3D&        bField,
    const int              i     = 0,
    const double           h     = 0,
    const Vector3D&        kprev = Vector3D())
  {
    if (i == 0) validExtensionForStep<stepper_state_t, auctioneer_t>(state);
    return impl::k(
        tuple(), state, knew, bField, validExtensions.cbegin(), i, h, kprev);
  }

  /// @brief This functions implies the call of the method finalize(). It
  /// collects all extensions and arguments and passes them forward for
  /// evaluation and returns a boolean.
  template <typename stepper_state_t, typename stepper_data_t>
  bool
  finalize(stepper_state_t&      state,
           const double          h,
           const stepper_data_t& data,
           ActsMatrixD<7, 7>& D)
  {
    return impl::finalize(tuple(), state, h, data, D, validExtensions.cbegin());
  }

  /// @brief This functions implies the call of the method finalize(). It
  /// collects all extensions and arguments and passes them forward for
  /// evaluation and returns a boolean.
  template <typename stepper_state_t>
  bool
  finalize(stepper_state_t& state, const double h)
  {
    return impl::finalize(tuple(), state, h, validExtensions.cbegin());
  }
};

}  // namespace Acts
