// This file is part of the Acts project.
//
// Copyright (C) 2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Detector/detail/GridAxisGenerators.hpp"
#include "Acts/Detector/detail/IndexedSurfacesGenerator.hpp"
#include "Acts/Geometry/Extent.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Navigation/NavigationDelegates.hpp"
#include "Acts/Navigation/SurfaceCandidatesUpdators.hpp"
#include "Acts/Plugins/Json/AlgebraJsonConverter.hpp"
#include "Acts/Plugins/Json/GridJsonConverter.hpp"
#include "Acts/Utilities/Enumerate.hpp"
#include "Acts/Utilities/Grid.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "Acts/Utilities/detail/AxisFwd.hpp"

#include <tuple>
#include <vector>

namespace Acts {

using namespace Experimental::detail::GridAxisGenerators;

namespace IndexedSurfacesJsonConverter {

// Generate the possible axes in this case
static auto s_possibleAxes =
    std::tuple<EqBound, EqOpen, EqClosed,
               // All 1D Var  options
               VarBound, VarOpen, VarClosed,
               // All 2D EqEq options
               EqBoundEqBound, EqBoundEqOpen, EqBoundEqClosed, EqOpenEqBound,
               EqOpenEqOpen, EqOpenEqClosed, EqClosedEqBound, EqClosedEqOpen,
               EqClosedEqClosed,
               // All 2D EqVar options
               EqBoundVarBound, EqBoundVarOpen, EqBoundVarClosed,
               EqOpenVarBound, EqOpenVarOpen, EqOpenVarClosed, EqClosedVarBound,
               EqClosedVarOpen, EqClosedVarClosed,
               // All 2D VarEq options
               VarBoundEqBound, VarBoundEqOpen, VarBoundEqClosed,
               VarOpenEqBound, VarOpenEqOpen, VarOpenEqClosed, VarClosedEqBound,
               VarClosedEqOpen, VarClosedEqClosed,
               // All 2D VarEq options
               VarBoundVarBound, VarBoundVarOpen, VarBoundVarClosed,
               VarOpenVarBound, VarOpenVarOpen, VarOpenVarClosed,
               VarClosedVarBound, VarClosedVarOpen, VarClosedVarClosed>{};

/// @brief The actual conversion method
template <typename index_grid>
nlohmann::json convertImpl(const index_grid& indexGrid) {
  nlohmann::json jIndexedSurfaces;

  // Fill the casts
  nlohmann::json jCasts;
  // 1D casts
  if constexpr (index_grid::grid_type::DIM == 1u) {
    jCasts.push_back(indexGrid.casts[0u]);
  }
  // 1D casts
  if constexpr (index_grid::grid_type::DIM == 2u) {
    jCasts.push_back(indexGrid.casts[0u]);
    jCasts.push_back(indexGrid.casts[1u]);
  }
  jIndexedSurfaces["casts"] = jCasts;
  jIndexedSurfaces["transform"] =
      Transform3JsonConverter::toJson(indexGrid.transform);
  jIndexedSurfaces["grid"] = GridJsonConverter::toJson(indexGrid.grid);
  return jIndexedSurfaces;
}

/// @brief Convert the single delegate if it is of the type of the reference
///
/// @note It will do nothing if the type does not match
///
/// @param jIndexedSurfaces the json object to be filled
/// @param delegate the delegate to be translated
/// @param refInstance is a reference instance of potential type casting
template <typename instance_type>
void convert(nlohmann::json& jIndexedSurfaces,
             const Experimental::SurfaceCandidatesUpdator& delegate,
             [[maybe_unused]] const instance_type& refInstance) {
  using GridType =
      typename instance_type::template grid_type<std::vector<size_t>>;
  // Defining a Delegate type
  using DelegateType = Experimental::IndexedSurfacesAllPortalsImpl<
      GridType, Experimental::IndexedSurfacesImpl>;
  using SubDelegateType = Experimental::IndexedSurfacesImpl<GridType>;

  // Get the instance
  const auto* instance = delegate.instance();
  auto castedDelegate = dynamic_cast<const DelegateType*>(instance);
  if (castedDelegate != nullptr) {
    // Get the surface updator
    auto indexedSurfaces = std::get<SubDelegateType>(castedDelegate->updators);
    jIndexedSurfaces = convertImpl<SubDelegateType>(indexedSurfaces);
  }
}

/// @brief Unrolling function for catching the right instance
///
/// @param jIndexedSurfaces the json object to be filled
/// @param delegate the delegate to be translated
/// @param axesTuple the tuple of axes to be unrolled
///
/// @note parameters are as of the `convertImpl` method
template <typename tuple_type, size_t... I>
void unrollConvert(nlohmann::json& jIndexedSurfaces,
                   const Experimental::SurfaceCandidatesUpdator& delegate,
                   const tuple_type& axesTuple,
                   std::index_sequence<I...> /*unused*/) {
  (convert(jIndexedSurfaces, delegate, std::get<I>(axesTuple)), ...);
}

/// Convert a surface array into needed constituents
///
/// @param delegate the delegate to be translated
///
/// @note this is the entry point of the conversion, i.e. top of the
/// unrolling loop
///
/// @return a collection of proto surface object and a grid, and associations
static inline nlohmann::json toJson(
    const Experimental::SurfaceCandidatesUpdator& delegate) {
  // Convert if dynamic cast happens to work
  nlohmann::json jIndexedSurfaces;
  unrollConvert(jIndexedSurfaces, delegate, s_possibleAxes,
                std::make_index_sequence<
                    std::tuple_size<decltype(s_possibleAxes)>::value>());
  // Return the newly filled ones
  jIndexedSurfaces["type"] = "IndexedSurfaces";
  return jIndexedSurfaces;
}

/// @brief Convert the single delegate if it is of the type of the reference
///
/// @param jSurfaceNavigation the json file to read from
///
/// @return the surface navigation delegate
Experimental::SurfaceCandidatesUpdator fromJson(
    const nlohmann::json& jSurfaceNavigation);

}  // namespace IndexedSurfacesJsonConverter
}  // namespace Acts
