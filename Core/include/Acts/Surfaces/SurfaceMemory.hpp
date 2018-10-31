// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SurfaceMemory.h, Acts project
///////////////////////////////////////////////////////////////////

#pragma once

namespace Acts {

namespace {

  /// @brief Clone helpers to deal correclty with
  /// Surface pointers in terms of clone, delete and move operations
  ///
  /// @tparam surface_t Type of the surface object
  ///
  /// @param source is the source surface
  template <typename surface_t>
  const surface_t*
  createPtr(const surface_t& source)
  {
    // Clone if the surface is free
    if (source.isFree()) {
      const surface_t* cSurface = source.clone();
      return cSurface;
    }
    return (&source);
  }

  /// @brief Delete helper to deal correctly with surface pointers
  ///
  /// @tparam surface_t Type of the surface object
  ///
  /// @param surface is the surface to be deleted, it can be nullptr
  /// as part of a move operation with running out of scope
  template <typename surface_t>
  void
  deletePtr(const surface_t*& surface)
  {
    if (surface and surface->isFree()) {
      delete (surface);
    }
    // Set the target surface to null ptr
    surface = nullptr;
  }

  /// @brief Move helper to deal correctly with surface pointers
  ///
  /// @tparam surface_t Type of the surface object
  ///
  /// @param surface is the surface to be move
  template <typename surface_t>
  const surface_t*
  movePtr(const surface_t*& rhs)
  {
    // @TODO FIX THIS
    /// still a problem with BOOST variant
    return createPtr(*rhs);
    // const surface_t* returnSf = rhs;
    // rhs = nullptr; // The move has to set it to zero
    // return returnSf;
  }
}
};
