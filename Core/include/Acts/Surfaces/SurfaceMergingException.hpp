// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include <exception>
#include <memory>
#include <string>

namespace Acts {
class Surface;

/// @brief Exception type failures to merge two surfaces
class SurfaceMergingException : public std::exception {
 public:
  /// Constructor for surface merging exception
  /// @param surfaceA First surface that failed to merge
  /// @param surfaceB Second surface that failed to merge
  /// @param reason Description of why the merge failed
  SurfaceMergingException(std::weak_ptr<const Surface> surfaceA,
                          std::weak_ptr<const Surface> surfaceB,
                          const std::string& reason)
      : m_surfaceA(std::move(surfaceA)),
        m_surfaceB(std::move(surfaceB)),
        m_message(std::string{"Failure to merge surfaces: "} + reason) {}

  /// Get exception description
  /// @return C-style string describing the exception
  const char* what() const throw() override { return m_message.c_str(); }

 private:
  std::weak_ptr<const Surface> m_surfaceA;
  std::weak_ptr<const Surface> m_surfaceB;
  std::string m_message;
};

}  // namespace Acts
