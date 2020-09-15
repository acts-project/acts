// This file is part of the Acts project.
//
// Copyright (C) 2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Surfaces/Surface.hpp"
#include "ActsFatras/EventData/Hit.hpp"

namespace ActsFatras {

/// A digitizeed hit on a surface.
///
/// This is the digitized hit, output of any digitization module. It serves
/// as source_link in the Surface based Measurement classes and thus has to
/// fullfill the SourceLinkConcept.
///
class DigitizedHit {
 public:
  /// @brief Nested IContent base class
  ///
  /// This will be overwritten by specialisations of the
  /// varius digitisers that will know their own content
  class IContent {
   public:
    virtual ~IContent() = default;
  };

  DigitizedHit() = default;

  /// Constructor with arguments
  ///
  /// @param hits The simulated hits at input
  /// @param surface The surface for the measuremetn representation
  /// @param content The digitization content to be filled by the digitizer
  DigitizedHit(std::vector<Hit>&& hits, const Acts::Surface& surface,
               std::unique_ptr<const DigitizedHit::IContent> content);

  virtual ~DigitizedHit() = default;

  /// Required copy constructor
  ///
  /// @param other The source
  DigitizedHit(const DigitizedHit& other);

  /// Required equality operator
  ///
  /// @param other The one to compare to
  bool operator==(const DigitizedHit& other) const;

  /// Return the reference Surface
  const Acts::Surface& referenceSurface() const;

  /// Return the simulated hits
  const std::vector<Hit>& simulatedHits() const;

  /// Return the Digitization content
  const IContent& content() const;

 private:
  std::vector<Hit> m_simulatedHits = {};
  std::shared_ptr<const Acts::Surface> m_surface = nullptr;
  std::unique_ptr<const IContent> m_content = nullptr;
};

inline const std::vector<Hit>& DigitizedHit::simulatedHits() const {
  return m_simulatedHits;
}

inline const Acts::Surface& DigitizedHit::referenceSurface() const {
  return (*m_surface.get());
}

inline const DigitizedHit::IContent& DigitizedHit::content() const {
  return (*m_content.get());
}

}  // namespace ActsFatras