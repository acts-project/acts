// This file is part of the Acts project.
//
// Copyright (C) 2016-2020 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/EventData/detail/PrintParameters.hpp"

#include "Acts/Surfaces/Surface.hpp"

#include <ostream>

void Acts::detail::printBoundParameters(std::ostream& os,
                                        const Acts::Surface& surface,
                                        const Acts::BoundVector& params,
                                        const Acts::BoundSymMatrix* cov) {
  // Set stream output format
  auto oldPrecision = os.precision(7);
  auto oldFlags = os.setf(std::ios::fixed);

  os << "BoundTrackParameters:\n";
  os << "  parameters: " << params.transpose() << '\n';
  if (cov) {
    os << "  covariance:\n";
    os << *cov << '\n';
  } else {
    os << "  no covariance stored\n";
  }
  os << "  on surface: " << surface.geometryId() << ' ' << surface.name()
     << '\n';

  // Reset stream format
  os.setf(oldFlags);
  os.precision(oldPrecision);
}

void Acts::detail::printFreeParameters(std::ostream& os,
                                       const Acts::FreeVector& params,
                                       const Acts::FreeMatrix* cov) {
  // Set stream output format
  auto oldPrecision = os.precision(7);
  auto oldFlags = os.setf(std::ios::fixed);

  os << "FreeTrackParameters:\n";
  os << "  parameters: " << params.transpose() << '\n';
  if (cov) {
    os << "  covariance:\n";
    os << *cov << '\n';
  } else {
    os << "  no covariance stored\n";
  }

  // Reset stream format
  os.setf(oldFlags);
  os.precision(oldPrecision);
}
