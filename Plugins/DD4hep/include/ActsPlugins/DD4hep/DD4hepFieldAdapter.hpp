// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/MagneticField/MagneticFieldProvider.hpp"

#include <memory>

namespace dd4hep {
class OverlayedField;
}

namespace ActsPlugins {
/// @addtogroup dd4hep_plugin
/// @{

/// @ingroup magnetic_field dd4hep_plugin
class DD4hepFieldAdapter : public Acts::MagneticFieldProvider {
  /// Cache object for DD4hep field adapter
  /// @note As DD4hep does not implement a caching mechanism, this struct is
  ///       empty.
  struct Cache {};

 public:
  explicit DD4hepFieldAdapter(dd4hep::OverlayedField field);

  Acts::MagneticFieldProvider::Cache makeCache(
      const Acts::MagneticFieldContext& mctx) const override;

  Acts::Result<Acts::Vector3> getField(
      const Acts::Vector3& position,
      Acts::MagneticFieldProvider::Cache& cache) const override;

 private:
  double m_fieldConversionFactor;
  double m_lengthConversionFactor;
  std::unique_ptr<dd4hep::OverlayedField> m_field;
};

/// @}
}  // namespace ActsPlugins
