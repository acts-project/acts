// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTS/Plugins/DD4hepPlugins/ActsExtension.hpp"

Acts::ActsExtension::ActsExtension(const Config& cfg) : Acts::IActsExtension()
{
  setConfiguration(cfg);
}

Acts::ActsExtension::ActsExtension(const ActsExtension& det,
                                   const DD4hep::Geometry::DetElement&)
  : Acts::IActsExtension(), m_cfg(det.m_cfg)
{
}

void
Acts::ActsExtension::setConfiguration(const Acts::ActsExtension::Config& config)
{
  // @todo check consistency
  // copy the configuration
  m_cfg = config;
}
