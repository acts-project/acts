// This file is part of the Acts project.
//
// Copyright (C) 2018 Acts project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

///////////////////////////////////////////////////////////////////
// SurfaceMaterialMapper.cpp, Acts project
///////////////////////////////////////////////////////////////////

#include "Acts/Plugins/MaterialMapping/SurfaceMaterialMapper.hpp"

Acts::SurfaceMaterialMapper::SurfaceMaterialMapper(
    const Config&                 cfg,
    StraightLinePropagator        propagator,
    std::unique_ptr<const Logger> slogger)
  : m_cfg(cfg)
  , m_propagator(std::move(propagator))
  , m_logger(std::move(slogger))
{
}
