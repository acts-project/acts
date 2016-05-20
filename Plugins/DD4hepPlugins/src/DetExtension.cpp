// This file is part of the ACTS project.
//
// Copyright (C) 2016 ACTS project team
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ACTS/Plugins/DD4hepPlugins/DetExtension.hpp"

Acts::DetExtension::DetExtension() :
Acts::IDetExtension(),
m_segmentation(nullptr),
m_shape(Acts::ShapeType::None),
m_supportStructure(),
m_modules()
{}

Acts::DetExtension::~DetExtension()
{}