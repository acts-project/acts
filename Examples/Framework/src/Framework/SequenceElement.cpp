// This file is part of the Acts project.
//
// Copyright (C) 2017-2023 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "ActsExamples/Framework/SequenceElement.hpp"

#include "ActsExamples/Framework/DataHandle.hpp"

namespace ActsExamples {

void SequenceElement::registerWriteHandle(const DataHandleBase& handle) {
  m_writeHandles.push_back(&handle);
}

void SequenceElement::registerReadHandle(const DataHandleBase& handle) {
  m_readHandles.push_back(&handle);
}

const std::vector<const DataHandleBase*>& SequenceElement::writeHandles()
    const {
  return m_writeHandles;
}

const std::vector<const DataHandleBase*>& SequenceElement::readHandles() const {
  return m_readHandles;
}

}  // namespace ActsExamples
