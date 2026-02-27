// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Io/Podio/PodioCollectionDataHandle.hpp"

#include "ActsExamples/Framework/WhiteBoard.hpp"

#include <memory>

#include <podio/CollectionBase.h>

namespace ActsExamples {

CollectionBaseWriteHandle::CollectionBaseWriteHandle(SequenceElement* parent,
                                                     const std::string& name)
    : WriteDataHandleBase{parent, name} {
  registerAsWriteHandle();
}

void CollectionBaseWriteHandle::store(
    WhiteBoard& wb, std::unique_ptr<podio::CollectionBase> collection) const {
  if (!isInitialized()) {
    throw std::runtime_error{"CollectionBaseHandle '" + fullName() +
                             "' not initialized"};
  }

  add(wb, std::move(collection));
}

bool CollectionBaseWriteHandle::isCompatible(
    const DataHandleBase& other) const {
  const auto baseCollectionHash =
      Acts::typeHash<std::unique_ptr<podio::CollectionBase>>();
  if (other.typeHash() == baseCollectionHash) {
    return true;
  }

  const auto* typedReadHandle =
      dynamic_cast<const detail::PodioCollectionTypedReadHandle*>(&other);
  return typedReadHandle != nullptr &&
         typedReadHandle->collectionTypeHash() == m_declaredTypeHash;
}

const std::type_info& CollectionBaseWriteHandle::typeInfo() const {
  return *m_declaredTypeInfo;
}

std::uint64_t CollectionBaseWriteHandle::typeHash() const {
  return m_declaredTypeHash;
}

}  // namespace ActsExamples
