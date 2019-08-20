// This file is part of the Acts project.
//
// Copyright (C) 2017-2018 CERN for the benefit of the Acts project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at http://mozilla.org/MPL/2.0/.

#include "Acts/Plugins/DD4hep/ActsExtension.hpp"

Acts::ActsExtension::ActsExtension(const std::string& axes) {
  addType("axes", "definitions", axes);
}

Acts::ActsExtension::ActsExtension(const ActsExtension& ext,
                                   const dd4hep::DetElement& /*elem*/)
    : m_flagStore(ext.m_flagStore), m_valueStore(ext.m_valueStore) {}

double Acts::ActsExtension::getValue(const std::string& tag,
                                     const std::string& category) const
    noexcept(false) {
  return getT(m_valueStore, tag, category);
}

void Acts::ActsExtension::addValue(double value, const std::string& tag,
                                   const std::string& category) {
  addT(m_valueStore, value, tag, category, 0.0);
}

bool Acts::ActsExtension::hasValue(const std::string& tag,
                                   const std::string& category) const {
  return hasT(m_valueStore, tag, category);
}

bool Acts::ActsExtension::hasType(const std::string& type,
                                  const std::string& category) const {
  return hasT(m_flagStore, type, category);
}

void Acts::ActsExtension::addType(const std::string& type,
                                  const std::string& category,
                                  const std::string& word) {
  std::string catDec = "<-- category -->";
  addT(m_flagStore, word, type, category, catDec);
}

const std::string Acts::ActsExtension::getType(
    const std::string& type, const std::string& category) const
    noexcept(false) {
  return getT(m_flagStore, type, category);
}

std::string Acts::ActsExtension::toString() const {
  std::string rString = "--------------- Acts::ActsExtension --------------- ";
  rString += '\n';
  rString += "- type store: ";
  rString += '\n';
  for (auto const& [key, value] : m_flagStore) {
    rString += key;
    rString += " : ";
    rString += value;
    rString += '\n';
  }
  return rString;
}
