// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.

#include "ActsExamples/Framework/DataHandle.hpp"

#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/Framework/Sequencer.hpp"

#include <regex>

#include <boost/algorithm/string.hpp>
#include <boost/core/demangle.hpp>

namespace ActsExamples {

namespace {
/// Shorten some common but lengthy C++ constructs
std::string demangleAndShorten(std::string name) {
  name = boost::core::demangle(name.c_str());

  // Remove std::allocator from vector
  const static std::regex vector_pattern(
      R"??(std::vector<(.*), std::allocator<(\1\s*)>\s*>)??");
  name = std::regex_replace(name, vector_pattern, "std::vector<$1>");

  // Shorten Acts::BoundVariantMeasurement
  const static std::regex variant_pattern(
      R"??(std::variant<(Acts::Measurement<Acts::BoundIndices, [0-9]ul>(,|)\s+)+>)??");
  name = std::regex_replace(name, variant_pattern,
                            "Acts::BoundVariantMeasurement");

  // strip namespaces
  boost::algorithm::replace_all(name, "std::", "");
  boost::algorithm::replace_all(name, "boost::container::", "");
  boost::algorithm::replace_all(name, "Acts::", "");
  boost::algorithm::replace_all(name, "ActsExamples::", "");
  boost::algorithm::replace_all(name, "ActsFatras::", "");

  return name;
}

std::string symbol(const char* in) {
  std::stringstream ss;
  std::string s = demangleAndShorten(in);
  std::size_t pos = 0;
  while (pos + 80 < s.size()) {
    ss << "   " + s.substr(pos, 80);
    pos += 80;
  }
  ss << s.substr(pos);
  return ss.str();
};

}  // namespace

void WriteDataHandleBase::initialize(std::string_view key) {
  if (key.empty()) {
    throw std::invalid_argument{"Write handle '" + fullName() +
                                "' cannot receive empty key"};
  }
  m_key = key;
}

void DataHandleBase::maybeInitialize(std::string_view key) {
  if (!key.empty()) {
    m_key = key;
  }
}

bool WriteDataHandleBase::isCompatible(const DataHandleBase& other) const {
  return dynamic_cast<const ReadDataHandleBase*>(&other) != nullptr &&
         typeInfo() == other.typeInfo();
}

void WriteDataHandleBase::emulate(StateMapType& state,
                                  WhiteBoard::AliasMapType& aliases,
                                  const Acts::Logger& logger) const {
  if (!isInitialized()) {
    return;
  }

  ACTS_INFO("-> " << name() << " '" << key() << "':");

  if (auto it = state.find(key()); it != state.end()) {
    const auto& source = *it->second;
    ACTS_ERROR("White board will already contain key '"
               << key() << "'. Source: '" << source.fullName()
               << "' (cannot overwrite)");
    throw SequenceConfigurationException{
        "Incompatible data handle type for key '" + key() +
        "': " + source.fullName()};
  }

  state.try_emplace(key(), this);

  if (auto it = aliases.find(key()); it != aliases.end()) {
    ACTS_DEBUG("Key '" << key() << "' aliased to '" << it->second << "'");
    state.try_emplace(it->second, this);
  }
}

void ReadDataHandleBase::initialize(std::string_view key) {
  if (key.empty()) {
    throw std::invalid_argument{"Read handle '" + fullName() +
                                "' cannot receive empty key"};
  }
  m_key = key;
}

bool ReadDataHandleBase::isCompatible(const DataHandleBase& other) const {
  return dynamic_cast<const WriteDataHandleBase*>(&other) != nullptr &&
         typeInfo() == other.typeInfo();
}

void ReadDataHandleBase::emulate(StateMapType& state,
                                 WhiteBoard::AliasMapType& /*aliases*/,
                                 const Acts::Logger& logger) const {
  if (!isInitialized()) {
    return;
  }

  ACTS_INFO("<- " << name() << " '" << key() << "':");
  ACTS_INFO("   " << symbol(typeInfo().name()));

  if (auto it = state.find(key()); it != state.end()) {
    const auto& source = *it->second;
    if (!source.isCompatible(*this)) {
      ACTS_ERROR(
          "Adding " << m_parent->typeName() << " " << m_parent->name() << ":"
                    << "\n-> white board will contain key '" << key() << "'"
                    << "\nat this point in the sequence (source: "
                    << source.fullName() << "),"
                    << "\nbut the type will be\n"
                    << "'" << demangleAndShorten(source.typeInfo().name())
                    << "'"
                    << "\nand not\n"
                    << "'" << demangleAndShorten(typeInfo().name()) << "'");
      throw SequenceConfigurationException{
          "Incompatible data handle type for key '" + key() +
          "': " + source.fullName()};
    }
  } else {
    ACTS_ERROR(
        "Adding " << m_parent->typeName() << " " << m_parent->name() << ":"
                  << "\n-> white board will not contain key"
                  << " '" << key() << "' at this point in the sequence."
                  << "\n   Needed for read data handle '" << name() << "'");
    throw SequenceConfigurationException{"Missing data handle for key '" +
                                         key() + "'"};
  }
}

void ConsumeDataHandleBase::emulate(StateMapType& state,
                                    WhiteBoard::AliasMapType& aliases,
                                    const Acts::Logger& logger) const {
  if (!isInitialized()) {
    return;
  }

  ACTS_INFO("<< " << name() << " '" << key() << "':");
  ACTS_INFO("   " << symbol(typeInfo().name()));

  if (auto it = state.find(key()); it != state.end()) {
    if (const auto& source = *it->second; !source.isCompatible(*this)) {
      ACTS_ERROR(
          "Adding " << m_parent->typeName() << " " << m_parent->name() << ":"
                    << "\n-> white board will contain key '" << key() << "'"
                    << "\nat this point in the sequence (source: "
                    << source.fullName() << "),"
                    << "\nbut the type will be\n"
                    << "'" << demangleAndShorten(source.typeInfo().name())
                    << "'"
                    << "\nand not\n"
                    << "'" << demangleAndShorten(typeInfo().name()) << "'");
      throw SequenceConfigurationException{
          "Incompatible data handle type for key '" + key() +
          "': " + source.fullName()};
    }
    // Remove the key from state since it will be consumed
    ACTS_VERBOSE("Removing key '" << key() << "' from state");
    state.erase(it);

    for (const auto& [source, target] : aliases) {
      if (target == key()) {
        ACTS_VERBOSE("Removing alias target '" << target << "' from state");
        state.erase(source);
      }
      if (source == key()) {
        ACTS_VERBOSE("Removing alias source '" << source << "' from state");
        state.erase(target);
      }
    }
  } else {
    ACTS_ERROR(
        "Adding " << m_parent->typeName() << " " << m_parent->name() << ":"
                  << "\n-> white board will not contain key"
                  << " '" << key() << "' at this point in the sequence."
                  << "\n   Needed for consume data handle '" << name() << "'");
    throw SequenceConfigurationException{"Missing data handle for key '" +
                                         key() + "'"};
  }
}

void DataHandleBase::registerAsWriteHandle() {
  m_parent->registerWriteHandle(*this);
}

void DataHandleBase::registerAsReadHandle() {
  m_parent->registerReadHandle(*this);
}

}  // namespace ActsExamples
