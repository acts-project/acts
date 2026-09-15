// This file is part of the ACTS project.
//
// Copyright (C) 2016 CERN for the benefit of the ACTS project
//
// This Source Code Form is subject to the terms of the Mozilla Public
// License, v. 2.0. If a copy of the MPL was not distributed with this
// file, You can obtain one at https://mozilla.org/MPL/2.0/.



#include "Acts/Utilities/Logger.hpp"

#include "GeoModelHelpers/MaterialManager.h"

namespace ActsPlugins {

/// GeoModelMaterialManager implementation reading a deserialized JSON file as
/// database. The dictionary holds the table of the chemical elements and of the
/// macroscopic materials. The first one is encoded in the "periodicTable" block
/// which is a list of dictionaries
///
/// {
///         "periodicTable" :  [
///                 { "name" : "Hydrogen", "symbol" : "H", "charge" : 1,
///                 "atomicNumber" : 1.2 },
///                 ....
///         ]
///
///    where 'charge' is the number of protons in the nuclei, and the
///    "atomicNumber" the average number of nucleons. The 'name' field
///    identifies the element uniquely in the database for later usage.
///
///  The actual materials used for the GeoModel tree are defined in the section
///  "materials"
///
///          "materials": [
///                  {"name" : "air" ,  "rho" : 0.001290,
///                   "composition": [
///                                    {"component" : "Nitrogen", "parts" :
///                                    0.7494},
///                                    {"component" : "Oxygen", "parts" :
///                                    0.2369},
///                                    {"component" : "Argon", "parts" :
///                                    0.0129},
///                                    {"component" : "Hydrogen", "parts" :
///                                    0.0008}
///                                   ]},
///                 {"name" : "CO2", "rho": 0.00184,
///                   "composition" : [{"component" : "Carbon",  "parts": 1.},
///                                    {"component"  : "Oxygen", "parts": 2.}
///                                    ]},
///                 {"name" : "ArCO2", "rho" : 0.0054, "namespace" : "Muon",
///                     "composition" : [{"component" : "Argon", "parts" :
///                     0.93},
///                                      {"component" : "std::CO2", 0.7}]}
///          ]
///     Each entry defines a material which can then be used to define the
///     GeoLogVol "name": is the Identifier of the material which is combined
///     with the current namespace state. By default the namespace is always
///     "std" but can be changed to anything using the "namespace" key in the
///     dictionary. The Material can then later be retrieved from the Manager
///     via
///                 auto arCo2 =
///                 MaterialManager::getManager()->getMaterial("Muon::ArCO2");
///
///     "rho" defines the material's density and "composition" lists all the
///     chemical elements
///      or other composed materials out of which the material is made of. Each
///      "component" is associated with "parts" which defines the relative
///      fractional component. Note: if the parts violate unity, the Manager
///      automatically re-normalizes the parts
///
///
class GeoModelJsonMaterialManager : public MaterialManager {
 public:
  /// Use the singleton access pattern to instantiate the
  /// GeoModelJsonMaterialManager the inherited class shares the same singleton
  /// pointer as the parent. In case of an uninstantiated singleton, the
  /// GeoModelJsonMaterialManager is instantiated. The manager pointer is
  /// dynamically casted and then returned.
  /// @return Return the pointer to the `GeoModelJsonMaterialManager` instance
  /// which may be a nullptr if another implementation of the MaterialManager
  /// was called before
  static GeoModelJsonMaterialManager* getManager();
  /// Load the database from a material JSON file
  /// @param filePath: Location of the JSON data file
  /// @return Whether the file exists and the parsing of the material DB was successful
  bool loadMaterialMap(const std::string& filePath);

 protected:
  GeoModelJsonMaterialManager() = default;

 private:
  std::unique_ptr<const Acts::Logger> m_logger{Acts::getDefaultLogger(
      "GeoModelMaterialMgr", Acts::Logging::Level::INFO)};
  /// @brief The Acts logger object
  /// @return The Acts logger object associated with this class
  const Acts::Logger& logger() const { return *m_logger; }
};
}  // namespace ActsPlugins
