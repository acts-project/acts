#include "Acts/Plugins/EDM4hep/test.hpp"

#include "edm4hep/SimCalorimeterHitCollection.h"
#include "edm4hep/SimTrackerHitCollection.h"
#include "podio/EventStore.h"
#include "podio/ROOTReader.h"

namespace Acts {

void open(std::string path) {
  podio::ROOTReader reader;
  reader.openFile(path);

  auto store = podio::EventStore();
  store.setReader(&reader);

  std::cout << "events " << reader.getEntries() << std::endl;

  reader.getCollectionIDTable()->print();

  for (const auto& name : reader.getCollectionIDTable()->names()) {
    auto& sths = store.get<edm4hep::SimTrackerHitCollection>(name);

    if (sths.isValid()) {
      std::cout << "tacker hits " << sths.size() << std::endl;
    }
  }
}

}  // namespace Acts
