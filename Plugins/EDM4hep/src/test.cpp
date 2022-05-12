#include "Acts/Plugins/EDM4hep/test.hpp"

#include "podio/EventStore.h"
#include "podio/ROOTReader.h"

namespace Acts {

void open(std::string path) {
  podio::ROOTReader reader;
  reader.openFile(path);

  auto store = podio::EventStore();
  store.setReader(&reader);

  std::cout << "events " << reader.getEntries() << std::endl;
}

}
