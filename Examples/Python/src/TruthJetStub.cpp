// This is a stub version when the FastJet plugin is not available
#include <pybind11/pybind11.h>

namespace py = pybind11;

namespace Acts::Examples::Python {
void exportTruthJet(py::module& m) {
  // Empty implementation when FastJet is not available
}
}