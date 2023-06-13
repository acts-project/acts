#pragma once

#include <memory>
#include <string>

// DD4hep
#include "Acts/Digitization/CartesianSegmentation.hpp"
#include "Acts/Digitization/DigitizationModule.hpp"

#include "DD4hep/DetFactoryHelper.h"
#include "DD4hep/Detector.h"
#include "DD4hep/Objects.h"
#include "DD4hep/Segmentations.h"
#include "DDSegmentation/BitField64.h"
#include "DDSegmentation/CartesianGridXY.h"
#include "DDSegmentation/CartesianGridXYZ.h"
#include "DDSegmentation/PolarGridRPhi.h"
#include "XML/XMLDetector.h"
#include "XML/XMLElements.h"

namespace Acts {
class DigitizationModule;
}  // namespace Acts

/** Given a xml element with several daughters with the same name, e.g.
 <detector> <layer name="1" /> <layer name="2"> </detector>
 this method returns the first daughter of type nodeName whose attribute has a
 given value
 e.g. returns <layer name="2"/> when called with (detector, "layer", "name",
 "1") */
namespace det {
namespace utils {

/// Global method to build an Acts::DigitizationModule with rectangular
/// segmentation.
/// @note This function should be used in order to create the input
/// needed for construction with
/// Acts::ActsExtension(std::shared_ptr<const DigitizationModule>)
/// @param halflengthX The half length in x of the detector module
/// @param halflengthZ The half length in z of the detector module
/// @param thickness The thickness of the detector module
/// @param segmentation the DD4hep segmentation
std::shared_ptr<const Acts::DigitizationModule> rectangleDigiModuleXZ(
    double halflengthX, double halflengthZ, double thickness,
    const dd4hep::Segmentation& segmentation);

/// Global method to build an Acts::DigitizationModule with rectangular
/// segmentation.
/// @note This function should be used in order to create the input
/// needed for construction with
/// Acts::ActsExtension(std::shared_ptr<const DigitizationModule>)
/// @param halflengthX The half length in x of the detector module
/// @param halflengthZ The half length in z of the detector module
/// @param thickness The thickness of the detector module
/// @param gridSizeX The grid size in x
/// @param gridSizeY The grid size in y
std::shared_ptr<const Acts::DigitizationModule> rectangleDigiModuleXZ(
    double halflengthX, double halflengthZ, double thickness, double gridSizeX,
    double gridSizeZ);

/// Global method to build an Acts::DigitizationModule with trapezoidal
/// segmentation.
/// @note This function should be used in order to create the input
/// needed for construction with
/// Acts::ActsExtension(std::shared_ptr<const DigitizationModule>)
/// @param minHalflengthX The half length in x of the detector module on the
/// negative side of z
/// @param maxHalflengthX The half length in x of the detector module on the
/// positive side of z
/// @param halflengthZ The half length in z of the detector module
/// @param thickness The thickness of the detector module
/// @param segmentation the DD4hep segmentation
std::shared_ptr<const Acts::DigitizationModule> trapezoidalDigiModuleXZ(
    double minHalflengthX, double maxHalflengthX, double halflengthZ,
    double thickness, const dd4hep::Segmentation& segmentation);

/// Global method to build an Acts::DigitizationModule with trapezoidal
/// segmentation.
/// @note This function should be used in order to create the input
/// needed for construction with
/// Acts::ActsExtension(std::shared_ptr<const DigitizationModule>)
/// @param minHalflengthX The half length in x of the detector module on the
/// negative side of z
/// @param maxHalflengthX The half length in x of the detector module on the
/// positive side of z
/// @param halflengthZ The half length in z of the detector module
/// @param thickness The thickness of the detector module
/// @param gridSizeX The grid size in x
/// @param gridSizeY The grid size in y
std::shared_ptr<const Acts::DigitizationModule> trapezoidalDigiModuleXZ(
    double minHalflengthX, double maxHalflengthX, double halflengthZ,
    double thickness, double gridSizeX, double gridSizeZ);

/// @brief Retrieves the node component from a mother by the string names
/// @param mother The Handle to the mother volume
/// @param nodeName The name of the note
/// @param attrName The name of the Atribute
/// @param attrValue The attribute value
dd4hep::xml::Component getNodeByStrAttr(const dd4hep::xml::Handle_t& mother,
                                        const std::string& nodeName,
                                        const std::string& attrName,
                                        const std::string& attrValue);

/// try to get attribute with double value, return defaultValue if attribute
/// not found
double getAttrValueWithFallback(const dd4hep::xml::Component& node,
                                const std::string& attrName,
                                const double& defaultValue);
}  // namespace utils
}  // namespace det
