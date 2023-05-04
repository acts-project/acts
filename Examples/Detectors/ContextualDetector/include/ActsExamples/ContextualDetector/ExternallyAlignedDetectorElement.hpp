#pragma once

#include "Acts/Definitions/Algebra.hpp"
#include "Acts/Geometry/GeometryIdentifier.hpp"
#include "Acts/Geometry/GeometryContext.hpp"
#include "Acts/Utilities/Logger.hpp"
#include "ActsExamples/ContextualDetector/AlignmentDecorator.hpp"
#include "ActsExamples/ContextualDetector/ExternallyAlignedDetectorElement.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "ActsExamples/Framework/IContextDecorator.hpp"
#include "Acts/Plugins/Identification/IdentifiedDetectorElement.hpp"
#include "Acts/Plugins/Identification/Identifier.hpp"
#include "Acts/Surfaces/Surface.hpp"
#include "ActsExamples/GenericDetector/GenericDetectorElement.hpp"

#include <map>
#include <vector>

namespace Acts {
class TrackingGeometry;
}

namespace ActsExamples {

namespace Contextual {

/// @brief A mockup service that rotates the modules in a
/// simple tracking geometry
///
/// It acts on the PayloadDetectorElement, i.e. the
/// geometry context carries the full transform store (payload)
class ExternalAlignmentDecorator : public AlignmentDecorator {
 public:
  /// @brief nested configuration struct
  struct Config : public AlignmentDecorator::Config {
    /// The trackng geometry
    std::shared_ptr<const Acts::TrackingGeometry> trackingGeometry = nullptr;
  };

  /// Constructor
  ///
  /// @param cfg Configuration struct
  /// @param logger The logging framework
  ExternalAlignmentDecorator(
      const Config& cfg,
      std::unique_ptr<const Acts::Logger> logger = Acts::getDefaultLogger(
          "ExternalAlignmentDecorator", Acts::Logging::INFO));
  ExternalAlignmentDecorator::ExternalAlignmentDecorator(const Config& cfg,
                                                       const std::string& logFileName,
                                                       std::unique_ptr<const Acts::Logger> logger)
    : AlignmentDecorator(cfg, std::move(logger)), m_cfg(cfg) {
  auto fileLogger = std::make_unique<Acts::FileLogger>(logFileName);
  logger().addLogger(std::move(fileLogger));
  }

  /// Virtual destructor
  ~ExternalAlignmentDecorator() override = default;

  /// @brief decorates (adds, modifies) the AlgorithmContext
  /// with a geometric rotation per event
  ///
  /// @note If decorators depend on each other, they have to be
  /// added in order.
  ///
  /// @param context the bare (or at least non-const) Event context
  ProcessCode decorate(AlgorithmContext& context) override;

  /// @brief decorator name() for screen output
  const std::string& name() const override { return m_name; }

 private:
  Config m_cfg;                                  ///< the configuration class
  std::unique_ptr<const Acts::Logger> m_logger;  ///!< the logging instance
  std::string m_name = "ExternalAlignmentDecorator";

  /// Map of nominal transforms
  std::vector<Acts::Transform3> m_nominalStore;

  std::unordered_map<
      unsigned int,
      std::shared_ptr<ExternallyAlignedDetectorElement::AlignmentStore>>
      m_activeIovs;

  std::mutex m_iovMutex;

  size_t m_eventsSeen{0};

  /// Private access to the logging instance
  const Acts::Logger& logger() const { return *m_logger; }

  /// Populate the nominal transforms
  /// this parses the TrackingGeometry and fills the nominal store
  ///
  /// @param tGeometry the tracking geometry
  void parseGeometry(const Acts::TrackingGeometry& tGeometry);
};
}  // namespace Contextual


template <typename SurfaceT>
void fillTransforms(std::map<size_t, Acts::Transform3>& aStore, 
		     const Acts::TrackingGeometry& trackingGeometry, 
		     const SurfaceT* surface, 
		     const Acts::GeometryContext& nominalCtx, 
		     std::vector<ResultList>& resultList)
{
	int layer = surface -> layer(); 
	int volume = surface -> volume(); 
	
	Acts::Transform3 result=surface->transform(nominalCtx); 
	
	struct ResultList{
		int layer; 
		int volume; 
		Acts::Transform3 result; 
	};
	
	
	ResultList list; 
	list.layer = layer; 
	list.volume = volume; 
	list.result = result; 
	resultList.push_back(list); 
	
	
	const auto* detectorElement = surface -> associatedDetectorElement(); 
	auto alignableElement = dynamic_cast<const Acts::ExternallyAlignedDetectorElement*>(detectorElement); 
	if(alignableElement) { 
	     aStore[alignableElement->identifier()] = result;
	 }
}  

class ExternallyAlignedDetectorElement
    : public Generic::GenericDetectorElement {
 public:
  struct AlignmentStore {
    // GenericDetector identifiers are sequential
    std::vector<Acts::Transform3> transforms;
    size_t lastAccessed = 0;
  };

  /// @class ContextType
  /// convention: nested to the Detector element
  struct ContextType {
    // GenericDetector identifiers are an integer sequence, so vector is fine!
    std::shared_ptr<const AlignmentStore> alignmentStore{nullptr};
  };

  using Generic::GenericDetectorElement::GenericDetectorElement;

  /// Return local to global transform associated with this identifier
  ///
  /// @param gctx The current geometry context object, e.g. alignment
  ///
  /// @note this is called from the surface().transform(gctx)
  const Acts::Transform3& transform(
      const Acts::GeometryContext& gctx) const override;
};

inline const Acts::Transform3& ExternallyAlignedDetectorElement::transform(
    const Acts::GeometryContext& gctx) const {
  if (!gctx.hasValue()) {  // Treating empty context => nominal alignment
    return GenericDetectorElement::transform(gctx);
  }
  // cast into the right context object
  const auto& alignContext = gctx.get<ContextType>();
  identifier_type idValue = identifier_type(identifier());

  if (alignContext.alignmentStore == nullptr) {
    // geometry construction => nominal alignment
    return GenericDetectorElement::transform(gctx);
  }

  // At this point, the alignment store should be populated
  assert(idValue < alignContext.alignmentStore->transforms.size());
  return alignContext.alignmentStore->transforms[idValue];
}

}  // end of namespace Contextual
}  // end of namespace ActsExamples

