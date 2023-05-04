#include "ActsExamples/ContextualDetector/ExternalAlignmentDecorator.hpp"
<<<<<<< HEAD

#include "Acts/Geometry/GeometryContext.hpp"
=======
>>>>>>> 2c318123d (adding hpp and cpp files-xternal alignment)
#include "Acts/Geometry/TrackingGeometry.hpp"
#include "Acts/Surfaces/SurfaceArray.hpp"
#include "ActsExamples/ContextualDetector/ExternallyAlignedDetectorElement.hpp"
#include "ActsExamples/Framework/AlgorithmContext.hpp"
#include "ActsExamples/Framework/RandomNumbers.hpp"

#include <cassert>
#include <ostream>
#include <thread>
#include <utility>


// ExternalAlignmentDecorator -- class that is responsible for applying the external alignment output to trackingGeometry
//                            -- basically based on the external alignment measurements this class can modify on the geometry of detector 
// TrackingGeometry  - geometry of the detector 

// this constructor ExternalAlignmentDecorator takes 2 parameters: Config and pointer to Logger object
// initalizing the variables: m_cfg  && m_logger (respect to some arguments)
// inside the if statement we're checking if m_cfg is not null, then we will parseGeometry function that will populates the member variable m_nominalStore  with transforms for the surfaces in the tracking geometry
ActsExamples::Contextual::ExternalAlignmentDecorator::
    ExternalAlignmentDecorator(const Config& cfg,
                               std::unique_ptr<const Acts::Logger> logger)
    : m_cfg(cfg), m_logger(std::move(logger)) {
  if (m_cfg.trackingGeometry != nullptr) {
    // parse and populate
    parseGeometry(*m_cfg.trackingGeometry.get());
  }
}


// decorate() - takes object AlgorithmContext as argument and returns an ActsExamples::ProcessCode enum
// m_iovMutex() - 
ActsExamples::ProcessCode
ActsExamples::Contextual::ExternalAlignmentDecorator::decorate(
    AlgorithmContext& context) {
  // Iov map access needs to be synchronized
  std::lock_guard lock{m_iovMutex};

  // In which iov batch are we?
  // IOV - Interval of Validity

  unsigned int iov = context.eventNumber / m_cfg.iovSize;
  ACTS_VERBOSE("IOV handling in thread " << std::this_thread::get_id() << ".");
  ACTS_VERBOSE("IOV resolved to " << iov << " - from event "
                                  << context.eventNumber << ".");

  m_eventsSeen++;

//  checking if m_cfg.randomNumberSvc is not null, we are doing check in if the current batch iov can be find in m_activeIovs map
// basically, if the IOV can be find 
  if (m_cfg.randomNumberSvc != nullptr) {
    if (auto it = m_activeIovs.find(iov); it != m_activeIovs.end()) {
      // Iov is already present, update last accessed
      it->second->lastAccessed = m_eventsSeen;
      context.geoContext =
          ExternallyAlignedDetectorElement::ContextType{it->second};
    } else {
      // Iov is not present yet, create it
      auto alignmentStore =
          std::make_unique<ExternallyAlignedDetectorElement::AlignmentStore>();
      alignmentStore->lastAccessed = m_eventsSeen;

      ACTS_VERBOSE("New IOV " << iov << " detected at event "
                              << context.eventNumber
                              << ", emulate new alignment.");

     

      // Create an algorithm local random number generator
      RandomEngine rng = m_cfg.randomNumberSvc->spawnGenerator(context);


      alignmentStore->transforms = m_nominalStore;  // copy nominal alignment
      for (auto& tForm : alignmentStore->transforms) {
        // Multiply alignment in place
        applyTransform(tForm, m_cfg, rng, iov);
      }

      auto [insertIterator, inserted] =
          m_activeIovs.emplace(iov, std::move(alignmentStore));
      assert(inserted && "Expected IOV to be created in map, but wasn't");

      // make context from iov pointer, address should be stable
      context.geoContext =
          ExternallyAlignedDetectorElement::ContextType{insertIterator->second};
    }
  }

  // Garbage collection
  if (m_cfg.doGarbageCollection) {
    for (auto it = m_activeIovs.begin(); it != m_activeIovs.end();) {
      auto& status = it->second;
      if (m_eventsSeen - status->lastAccessed > m_cfg.flushSize) {
        ACTS_DEBUG("IOV " << iov << " has not been accessed in the last "
                          << m_cfg.flushSize << " events, clearing");
        it = m_activeIovs.erase(it);
      } else {
        it++;
      }
    }
  }

  return ProcessCode::SUCCESS;
}

void ActsExamples::Contextual::ExternalAlignmentDecorator::parseGeometry(
    const Acts::TrackingGeometry& tGeometry) {
  // Double-visit - first count
  size_t nTransforms = 0;
  tGeometry.visitSurfaces([&nTransforms](const auto*) { ++nTransforms; });


// this is where initialization is done - object nominalCtx (type Acts::GeometryContext) with 
// with an empty context type ExternallyAlignedDetectorElement::ContextType

// GeometryContext - object where we can find an info about geometry -- basically the nominalCtx holds the context geomtry without applying external alignment 
  Acts::GeometryContext nominalCtx{
      ExternallyAlignedDetectorElement::ContextType{}};
//L111: creating a 'empty context type' based on the geometry that has been used 
  
  // Collect the surfacas into the nominal store -- > or in other words, vector that stores info about nominal alignment for detectors elements
  // typical declaration of a  std::vector object aStore with paramateres (nTransforms - number of elements, type Acts::Transform3)
  // for each element, we have init identity transformation
  std::vector<Acts::Transform3> aStore(nTransforms,
                                       Acts::Transform3::Identity());

// function that takes pointer to a surface object, and return void 
// aStore - vector of Acts::Transform3 objects and nominalCtx passed by reference because it will change 
  auto fillTransforms = [&aStore, &nominalCtx](const auto* surface) -> void {
    // it first casts the associatedDetectorElement() of the surface pointer to a const ExternallyAlignedDetectorElement*, and then uses the identifier of that element to index into the aStore vector and set the corresponding Transform3 object to the result of calling transform() on the surface object with the nominalCtx context. Essentially, this lambda function is used to populate the m_nominalStore vector with the nominal alignment of the surfaces in the tracking geometry.
        dynamic_cast<const ExternallyAlignedDetectorElement*>(
            surface->associatedDetectorElement());
    aStore[alignableElement->identifier()] = surface->transform(nominalCtx);
  };


// TOdo task: figure out how identifier() can access to specific layer and volume, an safe that information to new list where can access it whenever we want 
// first attemp: 
 auto fillTransforms = [&aStore, &nominalCtx, &resultList](const auto* surface) -> void {
    
    // we have to somehow retrieve the vol + sur info from surface/volume object
    int layer = surface->layer();
    int volume = surface->volume();

    // setting up the result of calling a method transform() on the surface object 
    Transform3 result = surface->transform(nominalCtx);

    // defining the structure in which we could store the all the relevant informations (surface, volume, etc)
    struct ResultList {
        int layer;
        int volume;
        Transform3 result;
    };

    // all the important informations regarding vol, layers, transformation resutls can be save into the list 
   ResultList list;
    list.layer = layer;
    list.volume = volume;
    list.result = result;

    // Adding new results to the list
    resultList.push_back(list);

    // Index into the aStore vector and set the corresponding Transform3 object.
    dynamic_cast<const ExternallyAlignedDetectorElement*>(
            surface->associatedDetectorElement());
    aStore[alignableElement->identifier()] = result;
};


  tGeometry.visitSurfaces(fillTransforms);
  m_nominalStore = std::move(aStore);
}
