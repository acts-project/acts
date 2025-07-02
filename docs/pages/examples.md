# Code Examples {#examples_page}

@tableofcontents

This page provides practical examples of using the Acts tracking software for common tasks.

## Basic Usage Examples

### 1. Setting Up a Propagator

A propagator is essential for track extrapolation:

@code{.cpp}
#include "Acts/Propagator/EigenStepper.hpp"
#include "Acts/Propagator/Propagator.hpp"
#include "Acts/MagneticField/ConstantBField.hpp"

// Setup magnetic field
Acts::Vector3 bField(0.0, 0.0, 2.0);  // 2T in z-direction
auto magneticField = std::make_shared<Acts::ConstantBField>(bField);

// Create stepper
Acts::EigenStepper<> stepper(magneticField);

// Create propagator
Acts::Propagator<Acts::EigenStepper<>> propagator(std::move(stepper));

// Usage example
Acts::PropagatorOptions<> options;
auto result = propagator.propagate(startParameters, targetSurface, options);
@endcode

### 2. Track Fitting

Fitting tracks to measurements:

```cpp
#include "Acts/TrackFitting/KalmanFitter.hpp"

// Setup Kalman fitter
Acts::KalmanFitterOptions<Acts::VoidOutlierFinder> kfOptions;
kfOptions.geoContext = geometryContext;
kfOptions.magFieldContext = magneticFieldContext;
kfOptions.calibrationContext = calibrationContext;

// Create fitter
Acts::KalmanFitter<Acts::Propagator<Acts::EigenStepper<>>, Acts::VoidOutlierFinder> 
    fitter(propagator, logger);

// Fit track
auto result = fitter.fit(measurements, initialParameters, kfOptions);

if (result.ok()) {
    const auto& track = result.value();
    std::cout << "Chi2: " << track.chi2() << std::endl;
    std::cout << "NDF: " << track.nDoF() << std::endl;
}
```

### 3. Seeding

Finding track seeds from space points:

@code{.cpp}
#include "Acts/Seeding/SeedFinder.hpp"
#include "Acts/Seeding/SpacePointGrid.hpp"

// Configure seed finder
Acts::SeedFinderConfig<SpacePoint> config;
config.rMax = 200.0;           // mm
config.rMin = 33.0;            // mm  
config.deltaRMin = 5.0;        // mm
config.deltaRMax = 60.0;       // mm
config.collisionRegionMin = -250.0;  // mm
config.collisionRegionMax = 250.0;   // mm
config.zMin = -2800.0;         // mm
config.zMax = 2800.0;          // mm
config.maxSeedsPerSpM = 5;
config.cotThetaMax = 7.40627;  // 2.7 eta

// Create seed finder
Acts::SeedFinder<SpacePoint> seedFinder(config);

// Find seeds
std::vector<Acts::Seed<SpacePoint>> seeds;
auto groupIt = spacePointsGrouped.begin();
auto endIt = spacePointsGrouped.end();

for (; groupIt != endIt; ++groupIt) {
    seedFinder.createSeedsForGroup(groupIt->second.begin(), 
                                   groupIt->second.end(),
                                   std::back_inserter(seeds));
}
@endcode

## Advanced Examples

### 1. Custom Material Service

Implementing a material service:

```cpp
class CustomMaterialService : public Acts::IMaterialDecorator {
public:
    void decorateSurface(Acts::Surface& surface) const override {
        // Get surface position and type
        auto center = surface.center();
        auto type = surface.type();
        
        if (type == Acts::Surface::Cylinder) {
            // Silicon for tracker
            double radius = center.norm();
            if (radius < 100.0) {  // Inner tracker
                auto material = Acts::Material::fromMolarDensity(
                    28.0855, 14, 2.329, 21.82  // Silicon properties
                );
                surface.assignMaterial(
                    std::make_shared<Acts::MaterialSlab>(material, 0.32)
                );
            }
        } else if (type == Acts::Surface::Plane) {
            // Lead for calorimeter
            auto material = Acts::Material::fromMolarDensity(
                207.2, 82, 11.34, 18.26  // Lead properties  
            );
            surface.assignMaterial(
                std::make_shared<Acts::MaterialSlab>(material, 2.0)
            );
        }
    }
};
```

### 2. Event Data Model Integration

Integrating with an event data model:

@code{.cpp}
#include "Acts/EventData/TrackParameters.hpp"
#include "Acts/EventData/Measurement.hpp"

// Convert from EDM to Acts
class EDMConverter {
public:
    static Acts::BoundTrackParameters convertTrack(const EDM::Track& edmTrack) {
        // Extract parameters
        auto pos = edmTrack.referencePoint();
        auto mom = edmTrack.momentum();
        auto charge = edmTrack.charge();
        
        // Create surface (example: perigee)
        auto perigeeSurface = Acts::Surface::makeShared<Acts::PerigeeSurface>(pos);
        
        // Build parameter vector
        Acts::BoundVector params;
        params[Acts::eBoundLoc0] = 0.0;     // d0
        params[Acts::eBoundLoc1] = 0.0;     // z0  
        params[Acts::eBoundTime] = edmTrack.time();
        params[Acts::eBoundPhi] = mom.phi();
        params[Acts::eBoundTheta] = mom.theta();
        params[Acts::eBoundQOverP] = charge / mom.norm();
        
        return Acts::BoundTrackParameters(
            perigeeSurface, params, edmTrack.covariance()
        );
    }
    
    static std::vector<Acts::SourceLink> convertHits(
        const std::vector<EDM::Hit>& edmHits) {
        
        std::vector<Acts::SourceLink> sourceLinks;
        for (const auto& hit : edmHits) {
            sourceLinks.emplace_back(EDMSourceLink{hit});
        }
        return sourceLinks;
    }
};
@endcode

### 3. Multi-Threading Example

Thread-safe track processing:

```cpp
#include <thread>
#include <future>

class ParallelTrackProcessor {
private:
    Acts::Propagator<Acts::EigenStepper<>> m_propagator;
    Acts::KalmanFitter<Acts::Propagator<Acts::EigenStepper<>>, Acts::VoidOutlierFinder> m_fitter;
    
public:
    std::vector<Acts::Result<FittedTrack>> processEvents(
        const std::vector<Event>& events) {
        
        std::vector<std::future<Acts::Result<FittedTrack>>> futures;
        
        // Launch async tasks
        for (const auto& event : events) {
            futures.push_back(
                std::async(std::launch::async, [this, &event]() {
                    return this->processEvent(event);
                })
            );
        }
        
        // Collect results
        std::vector<Acts::Result<FittedTrack>> results;
        for (auto& future : futures) {
            results.push_back(future.get());
        }
        
        return results;
    }
    
private:
    Acts::Result<FittedTrack> processEvent(const Event& event) {
        // Create thread-local contexts
        Acts::GeometryContext geoCtx;
        Acts::MagneticFieldContext magCtx;
        Acts::CalibrationContext calCtx;
        
        // Process tracks in this event
        Acts::KalmanFitterOptions<Acts::VoidOutlierFinder> options;
        options.geoContext = geoCtx;
        options.magFieldContext = magCtx;
        options.calibrationContext = calCtx;
        
        return m_fitter.fit(event.measurements, event.initialParams, options);
    }
};
```

## Utility Examples

### 1. Debugging and Visualization

Helper functions for debugging:

@code{.cpp}
#include "Acts/Utilities/Logger.hpp"

// Custom logger setup
auto logger = Acts::getDefaultLogger("TrackFitter", Acts::Logging::DEBUG);

// Debug track parameters
void printTrackInfo(const Acts::BoundTrackParameters& params) {
    const auto& p = params.parameters();
    
    ACTS_DEBUG("Track parameters:");
    ACTS_DEBUG("  d0     = " << p[Acts::eBoundLoc0] << " mm");
    ACTS_DEBUG("  z0     = " << p[Acts::eBoundLoc1] << " mm");
    ACTS_DEBUG("  phi    = " << p[Acts::eBoundPhi] << " rad");
    ACTS_DEBUG("  theta  = " << p[Acts::eBoundTheta] << " rad");
    ACTS_DEBUG("  q/p    = " << p[Acts::eBoundQOverP] << " 1/MeV");
    ACTS_DEBUG("  time   = " << p[Acts::eBoundTime] << " ns");
}

// Geometry validation
void validateGeometry(const Acts::TrackingGeometry& geometry) {
    Acts::GeometryValidator validator;
    validator.validate(geometry);
}
@endcode

### 2. Performance Monitoring

Track processing performance:

```cpp
#include <chrono>

class PerformanceMonitor {
private:
    std::chrono::high_resolution_clock::time_point m_start;
    std::string m_operation;
    
public:
    PerformanceMonitor(const std::string& operation) 
        : m_operation(operation)
        , m_start(std::chrono::high_resolution_clock::now()) {}
    
    ~PerformanceMonitor() {
        auto end = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::microseconds>(
            end - m_start).count();
        
        std::cout << m_operation << " took " << duration << " μs" << std::endl;
    }
};

// Usage
void processTrack(const TrackData& data) {
    PerformanceMonitor monitor("Track fitting");
    
    // ... track fitting code ...
}
```

## Best Practices

### Error Handling

Always check results from Acts operations:

```cpp
auto result = fitter.fit(measurements, initialParams, options);

if (!result.ok()) {
    ACTS_ERROR("Track fitting failed: " << result.error());
    return;
}

const auto& track = result.value();
// Use fitted track...
```

### Memory Management

Use smart pointers for geometry objects:

@code{.cpp}
// Good: automatic cleanup
auto surface = Acts::Surface::makeShared<Acts::PlaneSurface>(transform, bounds);

// Avoid: manual memory management  
// Acts::PlaneSurface* surface = new Acts::PlaneSurface(transform, bounds);
@endcode

### Configuration

Use configuration objects for complex setups:

```cpp
// Centralized configuration
struct TrackerConfig {
    double magneticField = 2.0;    // Tesla
    double materialBudget = 0.1;   // X0
    size_t maxIterations = 5;
    double convergenceCut = 1e-6;
};

TrackerConfig config;
// ... setup based on config
```

---

@par Related Resources
- @ref core_concepts "Core Concepts"
- @ref geometry_guide "Geometry System"
- @ref getting_started "Getting Started"