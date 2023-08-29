(vertexing)=
# ACTS Vertexing Tutorial - Example: Adaptive Multi-Vertex Finder (AMVF) - Pythia8

This tutorial sets up and runs the ACTS Adaptive Multi-Vertex Finder on smeared ACTS Pythia8 truth tracks.
*Note*: You have to have Pythia8 available on your machine for following this tutorial.

## Prerequisites

For setting up and building ACTS, please refer to the general how-to ACTS guide. For this example you will need to enable Pythia8 in the examples by adding
```bash
-DACTS_BUILD_EXAMPLES_PYTHIA8=ON
```
to your cmake command.

Note: Additionally adding
```bash
-DCMAKE_BUILD_TYPE=Release
```
to your cmake command will significantly speed up the execution time of the vertexing example algorithm below, it will however slow down the compilation process.

## Setting up an ACTS Adaptive Multi-Vertex Finder Algorithm

A template algorithm file with an (almost) empty `execute()` method to be filled in the following is provided here:
```
../Examples/Algorithms/Vertexing/src/TutorialVertexFinderAlgorithm.cpp
```
Open the file in your editor and let's start setting up the AMVF. We will start setting up all necessary components in the `execute()` method.
*Note:* You would normally **not** want to do all the following setup steps in the `execute()` method (that is run on every single event), but rather in e.g. the constructor. For the sake of this tutorial, however, everything will be set up and run in the `execute()` method.

### Setting up required tools: Stepper and propagator

We need the `Acts::Propagator` with the `Acts::EigenStepper`:
```cpp
// Set up EigenStepper
Acts::EigenStepper<> stepper(m_cfg.bField);
// Set up the propagator
using Propagator = Acts::Propagator<Acts::EigenStepper<>>;
auto propagator = std::make_shared<Propagator>(stepper);
```

### Setting up required tools for the vertex fitter

Now, set up an impact point estimator...
```cpp
// Set up ImpactPointEstimator
using IPEstimator = Acts::ImpactPointEstimator<Acts::BoundTrackParameters, Propagator>;
IPEstimator::Config ipEstimatorCfg(m_cfg.bField, propagator);
IPEstimator ipEstimator(ipEstimatorCfg);
```
... and track linearizer for helical track parameters:
```cpp
// Set up the helical track linearizer
using Linearizer = Acts::HelicalTrackLinearizer<Propagator>;
Linearizer::Config ltConfig(m_cfg.bField, propagator);
Linearizer linearizer(ltConfig);
```
Now, for the sake of this example, let's specify a user-defined annealing scheme for the AVMF:
```cpp
// Set up deterministic annealing with user-defined temperatures
std::vector<double> temperatures{8.0, 4.0, 2.0, 1.4142136, 1.2247449, 1.0};
Acts::AnnealingUtility::Config annealingConfig;
annealingConfig.setOfTemperatures = temperatures;
Acts::AnnealingUtility annealingUtility(annealingConfig);
```
The AMVF strongly interplays with its dedicated vertex fitter, the *Adaptive Multi-Vertex Fitter*. Let's configure and set it up with the annealing utility defined above:
```cpp
// Set up the vertex fitter with user-defined annealing
using Fitter = Acts::AdaptiveMultiVertexFitter<Acts::BoundTrackParameters, Linearizer>;
Fitter::Config fitterCfg(ipEstimator);
fitterCfg.annealingTool = annealingUtility;
Fitter fitter(fitterCfg);
```

### Setting up required tools: Vertex seed finder

The last tool we need to set up (before finally setting up the AMVF) is a vertex seed finder:
```cpp
// Set up the vertex seed finder
using SeedFinder = Acts::TrackDensityVertexFinder<Fitter, Acts::GaussianTrackDensity<Acts::BoundTrackParameters>>;
SeedFinder seedFinder;
```
### Setting up the AMVF tool

Now we are ready to set up the Adaptive Multi-Vertex Finder. ACTS vertex finders are templated on the vertex fitter and vertex seed finder type:
```cpp
// The vertex finder type
using Finder = Acts::AdaptiveMultiVertexFinder<Fitter, SeedFinder>;
```
We configure the vertex finder in such a way that we do *not* use a beam spot constraint here:
```cpp
Finder::Config finderConfig(std::move(fitter), seedFinder, ipEstimator, linearizer);
```
Create the AMVF instance and a finder state to be passed to the `find()` method below:
```cpp
// Instantiate the finder
Finder finder(finderConfig);
// The vertex finder state
Finder::State state;
```
Lastly, we need to provide vertexing options. Here, we could e.g. set a beam spot constraint to the vertexing.
```cpp
// Default vertexing options, this is where e.g. a constraint could be set
using VertexingOptions = Acts::VertexingOptions<Acts::BoundTrackParameters>;
VertexingOptions finderOpts(ctx.geoContext, ctx.magFieldContext);
 ```
### Deploying the vertex finder on the track collection

Now we're ready to actually use the AMVF tool that we have set up above to find vertices on our input track collection. The `find()` methods on ACTS vertex finders return an `Acts::Result` object that we can use to check if any errors occurred and to retrieve the vertex collection:
```cpp
// Find vertices
auto res = finder.find(inputTrackPointers, finderOpts, state);

if (res.ok()) {
  // Retrieve vertices found by vertex finder
  auto vertexCollection = *res;
  ACTS_INFO("Found " << vertexCollection.size() << " vertices in event.");
  
  unsigned int count = 0;
  for (const auto& vtx : vertexCollection) {
    ACTS_INFO("\t" << ++count << ". vertex at "
                   << "(" << vtx.position().x() << "," << vtx.position().y()
                   << "," << vtx.position().z() << ") with "
                   << vtx.tracks().size() << " tracks.");
  }
} else {
  ACTS_ERROR("Error in vertex finder: " << res.error().message());
}
```
For reference, the full tutorial code can also be found in a file called `AdaptiveMultiVertexFinderAlgorithm.cpp` in the same directory as `TutorialVertexFinderAlgorithm.cpp`.

## Running the example algorithm

In your build directory, recompile and run the example on three pileup-50 pythia events to get your first ACTS vertices:
```
make -j4
./bin/ActsTutorialVertexFinder --evg-pileup 50 -n 3
```
