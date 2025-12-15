# Boost-Histogram Refactoring Strategy

## Executive Summary

This document outlines a strategy to replace ROOT histogram usage in the ACTS Examples framework with boost-histogram, while maintaining backward compatibility and enabling multiple output formats (ROOT, CSV, HDF5, etc.).

## Current Architecture Analysis

### Core Components

1. **PlotHelpers Abstraction** (`Examples/Framework/include/ActsExamples/Utilities/Helpers.hpp`)
   - `PlotHelpers::Binning`: Abstraction for histogram binning (uniform, variable, logarithmic)
   - Helper functions: `bookHisto()`, `fillHisto()`, `bookEff()`, `fillEff()`, `bookProf()`, `fillProf()`, `anaHisto()`
   - Currently returns raw ROOT histogram pointers (TH1F*, TH2F*, TEfficiency*, TProfile*)

2. **Plot Tools** (`Examples/Framework/include/ActsExamples/Validation/`)
   - ResPlotTool: Residual and pull distributions
   - EffPlotTool: Tracking efficiency plots
   - FakePlotTool: Fake ratio plots
   - DuplicationPlotTool: Duplication ratio plots
   - TrackSummaryPlotTool: Track summary plots
   - TrackQualityPlotTool: Completeness/purity plots

   **Common Pattern**:
   ```cpp
   class XxxPlotTool {
     struct Config { /* binning configurations */ };
     struct Cache { /* raw ROOT histogram pointers */ };

     void book(Cache& cache) const;        // Create histograms
     void fill(Cache& cache, ...) const;   // Fill histograms
     void refinement(Cache& cache) const;  // Post-processing (optional)
     void write(const Cache& cache) const; // Write to ROOT file
     void clear(Cache& cache) const;       // Delete histograms
   };
   ```

3. **Performance Writers** (`Examples/Io/Root/`)
   - 28 ROOT writer implementations
   - Inherit from `WriterT<T>`
   - Own plot tools and caches
   - Manage TFile lifecycle
   - Call plot tool methods in constructor/writeT()/finalize()/destructor

### Dependencies

- **Examples/Framework**: Links against `ROOT::Core` and `ROOT::Hist`
- **Python**: Already has `boost-histogram` via `hist` package (requirements.txt)
- **C++ boost-histogram**: NOT currently a dependency (would need to be added)

### Key Histogram Types

| ROOT Type | Usage | Special Features |
|-----------|-------|------------------|
| TH1F | 1D distributions | Standard histogram |
| TH2F | 2D scatter plots | Standard 2D histogram |
| TProfile | Profile plots | Computes mean/width per bin automatically |
| TEfficiency | Efficiency plots | Proper Clopper-Pearson confidence intervals |
| TH1D | 1D projections | Used for Gaussian fitting in refinement |

### Critical Features to Preserve

1. **TEfficiency**: Proper efficiency error calculation (not just sqrt(N) statistics)
2. **TProfile**: Automatic mean/width computation in bins
3. **Refinement Step**: Gaussian fitting on 2D histogram projections (ResPlotTool)
4. **Thread Safety**: Mutex-protected histogram filling
5. **Memory Management**: Clear ownership and lifecycle

## Revised Architecture (Hybrid Approach)

### Core Principle: Separation of Concerns

**Data Collection**: boost::histogram (fast, type-safe, modern C++)
**Analytics**: ROOT (proven statistical methods, Gaussian fits, efficiency calculations)
**Serialization**: ROOT TFile (backward compatible, can add CSV later)

This hybrid approach:
- Minimizes risk by keeping ROOT for complex statistics
- Improves code quality with modern C++ for data collection
- Maintains backward compatibility
- Allows gradual migration

### Phase 1: Histogram Abstraction Layer

Create a histogram wrapper that uses boost::histogram internally but converts to ROOT for analytics:

```
Examples/Framework/include/ActsExamples/Histogram/
├── Histogram.hpp              # Abstract histogram interface
├── Histogram1D.hpp            # 1D histogram wrapper
├── Histogram2D.hpp            # 2D histogram wrapper
├── ProfileHistogram.hpp       # Profile histogram wrapper
├── EfficiencyHistogram.hpp    # Efficiency histogram wrapper
└── HistogramImpl.hpp          # Implementation details
```

**Design Principles**:
- Use boost-histogram internally for data storage
- Provide a unified API that works for all histogram types
- Support move semantics and modern C++ ownership (no raw pointers)
- Type-safe and exception-safe
- Enable serialization to multiple formats

**Key Classes**:

```cpp
namespace ActsExamples::Histogram {

// Binning configuration (already exists as PlotHelpers::Binning)
using Binning = PlotHelpers::Binning;

// Abstract base for all histogram types
class HistogramBase {
 public:
  virtual ~HistogramBase() = default;
  virtual std::string name() const = 0;
  virtual std::string title() const = 0;
};

// 1D histogram
class Histogram1D : public HistogramBase {
 public:
  Histogram1D(std::string name, std::string title, const Binning& binning);

  void fill(double value, double weight = 1.0);

  // Conversion to ROOT (for backward compatibility)
  std::unique_ptr<TH1F> toROOT() const;

  // Direct access to boost-histogram
  const auto& histogram() const { return m_hist; }

 private:
  std::string m_name;
  std::string m_title;
  // boost-histogram storage
  boost::histogram::histogram<
    boost::histogram::axis::variable<>,
    boost::histogram::unlimited_storage<>
  > m_hist;
};

// 2D histogram
class Histogram2D : public HistogramBase {
 public:
  Histogram2D(std::string name, std::string title,
              const Binning& xBinning, const Binning& yBinning);

  void fill(double xValue, double yValue, double weight = 1.0);

  // Extract 1D projection for refinement
  Histogram1D projectX(int yBin) const;

  std::unique_ptr<TH2F> toROOT() const;

 private:
  std::string m_name;
  std::string m_title;
  boost::histogram::histogram<
    std::tuple<boost::histogram::axis::variable<>,
               boost::histogram::axis::variable<>>,
    boost::histogram::unlimited_storage<>
  > m_hist;
};

// Profile histogram (delegates to 2D histogram + ROOT conversion)
class ProfileHistogram : public HistogramBase {
 public:
  ProfileHistogram(std::string name, std::string title,
                   const Binning& xBinning, const Binning& yBinning);

  void fill(double xValue, double yValue, double weight = 1.0);

  // Convert to ROOT TProfile for mean/width calculation
  std::unique_ptr<TProfile> toROOT() const;

 private:
  std::string m_name;
  std::string m_title;
  Binning m_xBinning;
  Binning m_yBinning;

  // Store as 2D histogram in boost, convert to TProfile for analytics
  boost::histogram::histogram<
    std::tuple<boost::histogram::axis::variable<>,
               boost::histogram::axis::variable<>>,
    boost::histogram::unlimited_storage<>
  > m_hist;
};

// Efficiency histogram (collects passed/total, delegates ratio to ROOT)
class EfficiencyHistogram : public HistogramBase {
 public:
  EfficiencyHistogram(std::string name, std::string title,
                      const Binning& binning);

  EfficiencyHistogram(std::string name, std::string title,
                      const Binning& xBinning, const Binning& yBinning);

  void fill(double value, bool passed);
  void fill(double xValue, double yValue, bool passed);  // 2D overload

  // Convert to ROOT TEfficiency (computes ratio + Clopper-Pearson errors)
  std::unique_ptr<TEfficiency> toROOT() const;

 private:
  std::string m_name;
  std::string m_title;
  bool m_is2D;

  // Store passed and total separately, let ROOT compute efficiency
  Histogram1D m_passed;
  Histogram1D m_total;

  // 2D variants (if m_is2D == true)
  std::optional<Histogram2D> m_passed2D;
  std::optional<Histogram2D> m_total2D;
};

}  // namespace ActsExamples::Histogram
```

### Phase 2: Update PlotHelpers Interface

Modify `Examples/Framework/include/ActsExamples/Utilities/Helpers.hpp` to return histogram objects:

```cpp
namespace ActsExamples::PlotHelpers {

// Binning stays the same (already good abstraction)
class Binning { /* ... existing ... */ };

// Updated booking functions return histogram objects
Histogram::Histogram1D bookHisto(const std::string& name,
                                 const std::string& title,
                                 const Binning& binning);

Histogram::Histogram2D bookHisto(const std::string& name,
                                 const std::string& title,
                                 const Binning& xBinning,
                                 const Binning& yBinning);

Histogram::ProfileHistogram bookProf(const std::string& name,
                                     const std::string& title,
                                     const Binning& xBinning,
                                     const Binning& yBinning);

Histogram::EfficiencyHistogram bookEff(const std::string& name,
                                       const std::string& title,
                                       const Binning& binning);

// Fill functions work directly with histogram objects
void fillHisto(Histogram::Histogram1D& hist, float value, float weight = 1.0);
void fillHisto(Histogram::Histogram2D& hist, float xValue, float yValue,
               float weight = 1.0);
void fillProf(Histogram::ProfileHistogram& prof, float xValue, float yValue,
              float weight = 1.0);
void fillEff(Histogram::EfficiencyHistogram& eff, float value, bool status);

// Refinement helpers
void anaHisto(const Histogram::Histogram1D& inputHist, int j,
              Histogram::Histogram1D& meanHist,
              Histogram::Histogram1D& widthHist);

}  // namespace ActsExamples::PlotHelpers
```

### Phase 3: Update Plot Tool Caches

Modify plot tools to store histogram objects instead of raw pointers:

**Before**:
```cpp
struct Cache {
  std::map<std::string, TH1F*> res;           // raw pointers!
  std::map<std::string, TH2F*> res_vs_eta;
  TEfficiency* trackEff_vs_eta{nullptr};
};
```

**After**:
```cpp
struct Cache {
  std::map<std::string, Histogram::Histogram1D> res;        // value semantics!
  std::map<std::string, Histogram::Histogram2D> res_vs_eta;
  Histogram::EfficiencyHistogram trackEff_vs_eta;
};
```

**Benefits**:
- Automatic memory management (no manual delete)
- Move semantics (efficient)
- Type-safe
- No null pointer issues

### Phase 4: Output Serialization Layer

Create pluggable serializers for different output formats:

```
Examples/Framework/include/ActsExamples/Histogram/
└── Serializers/
    ├── HistogramSerializer.hpp      # Abstract interface
    ├── RootSerializer.hpp           # ROOT TFile output
    ├── CsvSerializer.hpp            # CSV output
    └── Hdf5Serializer.hpp           # HDF5 output (future)
```

**Abstract Interface**:

```cpp
namespace ActsExamples::Histogram {

class HistogramSerializer {
 public:
  virtual ~HistogramSerializer() = default;

  virtual void open(const std::string& path) = 0;
  virtual void write(const HistogramBase& hist) = 0;
  virtual void close() = 0;
};

// ROOT implementation
class RootSerializer : public HistogramSerializer {
 public:
  void open(const std::string& path) override;
  void write(const HistogramBase& hist) override;
  void close() override;

 private:
  std::unique_ptr<TFile> m_file;
};

// CSV implementation
class CsvSerializer : public HistogramSerializer {
 public:
  void open(const std::string& path) override;
  void write(const HistogramBase& hist) override;
  void close() override;

 private:
  std::filesystem::path m_outputDir;
};

}  // namespace ActsExamples::Histogram
```

### Phase 5: Update Plot Tools

Modify the `write()` method to use serializers:

**Before**:
```cpp
void ResPlotTool::write(const Cache& cache) const {
  for (auto& [name, hist] : cache.res) {
    hist->Write();  // Direct ROOT call
  }
}
```

**After**:
```cpp
void ResPlotTool::write(const Cache& cache,
                        HistogramSerializer& serializer) const {
  for (auto& [name, hist] : cache.res) {
    serializer.write(hist);
  }
}
```

### Phase 6: Update Writers

Writers choose the serializer based on configuration:

```cpp
class RootTrackFitterPerformanceWriter {
  struct Config {
    std::string inputTracks;
    std::string inputParticles;
    std::string inputTrackParticleMatching;
    std::string filePath = "performance_track_fitter.root";

    // New: output format selection
    enum class OutputFormat { ROOT, CSV, HDF5 };
    OutputFormat outputFormat = OutputFormat::ROOT;

    ResPlotTool::Config resPlotToolConfig;
    EffPlotTool::Config effPlotToolConfig;
  };

  ProcessCode finalize() override {
    m_resPlotTool.refinement(m_resPlotCache);

    // Create serializer based on config
    auto serializer = createSerializer(m_cfg.outputFormat, m_cfg.filePath);
    serializer->open(m_cfg.filePath);

    m_resPlotTool.write(m_resPlotCache, *serializer);
    m_effPlotTool.write(m_effPlotCache, *serializer);

    serializer->close();
    return ProcessCode::SUCCESS;
  }
};
```

## Testing Strategy

### Unit Tests

Create comprehensive unit tests for the histogram abstraction layer:

```
Examples/Framework/test/Histogram/
├── Histogram1DTests.cpp
├── Histogram2DTests.cpp
├── ProfileHistogramTests.cpp
├── EfficiencyHistogramTests.cpp
└── SerializerTests.cpp
```

**Test Coverage**:
1. **Histogram Creation**: Test binning (uniform, variable, logarithmic)
2. **Filling**: Test fill operations with various weights
3. **Statistics**: Verify mean, variance, sum, underflow, overflow
4. **Projections**: Test 2D to 1D projections
5. **Conversions**: Verify boost-histogram ↔ ROOT conversions
6. **Thread Safety**: Multi-threaded fill operations
7. **Memory**: No leaks, proper RAII

### Conversion Tests

Critical tests to ensure boost-histogram → ROOT conversion preserves data:

```cpp
TEST(HistogramConversion, Histogram1D_MatchesROOT) {
  // Create boost histogram
  auto boostHist = PlotHelpers::bookHisto("test", "Test",
      Binning::Uniform("x", 100, 0, 10));

  // Fill with known data
  std::mt19937 rng(42);
  std::normal_distribution<double> dist(5.0, 1.0);
  for (int i = 0; i < 10000; ++i) {
    boostHist.fill(dist(rng));
  }

  // Convert to ROOT
  auto rootHist = boostHist.toROOT();

  // Verify bin-by-bin equality
  BOOST_CHECK_EQUAL(boostHist.histogram().size(), rootHist->GetNbinsX());
  for (int i = 1; i <= rootHist->GetNbinsX(); ++i) {
    double boostContent = boostHist.histogram().at(i-1);
    double rootContent = rootHist->GetBinContent(i);
    BOOST_CHECK_CLOSE(boostContent, rootContent, 1e-10);
  }

  // Verify statistics
  BOOST_CHECK_CLOSE(boostHist.mean(), rootHist->GetMean(), 1e-6);
  BOOST_CHECK_CLOSE(boostHist.stddev(), rootHist->GetStdDev(), 1e-6);
}

TEST(HistogramConversion, Efficiency_MatchesROOT_ClopperPearson) {
  // Create boost efficiency histogram
  auto boostEff = PlotHelpers::bookEff("test_eff", "Test Efficiency",
      Binning::Uniform("x", 10, 0, 10));

  // Fill with known pass/fail pattern
  std::vector<std::pair<double, bool>> data = {
    {0.5, true}, {1.5, false}, {2.5, true}, {2.5, true},
    {3.5, false}, {3.5, false}, {3.5, true}, // ... etc
  };
  for (auto [x, passed] : data) {
    boostEff.fill(x, passed);
  }

  // Convert to ROOT
  auto rootEff = boostEff.toROOT();

  // Verify efficiency values and Clopper-Pearson confidence intervals
  for (int i = 1; i <= 10; ++i) {
    double boostEffVal = boostEff.efficiency(i);
    double rootEffVal = rootEff->GetEfficiency(i);
    BOOST_CHECK_CLOSE(boostEffVal, rootEffVal, 1e-6);

    double boostErrLow = boostEff.errorLow(i);
    double rootErrLow = rootEff->GetEfficiencyErrorLow(i);
    BOOST_CHECK_CLOSE(boostErrLow, rootErrLow, 1e-6);

    double boostErrHigh = boostEff.errorHigh(i);
    double rootErrHigh = rootEff->GetEfficiencyErrorHigh(i);
    BOOST_CHECK_CLOSE(boostErrHigh, rootErrHigh, 1e-6);
  }
}

TEST(HistogramConversion, Profile_MatchesROOT) {
  // Test that profile histograms compute same mean/width as ROOT
  auto boostProf = PlotHelpers::bookProf("test_prof", "Test Profile",
      Binning::Uniform("x", 10, 0, 10),
      Binning::Uniform("y", 100, -5, 5));

  // Fill with known x→y relationship
  std::mt19937 rng(42);
  for (double x = 0; x < 10; x += 0.1) {
    std::normal_distribution<double> dist(x*0.5, 1.0);  // y = 0.5*x + noise
    for (int i = 0; i < 100; ++i) {
      boostProf.fill(x, dist(rng));
    }
  }

  auto rootProf = boostProf.toROOT();

  // Verify mean and spread in each bin
  for (int i = 1; i <= 10; ++i) {
    BOOST_CHECK_CLOSE(boostProf.mean(i), rootProf->GetBinContent(i), 1e-6);
    BOOST_CHECK_CLOSE(boostProf.error(i), rootProf->GetBinError(i), 1e-6);
  }
}
```

### Integration Tests

Test the full pipeline with existing performance writers:

```cpp
TEST(Integration, TrackFitterPerformanceWriter_ProducesIdenticalOutput) {
  // Run the same test with:
  // 1. Original ROOT-based implementation
  // 2. New boost-histogram-based implementation

  // Compare output ROOT files bin-by-bin
  // Should be bitwise identical
}
```

### Regression Tests

Use existing physmon workflows as regression tests:

```python
# CI/physmon/workflows/physmon_trackfitting_gsf.py
# Add a check that compares output histograms to reference
def test_histogram_consistency(output_file, reference_file):
    """Ensure boost-histogram produces same results as ROOT version"""
    import uproot

    output = uproot.open(output_file)
    reference = uproot.open(reference_file)

    for key in reference.keys():
        hist_out = output[key]
        hist_ref = reference[key]

        # Check bin contents
        assert np.allclose(hist_out.values(), hist_ref.values(), rtol=1e-10)

        # Check bin errors
        assert np.allclose(hist_out.errors(), hist_ref.errors(), rtol=1e-10)
```

## Concrete Implementation Plan (Non-Breaking)

### Key Insight: Preserve External API

The refactoring can be **100% non-breaking** by keeping the same API at the writer level:
- PlotTools still have `book()`, `fill()`, `refinement()`, `write()`, `clear()` methods
- Writers still call these methods in the same way
- Only internal Cache storage changes from `TH1F*` to `Histogram1D`
- Only internal implementation of `write()` changes (convert boost→ROOT before writing)

### Milestone 1: Histogram Wrapper Foundation (1-2 weeks)
**Goal**: Create boost::histogram wrappers with ROOT conversion

**Tasks**:
- [ ] Create `Examples/Framework/include/ActsExamples/Histogram/` directory structure
- [ ] Implement `Histogram1D` wrapper around `boost::histogram::histogram<axis::variable<>>`
- [ ] Implement `Histogram2D` wrapper around `boost::histogram::histogram<tuple<axis::variable<>, axis::variable<>>>`
- [ ] Implement `toROOT()` conversion methods (boost → TH1F/TH2F)
- [ ] Write comprehensive unit tests for:
  - Histogram creation with uniform/variable/logarithmic binning
  - Fill operations with weights
  - Underflow/overflow handling
  - Bin-by-bin comparison with ROOT after conversion
  - Statistical properties (mean, stddev, integral)

**Deliverable**: `Histogram1D` and `Histogram2D` classes with proven ROOT conversion

**Risk**: None - no existing code changes

### Milestone 2: Special Histogram Types (1 week)
**Goal**: Add ProfileHistogram and EfficiencyHistogram

**Tasks**:
- [ ] Implement `ProfileHistogram` (stores 2D boost::histogram, converts to TProfile)
- [ ] Implement `EfficiencyHistogram` (stores passed/total separately, converts to TEfficiency)
- [ ] Write conversion tests verifying:
  - TProfile mean/width matches ROOT computation
  - TEfficiency ratios and Clopper-Pearson errors match ROOT
- [ ] Performance benchmark: boost fill vs ROOT fill

**Deliverable**: All 4 histogram types working with ROOT conversion

**Risk**: None - still no existing code changes

### Milestone 3: PlotHelpers Integration (1 week)
**Goal**: Make PlotHelpers use new histogram classes internally

**Tasks**:
- [ ] Update `PlotHelpers::bookHisto()` implementations to use `Histogram1D`/`Histogram2D`
- [ ] Update `PlotHelpers::bookProf()` to use `ProfileHistogram`
- [ ] Update `PlotHelpers::bookEff()` to use `EfficiencyHistogram`
- [ ] Update `PlotHelpers::fillHisto()` etc. to work with new types
- [ ] Update `PlotHelpers::anaHisto()` to:
  - Take boost histogram as input
  - Convert to ROOT TH1D for Gaussian fitting
  - Extract fit results to output histograms

**API Strategy**: Two approaches -

**Option A - Keep Raw Pointers (Zero Breaking Changes)**:
```cpp
// PlotHelpers returns pointers to wrapper objects
Histogram::Histogram1D* bookHisto(const std::string& name, ...);
void fillHisto(Histogram::Histogram1D* hist, float value, float weight);
```
Cache structs change from `TH1F*` to `Histogram::Histogram1D*` - still pointers!

**Option B - Use Value Semantics (Clean but requires Cache updates)**:
```cpp
// PlotHelpers returns objects by value
Histogram::Histogram1D bookHisto(const std::string& name, ...);
void fillHisto(Histogram::Histogram1D& hist, float value, float weight);
```
Cache structs change from `map<string, TH1F*>` to `map<string, Histogram::Histogram1D>`

**Recommendation**: Start with Option A (pointers) for minimal disruption, migrate to Option B later

**Deliverable**: PlotHelpers working with boost::histogram internally

**Risk**: Low - PlotHelpers is an internal implementation detail

### Milestone 4: Single Plot Tool Migration (1 week)
**Goal**: Migrate one plot tool as proof of concept

**Tasks**:
- [ ] Choose simplest tool: **EffPlotTool** (only uses TEfficiency)
- [ ] Update `EffPlotTool::Cache` to use `EfficiencyHistogram` instead of `TEfficiency*`
- [ ] Update `EffPlotTool::book()` to use `PlotHelpers::bookEff()` (returns new type)
- [ ] Update `EffPlotTool::fill()` - no changes needed if using pointers
- [ ] Update `EffPlotTool::write()` to:
  ```cpp
  void EffPlotTool::write(const Cache& cache) const {
    // Convert boost histograms to ROOT
    auto rootHist = cache.trackEff_vs_eta.toROOT();
    rootHist->Write();
    // Repeat for all histograms
  }
  ```
- [ ] Update `EffPlotTool::clear()` to delete wrapper objects (if using pointers)
- [ ] Run existing tests - should pass with no changes to test code

**Deliverable**: One plot tool fully migrated, tests passing

**Risk**: Medium - first real code change, but isolated to one tool

### Milestone 5: Remaining Plot Tools Migration (2 weeks)
**Goal**: Migrate all 6 plot tools

**Tasks**:
- [ ] Migrate **ResPlotTool** (most complex - has refinement step)
  - Update Cache to use `Histogram1D` and `Histogram2D`
  - Update `refinement()` method to convert boost→ROOT for Gaussian fitting
- [ ] Migrate **TrackSummaryPlotTool** (uses TProfile)
- [ ] Migrate **FakePlotTool**
- [ ] Migrate **DuplicationPlotTool**
- [ ] Migrate **TrackQualityPlotTool**
- [ ] Run full test suite after each migration
- [ ] Update unit tests for plot tools if needed

**Deliverable**: All plot tools using boost::histogram, tests passing

**Risk**: Medium-High - touches core validation code

### Milestone 6: Performance Writers (No Changes Needed!)
**Goal**: Verify writers work unchanged

**Tasks**:
- [ ] Run `RootTrackFitterPerformanceWriter` - should work unchanged
- [ ] Run `RootTrackFinderPerformanceWriter` - should work unchanged
- [ ] Spot-check several other writers

**Deliverable**: Writers confirmed working with no code changes

**Risk**: None - writers should be unaffected by internal plot tool changes

### Milestone 7: Integration Testing (1 week)
**Goal**: Validate full physics workflows

**Tasks**:
- [ ] Run all physmon workflows:
  - `physmon_trackfitting_gsf.py`
  - `physmon_trackrefitting_gsf.py`
  - All CKF, vertexing workflows
- [ ] Compare output ROOT files to reference (should be identical)
- [ ] Run full CI test suite
- [ ] Performance benchmarking (expect improvement in fill performance)

**Deliverable**: All workflows passing, performance characterized

**Risk**: Low - if plot tool tests pass, workflows should pass

### Milestone 8: Optional Enhancements (1-2 weeks)
**Goal**: Add CSV serialization (if desired)

**Tasks**:
- [ ] Implement `HistogramSerializer` abstract interface
- [ ] Implement `CsvSerializer` that writes histograms as CSV files
- [ ] Add config option to writers: `enum class OutputFormat { ROOT, CSV }`
- [ ] Update plot tool `write()` methods to take serializer parameter
- [ ] Add tests for CSV output

**Deliverable**: Optional CSV output capability

**Risk**: None - purely additive feature

### Milestone 9: API Cleanup (1 week, optional)
**Goal**: Migrate from pointers to value semantics

**Tasks**:
- [ ] Change PlotHelpers to return by value instead of pointer
- [ ] Update all Cache structs to use value semantics
- [ ] Remove manual `delete` calls in `clear()` methods
- [ ] Update tests

**Deliverable**: Modern C++ API with RAII

**Risk**: Medium - requires coordinated changes across all plot tools

### Total Timeline: 8-12 weeks

**Critical Path**:
- Weeks 1-2: Foundation
- Week 3: Special types
- Week 4: PlotHelpers
- Week 5: First plot tool
- Weeks 6-7: Remaining plot tools
- Week 8: Integration testing

**Optional**:
- Weeks 9-10: CSV serialization
- Week 11: API cleanup

## Migration Strategy

### Backward Compatibility

During migration, support both interfaces:

```cpp
// Old interface (deprecated)
TH1F* PlotHelpers::bookHisto_deprecated(const std::string& name, ...);

// New interface
Histogram::Histogram1D PlotHelpers::bookHisto(const std::string& name, ...);

// Conversion helper
TH1F* toROOT(const Histogram::Histogram1D& hist) {
  return hist.toROOT().release();
}
```

### Gradual Migration

1. Keep old implementation alongside new
2. Migrate one plot tool at a time
3. Run comparison tests to verify equivalence
4. Only remove old code after all migrations complete

## Benefits

### Immediate Benefits
1. **Decoupling**: Framework no longer depends on ROOT for histogram creation
2. **Type Safety**: Compile-time type checking, no raw pointers
3. **Memory Safety**: RAII, no manual memory management
4. **Performance**: boost-histogram is faster than ROOT for filling
5. **Flexibility**: Easy to add new output formats (CSV, HDF5, JSON)

### Long-term Benefits
1. **Maintainability**: Cleaner code, fewer dependencies
2. **Testability**: Easier to write unit tests without ROOT
3. **Portability**: Reduced ROOT dependency surface
4. **Interoperability**: Can process histograms in Python without ROOT

## Risks and Mitigations

### Risk 1: Incomplete Feature Parity
**Risk**: boost-histogram might not support all ROOT features
**Mitigation**: Identify all required features upfront, implement custom logic where needed (e.g., Clopper-Pearson for efficiency)

### Risk 2: Performance Regression
**Risk**: New implementation might be slower
**Mitigation**: Benchmark critical paths, optimize hot loops, boost-histogram is generally faster for filling

### Risk 3: Numerical Differences
**Risk**: Floating-point differences between implementations
**Mitigation**: Comprehensive comparison tests, define acceptable tolerance levels

### Risk 4: Breaking Changes
**Risk**: Refactoring breaks downstream users
**Mitigation**: Maintain backward compatibility during transition, provide deprecation warnings

## Technical Decisions (CONFIRMED)

1. **Histogram Library**: Use `boost::histogram` (header-only, part of Boost)
   - Already have Boost dependency (>= 1.77)
   - boost::histogram available since Boost 1.70
   - No additional dependencies needed

2. **Gaussian Fitting**: Keep ROOT's TFitResult for refinement step
   - Collect 2D histograms in boost::histogram
   - Convert to ROOT for Gaussian fitting
   - Minimizes risk, leverages proven ROOT implementation

3. **Efficiency Calculations**: Hybrid approach
   - Collect "passed" and "total" histograms separately in boost::histogram
   - Convert both to ROOT TH1
   - Use ROOT's TEfficiency for Clopper-Pearson confidence intervals
   - Avoids reimplementing complex statistical methods

4. **Profile Histograms**: Delegate to ROOT
   - Collect 2D histogram in boost::histogram
   - Convert to ROOT TProfile for mean/width computation
   - Alternative: Could use boost::histogram's weighted_mean accumulator (future optimization)

5. **CMake Dependencies**:
   - Require Boost >= 1.77 (current requirement)
   - boost::histogram available automatically
   - Keep ROOT dependency for analytical functions

6. **Output Format Configuration**:
   - Phase 1: Always output ROOT (no breaking changes)
   - Phase 2: Add optional CSV serializer
   - Phase 3: Make serialization pluggable

## References

- boost::histogram documentation: https://www.boost.org/doc/libs/release/libs/histogram/
- ROOT TEfficiency: https://root.cern.ch/doc/master/classTEfficiency.html
- Clopper-Pearson intervals: https://en.wikipedia.org/wiki/Binomial_proportion_confidence_interval
