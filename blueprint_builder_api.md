# Blueprint Builder API — Design Notes

## Status
- Renaming executed (LayerHelper→LayerAssembler etc.) ✓
- Declarative option: planned, not yet implemented

---

## Declarative Assembly Description — Plan

### Motivation

The current `onLayer` customizer in experiment code (`makeCustomizer` in ODD) is
mostly disguised configuration data: which binning to use, which envelope, which
navigation policy. A declarative `BarrelEndcapAssemblyDesc` struct exposes that
structure directly, eliminates boilerplate, and opens a path to serialization
(JSON/YAML detector descriptions).

### The struct

```cpp
// Pure data — no builder dependency, no construction logic
template<GeometryBackend B>
struct BarrelEndcapAssemblyDesc {
  std::string          assembly;       // name looked up via backend.findByName
  GeometryAxes         barrelAxes;
  GeometryAxes         endcapAxes;
  std::regex           layerFilter;

  Acts::ExtentEnvelope layerEnvelope = Acts::ExtentEnvelope::Zero();

  // BinCount: literal int or function of layer index (extracted from layerFilter
  // capture group 1 by default; see layerIndexExtractor below)
  using BinCount = std::variant<int, std::function<int(int layerIdx)>>;
  struct BarrelBins { BinCount phi; BinCount z; };
  struct EndcapBins { BinCount r;   BinCount phi; };
  BarrelBins barrelBins;
  EndcapBins endcapBins;

  Acts::VolumeAttachmentStrategy attachmentStrategy = VolumeAttachmentStrategy::Gap;
  Acts::VolumeResizeStrategy     resizeStrategy     = VolumeResizeStrategy::Gap;

  // Optional: override layer index extraction (default: regex capture group 1)
  std::optional<std::function<int(const typename B::Element& element)>>
      layerIndexExtractor;

  // Escape hatch: if set, replaces the generated onLayer customizer entirely
  std::optional<typename BarrelEndcapAssembler<B>::LayerCustomizer> onLayer;
};
```

### Call site (ODD example)

```cpp
builder.addBarrelEndcapAssembly(outer, {
    .assembly      = "Pixels",
    .barrelAxes    = "XYZ",
    .endcapAxes    = "XZY",
    .layerFilter   = std::regex{"(?:PixelLayer|PixelEndcap[NP])(\\d)"},
    .layerEnvelope = ExtentEnvelope::Zero().set(AxisZ, {2,2}).set(AxisR, {2,2}),
    .barrelBins    = { .phi = [&](int n){ return constant("pix_b{}_sf_b_phi", n); },
                       .z   = constant("pix_b_sf_b_z") },
    .endcapBins    = { .r   = constant("pix_e_sf_b_r"),
                       .phi = constant("pix_e_sf_b_phi") },
});
```

### Implementation approach — no duplication

The desc is a thin adapter over the fluent API. Construction logic lives once,
in `BarrelEndcapAssembler::build()`.

```
BarrelEndcapAssemblyDesc  ──wraps──>  BarrelEndcapAssembler  ──builds──>  node
      (pure data)                        (fluent API)              (construction)
```

`BlueprintBuilder` gets one new method:

```cpp
template<GeometryBackend B>
void BlueprintBuilder<B>::addBarrelEndcapAssembly(
    Acts::Experimental::BlueprintNode& parent,
    const BarrelEndcapAssemblyDesc<B>& desc)
{
  auto resolveBin = [](const typename BarrelEndcapAssemblyDesc<B>::BinCount& bc,
                       int idx) {
    return std::visit(Acts::overloaded{
        [](int v)           { return v; },
        [idx](const auto& f){ return f(idx); }
    }, bc);
  };

  auto extractIdx = [&desc](const auto& elem) -> int {
    if (desc.layerIndexExtractor)
      return (*desc.layerIndexExtractor)(elem);
    std::cmatch m;
    std::regex_search(elem.name(), m, desc.layerFilter);
    return m.size() > 1 ? std::stoi(m[1]) : 0;
  };

  LayerCustomizer layerCustomizer =
      desc.onLayer.value_or(
          [&](const auto& elem, auto& node) {   // void/mutable-ref form
            node->setEnvelope(desc.layerEnvelope);
            int idx = extractIdx(elem);
            SrfArrayNavPol::Config navCfg;
            if (m_backend.isBarrel(elem)) {
              navCfg.layerType = SrfArrayNavPol::LayerType::Cylinder;
              navCfg.bins = { resolveBin(desc.barrelBins.phi, idx),
                              resolveBin(desc.barrelBins.z,   idx) };
            } else {
              navCfg.layerType = SrfArrayNavPol::LayerType::Disc;
              navCfg.bins = { resolveBin(desc.endcapBins.r,   idx),
                              resolveBin(desc.endcapBins.phi, idx) };
            }
            node->setNavigationPolicyFactory(
                NavigationPolicyFactory{}
                    .add<CylinderNavigationPolicy>()
                    .add<SrfArrayNavPol>(navCfg)
                    .asUniquePtr());
          });

  barrelEndcap()
      .setAssembly(desc.assembly)
      .setAxes(desc.barrelAxes, desc.endcapAxes)
      .setLayerFilter(desc.layerFilter)
      .onLayer(std::move(layerCustomizer))
      .onContainer([&desc](const auto&, auto& node) {
          node->setAttachmentStrategy(desc.attachmentStrategy);
          node->setResizeStrategies(desc.resizeStrategy, desc.resizeStrategy);
      })
      .addTo(parent);
}
```

### Further step: backend constant specs (optional, Level 2)

Eliminates the remaining bin lambdas for the common DD4hep case:

```cpp
struct DD4hepConstant { std::string key; };  // {} replaced by layer index
using BinCount = std::variant<int, std::function<int(int)>, DD4hepConstant>;
```

Call site becomes fully literal (no lambdas):
```cpp
.barrelBins = { .phi = DD4hepConstant{"pix_b{}_sf_b_phi"},
                .z   = DD4hepConstant{"pix_b_sf_b_z"} },
```

Backend concept would need: `b.resolveConstant(key) -> int`.
Serialization (JSON/YAML) becomes straightforward once lambdas are gone.

---

## Naming decisions (applied)

| Old | New |
|-----|-----|
| `LayerHelper` | `LayerAssembler` |
| `BarrelEndcapAssemblyHelper` | `BarrelEndcapAssembler` |
| `layerHelper()` | `layers()` |
| `barrelEndcapAssemblyHelper()` | `barrelEndcap()` |
| `LayerAssembler::customize` | `onLayer` |
| `BarrelEndcapAssembler::customizeLayer` | `onLayer` |
| `BarrelEndcapAssembler::customize` | `onContainer` |
| `setLayerPattern` | `setLayerFilter` |

## Backend abstraction plan

### Status
- Not yet implemented; DD4hep-only first step agreed

### Motivation

`LayerAssembler` and `BarrelEndcapAssembler` currently hard-code `dd4hep::DetElement`
as the element type and `TGeoAxes` as the axis spec. Templating over a
`GeometryBackend` concept lets TGeo and GeoModel reuse the same assembler logic
with zero virtual dispatch and zero loss of native element types in customizers.

### Why explicit template instantiation, not type erasure

Type erasure (virtual interface + `std::any`) would require customizer
signatures to use `std::any` for the element parameter, destroying the type
safety that makes the API useful. Explicit template instantiation keeps
`Customizer` strongly typed (`const dd4hep::DetElement&` etc.) while still
hiding the implementation in `.cpp` files.

### File layout

```
Core/include/Acts/Geometry/
  GeometryBackend.hpp          // concept + assembler template DECLARATIONS only

Plugins/DD4hep/include/ActsPlugins/DD4hep/
  DD4hepBackend.hpp            // DD4hepBackend struct + extern template + using aliases
  BlueprintBuilder.hpp         // thin logger-aware wrapper (unchanged public API)

Plugins/DD4hep/src/
  BlueprintBuilder.cpp         // #include "GeometryBackend_impl.hpp"
                               // + explicit instantiation: template class LayerAssembler<DD4hepBackend>;

Plugins/Root/include/ActsPlugins/Root/   (or a shared header)
  GeometryAxes.hpp             // TGeoAxes lifted here, renamed GeometryAxes
```

`GeometryBackend_impl.hpp` (header included **only** from plugin `.cpp` files):
contains the `build()` method bodies for `LayerAssembler<B>` and
`BarrelEndcapAssembler<B>`. Never installed as a public include.

### GeometryBackend concept

```cpp
template<typename B>
concept GeometryBackend = requires(const B& b,
                                   const B::Element& elem,
                                   const std::regex& pat) {
  // Element type — native handle (dd4hep::DetElement, TGeoNode*, PVConstLink…)
  typename B::Element;

  // Axis spec — empty struct for GeoModel, {axes, layerAxes} for DD4hep/TGeo
  typename B::LayerSpec;

  // Search
  { b.world() }                              -> std::same_as<B::Element>;
  { b.nameOf(elem) }                         -> std::convertible_to<std::string_view>;
  { b.findByName(std::string_view{}) }       -> std::same_as<std::optional<B::Element>>;
  { b.findByPattern(elem, pat) }             -> std::same_as<std::vector<B::Element>>;
  { b.children(elem) }                       -> std::same_as<std::vector<B::Element>>;
  { b.isValid(elem) }                        -> std::same_as<bool>;

  // Classification
  { b.isBarrel(elem) }                       -> std::same_as<bool>;
  { b.isEndcap(elem) }                       -> std::same_as<bool>;
  { b.isSensitive(elem) }                    -> std::same_as<bool>;
};
```

### Handle semantics and lifetime model

This is the critical cross-backend constraint:

- DD4hep uses **value handles** (`dd4hep::DetElement`).
- TGeo will likely use **raw/non-owning pointers** (`TGeoNode*`-like handles).
- GeoModel uses **intrusive smart pointers** (`*ConstLink` handles).

To keep one assembler implementation working for all three, `B::Element` must
be treated as an opaque handle, never as an owning value object.

Rules for assembler templates:

- Store and pass `Element` by value/`const&` only. Do not assume ownership model.
- Do not dereference in generic code; all access goes through backend methods.
- Never compare handles directly in generic code (`==`, pointer identity, etc.).
  If identity is needed later, backend must provide it explicitly.
- `std::optional<Element>` is the only "not found" channel. For pointer
  backends, `findByName` must return `std::nullopt` instead of `optional{nullptr}`.
- Backend methods taking `const Element&` assume a valid handle; invalid/null
  checks happen at backend boundaries.

Practical mapping:

- `DD4hepBackend::Element = dd4hep::DetElement` (already value-handle style)
- `TGeoBackend::Element = const TGeoNode*` (or wrapper around it)
- `GeoModelBackend::Element = GeoModel::*ConstLink` (intrusive handle)

This keeps customizers strongly typed while avoiding a DD4hep-specific handle
assumption in the templated assembler logic.

`makeLayer` is **excluded** from the concept: DD4hep/TGeo take `(elem, LayerSpec)`
while GeoModel derives orientation from shape and takes no axis spec. Each
backend provides its own `makeLayer` overload outside the concept.

### DD4hepBackend struct

```cpp
// DD4hepBackend.hpp
struct DD4hepBackend {
  struct LayerSpec {
    TGeoAxes axes;
    std::optional<TGeoAxes> layerAxes;
  };
  using Element = dd4hep::DetElement;

  struct Config {
    const dd4hep::Detector* detector = nullptr;
    double lengthScale = 1.0;
    std::reference_wrapper<const Acts::GeometryContext> gctx;
    BlueprintBuilder::ElementFactory elementFactory
        = BlueprintBuilder::defaultElementFactory;
  };

  explicit DD4hepBackend(Config cfg);

  // Concept interface
  Element world() const;
  std::string_view nameOf(const Element& e) const;
  std::optional<Element> findByName(std::string_view name) const;
  std::vector<Element> findByPattern(const Element& parent,
                                     const std::regex& pat) const;
  std::vector<Element> children(const Element& parent) const;
  bool isBarrel(const Element& e) const;
  bool isEndcap(const Element& e) const;
  bool isSensitive(const Element& e) const;

  // DD4hep-specific (outside concept)
  std::shared_ptr<Acts::Experimental::LayerBlueprintNode>
  makeLayer(const Element& parent,
            std::span<const Element> sensitives,
            const LayerSpec& spec,
            std::optional<std::string> name = std::nullopt) const;

  std::shared_ptr<Acts::Experimental::StaticBlueprintNode>
  makeBeampipe() const;

 private:
  Config m_cfg;
};

// Suppress implicit instantiation — definitions live in BlueprintBuilder.cpp
extern template class LayerAssembler<DD4hepBackend>;
extern template class BarrelEndcapAssembler<DD4hepBackend>;

// Callers keep the existing unqualified names
using LayerAssembler        = ::ActsPlugins::DD4hep::LayerAssembler<DD4hepBackend>;
using BarrelEndcapAssembler = ::ActsPlugins::DD4hep::BarrelEndcapAssembler<DD4hepBackend>;
```

### BlueprintBuilder — thin wrapper

After the split, `BlueprintBuilder` retains only:
- `Config` (same fields, delegates to `DD4hepBackend::Config`)
- Logger ownership
- Forwarding methods `layers()` / `barrelEndcap()` that inject the logger

All search/classification/construction logic moves into `DD4hepBackend`.

### Explicit instantiation in BlueprintBuilder.cpp

```cpp
// BlueprintBuilder.cpp

// Pull in template bodies (NOT a public header)
#include "GeometryBackend_impl.hpp"

// Instantiate for DD4hep — one TU, link once
template class ActsPlugins::DD4hep::LayerAssembler<DD4hepBackend>;
template class ActsPlugins::DD4hep::BarrelEndcapAssembler<DD4hepBackend>;
```

### TGeoAxes → GeometryAxes

`TGeoAxes` has zero external dependencies (std lib only). Steps:
1. Copy to a shared header (e.g. `Acts/Geometry/GeometryAxes.hpp` or
   `ActsPlugins/GeometryAxes.hpp`), rename class to `GeometryAxes`.
2. In `ActsPlugins/Root/TGeoAxes.hpp` add `using TGeoAxes = GeometryAxes;`
   (backwards-compatible alias).
3. DD4hep plugin uses `GeometryAxes` natively; TGeo plugin keeps `TGeoAxes`
   as-is via the alias.

### Migration path

1. Lift `TGeoAxes` → `GeometryAxes` (shared header, alias in Root plugin)
2. Write `DD4hepBackend` struct + concept check
3. Template `LayerAssembler` / `BarrelEndcapAssembler` over `GeometryBackend`
4. Add `_impl.hpp`, explicit instantiation in `.cpp`
5. Make `BlueprintBuilder` a thin wrapper delegating to `DD4hepBackend`
6. Update `using` aliases in `DD4hepBackend.hpp` — callers unchanged
7. (Later) Add `TGeoBackend` and/or `GeoModelBackend` following same pattern
