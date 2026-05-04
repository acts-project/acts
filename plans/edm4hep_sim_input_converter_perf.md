# EDM4hepSimInputConverter — performance ROI plan

Target: [Examples/Io/EDM4hep/src/EDM4hepSimInputConverter.cpp](../Examples/Io/EDM4hep/src/EDM4hepSimInputConverter.cpp)

Validation command:
```bash
./run.sh $(which python3) Examples/Scripts/Python/full_chain_odd.py \
    -n 13 --output-parquet \
    --edm4hep '/Users/pagessin/Downloads/edm4hep.root' -j1
```
The chain prints per-step timing at the end; this converter is one of the early steps and its wall time is the metric to track. Re-run after each change and record the delta.

Ranking is by **(expected speedup or memory drop) ÷ (implementation effort + risk)**. "Effort" is rough hours-of-work; "risk" reflects how much surrounding code has to be reasoned about.

---

## Tier 1 — high ROI, attack first

### 1. Memoize "subtree has hits"
- **Where:** [`particleOrDescendantsHaveHits`](../Examples/Io/EDM4hep/src/EDM4hepSimInputConverter.cpp#L172-L189) called from [`processChildren`](../Examples/Io/EDM4hep/src/EDM4hepSimInputConverter.cpp#L949-L950).
- **Problem:** Walks the entire descendant subtree once per daughter visited by `processChildren`. Worst case Θ(N·depth); for typical Geant gamma/π⁰ chains this is the dominant cost.
- **Fix:** Single bottom-up pass over `mcParticleCollection` that fills a `std::vector<std::uint16_t> subtreeHits` (sum of own + children's `numSimHits`). Predicate becomes `subtreeHits[idx] > 0`, O(1).
- **Bonus:** Fold the pass into the existing `numSimHits` counting loop — no extra allocation.
- **Effort:** ~1 h. **Risk:** low (pure replacement of a predicate).
- **Expected gain:** large on busy events; should be visible in the converter's wall time.

### 2. Drop per-hit `find` on `SimParticleContainer`
- **Where:** [lines 656-666](../Examples/Io/EDM4hep/src/EDM4hepSimInputConverter.cpp#L656-L666) in the simhit loop.
- **Problem:** `flat_set::find` is O(log P) per hit (cache-friendly but still ~17 comparisons over a non-contiguous access pattern). H·log(P) work for something we already know.
- **Fix:** Use `numSimHits` (already built). After constructing `particlesGenerator`/`particlesSimulated`, walk each container once and set `numberOfHits` directly. Either invert `edm4hepParticleMap` to index from final-container position, or index `numSimHits` by `unorderedParticlesInitial` slot instead of edm4hep index.
- **Effort:** ~1–2 h (touch the bookkeeping in two places). **Risk:** low–medium (need to confirm `numberOfHits` is not part of the flat_set ordering key — it shouldn't be, ordering is by `particleId()`).
- **Expected gain:** removes a hot per-hit operation entirely; speedup scales with hit count.

### 3. Replace `multimap` in time-sort with sorted vector
- **Where:** [lines 743-790](../Examples/Io/EDM4hep/src/EDM4hepSimInputConverter.cpp#L743-L790).
- **Problem:** `std::multimap` does one node allocation per hit (10⁶ allocations on a busy event, fragmented heap, big peak memory). Also `count()`+`equal_range()`+`upper_bound()` redundantly walks the tree.
- **Fix:** `std::vector<std::pair<SimBarcode, std::size_t>>`, sort once, group with adjacent iterators.
- **Effort:** ~1 h. **Risk:** low.
- **Expected gain:** wall-time and a meaningful drop in peak RSS.

### 4. Reserve `unorderedParticlesInitial` up front
- **Where:** [line 233](../Examples/Io/EDM4hep/src/EDM4hepSimInputConverter.cpp#L233).
- **Problem:** Geometric reallocation; transient peak ~2× the final size.
- **Fix:** `unorderedParticlesInitial.reserve(mcParticleCollection.size())` immediately after the count is known.
- **Effort:** 5 min. **Risk:** none.
- **Expected gain:** small but free; helps high-water mark.

---

## Tier 2 — solid wins, slightly more work

### 5. Hash-bucket `maybeAddVertex` by `SimVertexBarcode`
- **Where:** [lines 803-847](../Examples/Io/EDM4hep/src/EDM4hepSimInputConverter.cpp#L803-L847).
- **Problem:** Linear scan over all vertices for every simulated particle with hits — O(N²).
- **Fix:** `std::unordered_map<SimVertexBarcode, std::vector<std::size_t>>` keyed on id; only scan the (typically tiny) bucket. Also gate `getMinDistance()` on `logger().doPrint(VERBOSE)` — currently the function call itself runs even when the message is suppressed.
- **Effort:** ~1.5 h. **Risk:** low (need a hash for `SimVertexBarcode` — likely already exists or trivial).
- **Expected gain:** scales with particles-with-hits; on pile-up events this can be substantial.

### 6. De-`std::function` the hot callbacks
- **Where:** `getNumHits` ([line 335](../Examples/Io/EDM4hep/src/EDM4hepSimInputConverter.cpp#L335)), `geometryMapper` / `particleMapper` ([lines 560, 594](../Examples/Io/EDM4hep/src/EDM4hepSimInputConverter.cpp#L560)).
- **Problem:** `std::function` blocks inlining and forces an indirect call per hit / per descendant.
- **Fix:** Template `processChildren` and `EDM4hepUtil::readSimHit` on the callable type, pass lambdas by `auto&&`. While there, fix the `uint8_t` (line 174) vs `uint16_t` (line 335) truncation — currently silently wraps for particles with ≥256 hits.
- **Effort:** ~2 h (header churn). **Risk:** medium (signature change ripples to header/util).
- **Expected gain:** modest individually, but compounds with other hot-loop fixes.

### 7. Hit-association: prefix-sum the collection sizes
- **Where:** [lines 713-731](../Examples/Io/EDM4hep/src/EDM4hepSimInputConverter.cpp#L713-L731).
- **Problem:** For every output hit, walks `simHitCollections` to find which input collection a global index belongs to. O(H·C).
- **Fix:** Precompute `std::vector<size_t> prefix` of cumulative sizes; binary-search with `std::ranges::upper_bound`. O(H·log C).
- **Effort:** ~30 min. **Risk:** low. Only runs when `m_outputSimHitAssociation` is initialized.
- **Expected gain:** small unless C is large; cheap enough to do anyway.

---

## Tier 3 — nice-to-have, low priority

### 8. Primary-vertex grouping uses a hash
- **Where:** [lines 257-259](../Examples/Io/EDM4hep/src/EDM4hepSimInputConverter.cpp#L257-L259).
- **Problem:** Linear scan with exact `Vector3` equality. Fragile and O(V·P) per parentless particle. For typical luminosities V is small so this is rarely the bottleneck — but exact float equality is a latent bug.
- **Fix:** If EDM4hep exposes a parent-vertex index, key on that. Otherwise quantize position to a hash key.
- **Effort:** 30 min–1 h depending on EDM4hep API. **Risk:** low.

### 9. `findGeneratorStableParticles` deduplication
- **Where:** [line 211](../Examples/Io/EDM4hep/src/EDM4hepSimInputConverter.cpp#L211).
- **Problem:** O(k²) linear-scan dedup per primary vertex.
- **Fix:** A `std::vector<bool>`-or-`uint8_t` `seen` array sized to `mcParticleCollection.size()`, allocated once outside the PV loop.
- **Effort:** 20 min. **Risk:** none.
- **Expected gain:** minor for low-multiplicity events; matters for pile-up.

### 10. Loop-fuse the generator-particle counter
- **Where:** [lines 301-306](../Examples/Io/EDM4hep/src/EDM4hepSimInputConverter.cpp#L301-L306).
- **Problem:** Whole-collection scan whose only purpose is one counter; can fold into the simhit-counting loop.
- **Effort:** 5 min. **Risk:** none.
- **Expected gain:** trivial; do it for tidiness while editing nearby.

### 11. Drop unused `position` in surface-map setup
- **Where:** [lines 113-119](../Examples/Io/EDM4hep/src/EDM4hepSimInputConverter.cpp#L113-L119).
- One-time cost so irrelevant to runtime, but it's dead code.

---

## Suggested execution order

Do them in this order, profiling between each (the validation command above will tell you):

1. **#1 + #4 + #10** as one commit — same area of code, all in the early counting/walk passes.
2. **#2** — independent, drops a hot per-hit operation.
3. **#3** — independent, biggest peak-memory win.
4. **#5** — wall-time on busy events.
5. **#6** — once everything else is settled, since it touches headers.
6. **#7, #8, #9, #11** — opportunistic cleanup.

Stop when the converter falls out of the top-N steps in the timing breakdown. Don't pre-emptively do tier 3 if tier 1+2 already make the converter cheap relative to fitting/finding.

---

## Things explicitly *not* on this list

- `std::vector<bool>` for the subtree-has-hits memo. Bit-packing saves ~12 KB at this scale, costs branchier codegen, and the proxy-reference semantics are awkward. Use `std::vector<std::uint16_t>` (folded into `numSimHits`) or `std::vector<std::uint8_t>`.
- Parallelizing the per-event passes. The chain is `-j1` and event-level parallelism already exists upstream; intra-event threading would fight the framework's expectations.
- Caching the `dd4hepDetector` volume manager lookups. The data is already a hash; layering another cache on top rarely helps and adds invalidation surface.
