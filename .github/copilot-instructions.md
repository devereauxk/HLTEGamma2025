# HLTEgamma2025 repository guidance

## Repository summary

This repository contains ROOT-based trigger performance studies for heavy-ion and pp reference EGamma analyses. The main workflows are:

- forest production and configuration notes in `README.md` and `ForestSubmission/`
- standalone ROOT macros such as `triggerAnalysis_ele_mc.C`, `triggerAnalysis_ele_PbPb.C`, `triggerAnalysis_pho_mc.C`
- compiled C++ studies such as `HLTperformance_*.cpp` and `HLTrate_ppref.cpp`
- generated plots in `plots/` and ROOT histogram outputs in `output/`

The code is analysis-oriented rather than library-oriented: many scripts hardcode default EOS inputs and produce plots directly.

## Build and run conventions

- Existing compiled tools are built with plain `g++` plus ``root-config --cflags --libs`` and `-lASImage -std=c++17`.
- The top-level `makefile` is intentionally simple and is the right place to add new executable targets.
- `triggerAnalysis_ele_mc.cpp` is a compiled wrapper around `triggerAnalysis_ele_mc.C`.
- Build it with:

```bash
make Execute_triggerAnalysis_ele_mc
```

- Run the compiled electron MC executable with flags like:

```bash
./Execute_triggerAnalysis_ele_mc \
  --inputForest <forest-dir> \
  --outputBase <tag> \
  --nametag <plot-label> \
  --matchToSelf <true|false> \
  --nfiles <N or -1> \
  --minHiBin <min> \
  --maxHiBin <max>
```

- Batch production for the 2025 ZToEE and JpsiToEE electron MC plots is handled by:

```bash
./run-ele-mc.sh
```

## Important analysis conventions

- `hiBin = 2 * centrality percentage`. For example:
  - `0-100%` -> `minHiBin=0`, `maxHiBin=200`
  - `0-30%` -> `minHiBin=0`, `maxHiBin=60`
  - `30-100%` -> `minHiBin=60`, `maxHiBin=200`
- Trigger plots are expected to be written to `plots/`.
- Histograms are expected to be written to `output/`.
- For the electron MC workflow, the current intended efficiency definition is:
  - select reco electrons passing the offline cuts
  - require the appropriate L1 seed for the denominator
  - if reco-to-HLT matching is enabled, match reco electrons to HLT objects by minimum `deltaR`
  - numerator = reco electron matched to HLT object **and** HLT path fired

## Known sample-specific facts

- The 2025 **ZToEE** forest files in the default EOS path include `hltobject/HLT_HIEle{15,20,30,40,50}Gsf_v` trees, so `matchToSelf=true` works there.
- The 2025 **JpsiToEE** forest files in the default EOS path do **not** include those `hltobject/...` trees. They still contain `hltanalysis/HltTree` and trigger decision branches, but not HLT object collections.
- Because of that, JpsiToEE jobs must run with `matchToSelf=false` unless a different input source with HLT object trees is provided.

## ROOT/TChain gotchas learned here

- Do not assume `hltobject` trees are present in every forest sample just because `hltanalysis/HltTree` exists.
- For the current `triggerAnalysis_ele_mc.C`, reading HLT object trees through a multi-file `TChain` was unreliable across file boundaries. The working approach is to:
  - use the main `HltTree` chain to detect file changes
  - fetch the `hltobject/...` trees from the current underlying file
  - read them with the local entry number
- If `matchToSelf=false`, do not dereference HLT object vectors. In that mode, the code should fall back to pure HLT decision bits.

## Code patterns to preserve

- Keep changes surgical and analysis-focused; many scripts are used as one-off studies with specific default samples.
- Prefer updating labels, defaults, and output naming to match the actual analysis logic when behavior changes.
- When adding new batch workflows, preserve the pattern of:
  - executable target in `makefile`
  - EOS defaults in the script
  - explicit output tags per sample and centrality bin

## Files future agents should read first for electron MC work

- `README.md`
- `triggerAnalysis_ele_mc.C`
- `triggerAnalysis_ele_mc.cpp`
- `run-ele-mc.sh`
- `makefile`
