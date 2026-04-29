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
- `triggerRatio_ele_mc.cpp` compares two electron `triggerAnalysis` ROOT outputs and plots menu-to-menu efficiency ratios.
- Build it with:

```bash
make Execute_triggerAnalysis_ele_mc
```

- Build the ratio plotter with:

```bash
make Execute_triggerRatio_ele_mc
```

- Run the compiled electron MC executable with flags like:

```bash
./Execute_triggerAnalysis_ele_mc \
  --inputForest <forest-dir> \
  --inputHLT <emulation-dir or empty> \
  --outputBase <tag> \
  --nametag <plot-label> \
  --matchToSelf <true|false> \
  --matchToEmulation <true|false> \
  --nfiles <N or -1> \
  --minHiBin <min> \
  --maxHiBin <max>
```

- Batch production for the 2025 ZToEE and JpsiToEE electron MC plots is handled by:

```bash
./run-ele-mc.sh
```

- The current `run-ele-mc.sh` also drives the official-emulation vs `HOverE0p5` comparison workflow and writes ratio plots named like `TriggerEfficiencyRatio_ele_<tag>_{barrel,endcap}.png`.

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
- The 2025 ZToEE reemulation files in `Reemulation/2025PbPb_ZtoEE` and `Reemulation/2025PbPb_ZtoEE_thrOverEEE_0p5` store trigger decisions in `hltanalysis/HltTree` with unversioned branches like `HLT_HIEle20Gsf_v`, plus matching `hltobject/HLT_HIEle{15,20,30,40,50}Gsf_v` trees.
- In the current reduced 2025 ZToEE reemulation outputs, the studied HLT object trees and HLT decision branches are populated, but `L1_SingleEG{7,15,21}_BptxAND` are empty. For this workflow, keep the forest L1 denominator and only replace the HLT decisions and object collections in emulation mode.

## ROOT/TChain gotchas learned here

- Do not assume `hltobject` trees are present in every forest sample just because `hltanalysis/HltTree` exists.
- For the current `triggerAnalysis_ele_mc.C`, reading HLT object trees through a multi-file `TChain` was unreliable across file boundaries. The working approach is to:
  - use the main `HltTree` chain to detect file changes
  - fetch the `hltobject/...` trees from the current underlying file
  - read them with the local entry number
- If `matchToSelf=false`, do not dereference HLT object vectors. In that mode, the code should fall back to pure HLT decision bits.
- If `matchToEmulation=true`, do not use same-file `hltobject` trees. Instead key the emulation lookup on `(Run, LumiBlock, Event)` from `hltanalysis/HltTree`, use the emulation HLT decisions plus `hltobject/...` collections from the separate ROOT files, keep the forest L1 denominator, and skip forest events that have no emulation entry.

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
