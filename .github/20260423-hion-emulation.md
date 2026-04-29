# 2026-04-23 HIon HLT emulation status

This note records the current status of the `HLTClayton` Heavy-Ion 2025 HLT emulation work as validated from the `HLTEgamma2025` task flow.

## Working menu

- Current working point: `/users/fdamas/2025/HIon/151X/HLT/V5`
- Release used: `CMSSW_15_1_0`
- L1 XML used: `L1Menu_CollisionsHeavyIons2025_v1_0_3.xml`
- Era: `Run3_pp_on_PbPb_2025`
- Input used for validation: 2025 Zee DIGIRAW sample

This menu now runs in the **full untrimmed configuration** after replacing the HI PPOnAA tracking `GBRWrapperRcd`-based classifiers with the TF/DNN classifier already present in the generated menu.

## Fixes that were validated

### 1. Replace missing HI tracking GBR payload users with TF/DNN classifiers

Problem:

- `TrackMVAClassifierPrompt/hltFullIter0TrackMVAClassifierPPOnAA`
- `TrackMVAClassifierPrompt/hltFullIter1TrackMVAClassifierPPOnAA`
- `TrackMVAClassifierPrompt/hltFullIter2TrackMVAClassifierPPOnAA`

were requesting:

- `HIMVASelectorInitialStep_Phase1`
- `HIMVASelectorLowPtQuadStep_Phase1`
- `HIMVASelectorHighPtTripletStep_Phase1`

from `GBRWrapperRcd`, and the job failed with:

- `NoProductResolverException` for missing `GBRForest`

Validated replacement:

- use `TrackTfClassifier`
- use the existing `hltESPTrackSelectionTfCKF`
- keep PPOnAA step-specific quality cuts

This replacement is injected through `setup_hltConfig.sh` with:

- `HLT_FORCE_HI_TRACK_DNN=1`

### 2. Clamp invalid HIon pixel morphing config in dev V9

Problem:

- `/dev/CMSSW_15_1_0/HIon/V9` configured
  - `hltSiPixelClustersPPOnAASoA`
  - `hltSiPixelClustersPPOnAASoASerialSync`

with:

- `MaxFakesInModule = 4000`

but `RecoLocalTracker/SiPixelClusterizer/plugins/alpaka/SiPixelRawToCluster.cc` in `CMSSW_15_1_0` enforces:

- `TrackerTraits::maxPixInModuleForMorphing <= 1000`

Validated workaround:

- set `MaxFakesInModule = 1000`

This is injected through `setup_hltConfig.sh` with:

- `HLT_CAP_PIXEL_MAX_FAKES=1`

## Menu-by-menu status

### `/users/fdamas/2025/HIon/151X/HLT/V5`

Status:

- **works in full untrimmed mode** with `HLT_FORCE_HI_TRACK_DNN=1`

Validated output:

- `hltanalysis/HltTree`
- `hltobject/HLT_HIEle20Gsf_v`
- `hltobject/HLT_HIDoubleEle10Gsf_v`

Validated 100-event rates included:

- `HLT_HIEle10Gsf_v16 = 89`
- `HLT_HIEle20Gsf_v16 = 86`
- `HLT_HIDoubleEle10Gsf_v16 = 74`
- `HLT_HIL1SingleMu5_SingleEG20Gsf_v9 = 1`

### `/dev/CMSSW_15_1_0/HIon/V9`

Status:

- **not usable as a full menu** with the current L1 XML

After fixing the old configuration/runtime blockers, the full menu fails first in:

- path `HLT_HIUPC_HFafterglowCombined_v2`
- module `HLTL1TSeed/hltL1sUPCHFafterglowCombined`

because the menu requests L1 algorithms such as:

- `L1_SingleJet8_ZDC1n_XOR_NotPreBptx_BptxAND`

which are not present in:

- `L1Menu_CollisionsHeavyIons2025_v1_0_3.xml`

Also validated:

- the **electron-only untrimmed** dev V9 still segfaults in
  - `alpaka_serial_sync::CAHitNtupletAlpakaHIonPhase1:hltPixelTracksPPOnAASoA`

- the **electron-only trimmed** fallback works end-to-end with:
  - `HLT_FORCE_HI_TRACK_DNN=1`
  - `HLT_CAP_PIXEL_MAX_FAKES=1`
  - `HLT_TRIM_ELECTRON_TRACKING=1`
  - a path filter keeping only electron paths

## Practical run settings

For the working full fdamas menu, the key setting is:

- `HLT_FORCE_HI_TRACK_DNN=1`

For dev V9 studies, also use:

- `HLT_CAP_PIXEL_MAX_FAKES=1`

For dev V9 electron-only fallback studies, additionally use:

- `HLT_TRIM_ELECTRON_TRACKING=1`
- `HLT_KEEP_PATHS_REGEX='HLTrigger(First|Final)Path|^HLT_HI(?!.*Mass50)(Ele|DoubleEle|L1SingleMu5_SingleEG)'`

## Zee EOS production status

- Dataset size: 83 DIGIRAW files
- EOS target for the working fdamas menu:
  - `/eos/cms/store/group/phys_heavyions/kdeverea/Reemulation/2025PbPb_ZtoEE`
- EOS target for the thrOverEEE menu:
  - `/eos/cms/store/group/phys_heavyions/kdeverea/Reemulation/2025PbPb_ZtoEE_thrOverEEE_0p5`

Reduced validation completed before launching the full productions:

- fdamas output prefix:
  - `fdamas-v5-`
- thrOverEEE output prefix:
  - `fdamas-v5-thrOverEEE0p5-`
- validated files:
  - `fdamas-v5-0.root`
  - `fdamas-v5-1.root`
  - `fdamas-v5-thrOverEEE0p5-0.root`
  - `fdamas-v5-thrOverEEE0p5-1.root`
- validated content for each reduced output:
  - `hltanalysis/HltTree` with 100 entries
  - `hltobject/HLT_HIEle20Gsf_v` with 100 entries
  - `hltobject/HLT_HIDoubleEle10Gsf_v` with 100 entries

Full detached productions were then launched from isolated work directories to avoid collisions on
shared files such as `test_pset.py`, `openHLT.root`, `logs/log_i.txt`, and `myGets/myGet_i.txt`:

- fdamas workdir:
  - `~/HLTClayton/CMSSW_15_1_0/src/HLTrigger/Configuration/test/workstation/HLT_emulation/scripts/PbPb_fdamas_prod`
- thrOverEEE workdir:
  - `~/HLTClayton/CMSSW_15_1_0/src/HLTrigger/Configuration/test/workstation/HLT_emulation/scripts/PbPb_thrOverEEE_prod`
- detached driver logs:
  - `/tmp/fdamas-full-zee.driver.log`
  - `/tmp/throvereee-full-zee.driver.log`

## Remaining known issues

1. No validated replacement L1 XML has yet been found locally for the missing dev V9 UPC/ZDC/jet seed names.
2. The untrimmed dev V9 electron path still crashes in `hltPixelTracksPPOnAASoA`.
3. The long Zee productions above were launched successfully, but they were still running at the time this note was updated.
