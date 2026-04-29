# 2026-04-23 Zee EOS production status

This note records the current state of the 2025 Zee DIGIRAW HLT re-emulation outputs written to EOS.

## EOS destinations

- fdamas menu:
  - `/eos/cms/store/group/phys_heavyions/kdeverea/Reemulation/2025PbPb_ZtoEE`
- thrOverEEE menu:
  - `/eos/cms/store/group/phys_heavyions/kdeverea/Reemulation/2025PbPb_ZtoEE_thrOverEEE_0p5`

## Menus used

- fdamas:
  - `/users/fdamas/2025/HIon/151X/HLT/V5`
- thrOverEEE:
  - `/users/kdeverea/HIon_2025/fdamas_151X_V5_thrOverEEE_0p5/V2`

Both were run from `CMSSW_15_1_0` with:

- `L1Menu_CollisionsHeavyIons2025_v1_0_3.xml`
- `Run3_pp_on_PbPb_2025`
- `HLT_FORCE_HI_TRACK_DNN=1`

## Reduced validation that passed

Before launching larger productions, reduced 100-event checks were validated for both menus.

Validated reduced outputs:

- `fdamas-v5-0.root`
- `fdamas-v5-1.root`
- `fdamas-v5-thrOverEEE0p5-0.root`
- `fdamas-v5-thrOverEEE0p5-1.root`

Each of those contained nontrivial trigger content, including:

- `hltanalysis/HltTree`
- `hltobject/HLT_HIEle20Gsf_v`
- `hltobject/HLT_HIDoubleEle10Gsf_v`

with 100 entries in the reduced tests.

## Full-production coverage audit

An intermediate audit found only 67 out of 83 expected index-labeled ROOT files in each EOS destination.

The matching missing-index pattern at that stage was:

- `7 8 16 17 25 26 34 35 43 44 52 53 61 62 70 71`

This suggested that the remaining work should be treated as gap-filling rather than rerunning the whole sample blindly.

## Gap-filling validation

The first missing index (`7`) was revalidated with a 100-event run for both menus in isolated work directories.

Temporary validation outputs:

- `/tmp/gapcheck-fdamas-7.root`
- `/tmp/gapcheck-thr-7.root`

Validated content:

- `gapcheck-fdamas-7`: `HltTree=100`, `HLT_HIEle20Gsf_v=100`, `HLT_HIDoubleEle10Gsf_v=100`
- `gapcheck-thr-7`: `HltTree=100`, `HLT_HIEle20Gsf_v=100`, `HLT_HIDoubleEle10Gsf_v=100`

So the missing index itself is not intrinsically broken under a reduced test.

## First full gap fill validation

The first missing pair (`7..8`) was then written to the final EOS destinations for both menus.

Validated new full outputs:

- `fdamas-v5-8.root`
- `fdamas-v5-thrOverEEE0p5-8.root`

Validated content:

- `fdamas-v5-8.root`:
  - `hltanalysis/HltTree = 2101`
  - `hltobject/HLT_HIEle20Gsf_v = 2101`
  - `hltobject/HLT_HIDoubleEle10Gsf_v = 2101`
- `fdamas-v5-thrOverEEE0p5-8.root`:
  - `hltanalysis/HltTree = 2101`
  - `hltobject/HLT_HIEle20Gsf_v = 2101`
  - `hltobject/HLT_HIDoubleEle10Gsf_v = 2101`

The corresponding `7.root` files were also present on EOS with nontrivial file sizes.

## Notes on prior detached attempts

Separate isolated work-directory launches were used because `emulate_HLT_HIon_2025.sh` writes fixed filenames such as:

- `test_pset.py`
- `output_130.root`
- `openHLT.root`
- `logs/log_i.txt`
- `myGets/myGet_i.txt`

so concurrent runs in the same directory are not safe.

Earlier detached launches from copied work directories did not provide a clean end-to-end full-fill solution for the missing coverage; one observed failure mode was a crash in:

- `alpaka_serial_sync::CAHitNtupletAlpakaHIonPhase1:hltPixelTracksPPOnAASoA`

Because of that, the current strategy is to validate and fill the missing EOS indices explicitly.

## Final EOS state

After the later updates visible in the shared environment, both EOS destinations reached full index coverage:

- fdamas:
  - `83 / 83` files present for indices `0..82`
- thrOverEEE:
  - `83 / 83` files present for indices `0..82`

Additional spot checks beyond the original reduced tests:

- `fdamas-v5-82.root`:
  - `hltanalysis/HltTree = 2166`
  - `hltobject/HLT_HIEle20Gsf_v = 2166`
  - `hltobject/HLT_HIDoubleEle10Gsf_v = 2166`
- `fdamas-v5-thrOverEEE0p5-82.root`:
  - `hltanalysis/HltTree = 2166`
  - `hltobject/HLT_HIEle20Gsf_v = 2166`
  - `hltobject/HLT_HIDoubleEle10Gsf_v = 2166`
- `fdamas-v5-71.root`:
  - `hltanalysis/HltTree = 2083`
  - `hltobject/HLT_HIEle20Gsf_v = 2083`
  - `hltobject/HLT_HIDoubleEle10Gsf_v = 2083`
- `fdamas-v5-thrOverEEE0p5-71.root`:
  - `hltanalysis/HltTree = 2083`
  - `hltobject/HLT_HIEle20Gsf_v = 2083`
  - `hltobject/HLT_HIDoubleEle10Gsf_v = 2083`

So the final EOS state is complete and nontrivial for both menus.
