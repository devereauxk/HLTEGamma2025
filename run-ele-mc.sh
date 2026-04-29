#!/usr/bin/env bash

set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${repo_dir}"

make Execute_triggerAnalysis_ele_mc Execute_triggerRatio_ele_mc
mkdir -p plots output

ztoee_input="/eos/cms/store/group/phys_heavyions/kdeverea/Run3_PbPb_2025MC/JpsiDielectron_pTHatMin4_HydjetEmbedded_Pythia8_TuneCP5_1510pre6/crab_Run3_PbPb_2025MC_ZToEE/251115_102359/0000/"
official_emulation="/eos/cms/store/group/phys_heavyions/kdeverea/Run3_PbPb_2025MC/JpsiDielectron_pTHatMin4_HydjetEmbedded_Pythia8_TuneCP5_1510pre6/crab_Run3_PbPb_2025MC_ZToEE/emulation/dev_CMSSW16_0_0_HIon_V56/"
HOverE0p5_emulation="/eos/cms/store/group/phys_heavyions/kdeverea/Run3_PbPb_2025MC/JpsiDielectron_pTHatMin4_HydjetEmbedded_Pythia8_TuneCP5_1510pre6/crab_Run3_PbPb_2025MC_ZToEE/emulation/dev_CMSSW16_0_0_HIon_V56_HOverE0p5/"

run_case() {
  local input_forest="$1"
  local input_hlt="$2"
  local output_base="$3"
  local nametag="$4"
  local min_hi_bin="$5"
  local max_hi_bin="$6"
  local match_to_self="${7:-true}"
  local match_to_emulation="${8:-false}"

  echo "Running ${output_base} (hiBin ${min_hi_bin}-${max_hi_bin})"
  ./Execute_triggerAnalysis_ele_mc \
    --inputForest "${input_forest}" \
    --inputHLT "${input_hlt}" \
    --outputBase "${output_base}" \
    --nametag "${nametag}" \
    --matchToSelf "${match_to_self}" \
    --matchToEmulation "${match_to_emulation}" \
    --nfiles -1 \
    --minHiBin "${min_hi_bin}" \
    --maxHiBin "${max_hi_bin}"
}

run_ratio_case() {
  local input_numerator="$1"
  local input_denominator="$2"
  local output_base="$3"
  local nametag="$4"
  local ratio_label="$5"
  local min_hi_bin="$6"
  local max_hi_bin="$7"

  echo "Running ratio ${output_base} (hiBin ${min_hi_bin}-${max_hi_bin})"
  ./Execute_triggerRatio_ele_mc \
    --inputNumerator "${input_numerator}" \
    --inputDenominator "${input_denominator}" \
    --outputBase "${output_base}" \
    --nametag "${nametag}" \
    --ratioLabel "${ratio_label}" \
    --minHiBin "${min_hi_bin}" \
    --maxHiBin "${max_hi_bin}"
}

#run_case "${ztoee_input}"    ""                           "MCZee0_100_2023PbPbcuts"                 "2025 Pythia8+Hydjet Z->EE"                               0   200  "true"  "false"
#run_case "${ztoee_input}"    ""                           "MCZee0_30_2023PbPbcuts"                  "2025 Pythia8+Hydjet Z->EE"                               0    60  "true"  "false"
#run_case "${ztoee_input}"    ""                           "MCZee30_100_2023PbPbcuts"                "2025 Pythia8+Hydjet Z->EE"                              60   200  "true"  "false"

run_case "${ztoee_input}" "${official_emulation}" "MCZee0_100_officialEmulation" "2025 Pythia8+Hydjet Z->EE" 0 200 "false" "true"
run_case "${ztoee_input}" "${official_emulation}" "MCZee0_30_officialEmulation" "2025 Pythia8+Hydjet Z->EE" 0 60 "false" "true"
run_case "${ztoee_input}" "${official_emulation}" "MCZee30_100_officialEmulation" "2025 Pythia8+Hydjet Z->EE" 60 200 "false" "true"

run_case "${ztoee_input}" "${HOverE0p5_emulation}" "MCZee0_100_HOverE0p5" "2025 Pythia8+Hydjet Z->EE" 0 200 "false" "true"
run_case "${ztoee_input}" "${HOverE0p5_emulation}" "MCZee0_30_HOverE0p5" "2025 Pythia8+Hydjet Z->EE" 0 60 "false" "true"
run_case "${ztoee_input}" "${HOverE0p5_emulation}" "MCZee30_100_HOverE0p5" "2025 Pythia8+Hydjet Z->EE" 60 200 "false" "true"

run_ratio_case \
  "output/output_ele_MCZee0_100_HOverE0p5.root" \
  "output/output_ele_MCZee0_100_officialEmulation.root" \
  "MCZee0_100_HOverE0p5_over_officialEmulation" \
  "2025 Pythia8+Hydjet Z->EE" \
  "HOverE0p5 emulation / official emulation" \
  0 \
  200

run_ratio_case \
  "output/output_ele_MCZee0_30_HOverE0p5.root" \
  "output/output_ele_MCZee0_30_officialEmulation.root" \
  "MCZee0_30_HOverE0p5_over_officialEmulation" \
  "2025 Pythia8+Hydjet Z->EE" \
  "HOverE0p5 emulation / official emulation" \
  0 \
  60

run_ratio_case \
  "output/output_ele_MCZee30_100_HOverE0p5.root" \
  "output/output_ele_MCZee30_100_officialEmulation.root" \
  "MCZee30_100_HOverE0p5_over_officialEmulation" \
  "2025 Pythia8+Hydjet Z->EE" \
  "HOverE0p5 emulation / official emulation" \
  60 \
  200
