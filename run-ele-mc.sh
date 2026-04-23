#!/usr/bin/env bash

set -euo pipefail

repo_dir="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${repo_dir}"

make Execute_triggerAnalysis_ele_mc
mkdir -p plots output

ztoee_input="/eos/cms/store/group/phys_heavyions/kdeverea/Run3_PbPb_2025MC/JpsiDielectron_pTHatMin4_HydjetEmbedded_Pythia8_TuneCP5_1510pre6/crab_Run3_PbPb_2025MC_ZToEE/251115_102359/0000/"
jpsitoee_input="/eos/cms/store/group/phys_heavyions/kdeverea/Run3_PbPb_2025MC/JpsiDielectron_pTHatMin4_HydjetEmbedded_Pythia8_TuneCP5_1510pre6/crab_Run3_PbPb_2025MC_JpsiToEE/251113_150758/0000/"

run_case() {
  local input_forest="$1"
  local output_base="$2"
  local nametag="$3"
  local min_hi_bin="$4"
  local max_hi_bin="$5"
  local match_to_self="${6:-true}"

  echo "Running ${output_base} (hiBin ${min_hi_bin}-${max_hi_bin})"
  ./Execute_triggerAnalysis_ele_mc \
    --inputForest "${input_forest}" \
    --outputBase "${output_base}" \
    --nametag "${nametag}" \
    --matchToSelf "${match_to_self}" \
    --nfiles -1 \
    --minHiBin "${min_hi_bin}" \
    --maxHiBin "${max_hi_bin}"
}

run_case "${ztoee_input}"   "MCZee0_100_2023PbPbcuts"      "2025 Pythia8+Hydjet Z->EE"      0   200     "true"
run_case "${ztoee_input}"   "MCZee0_30_2023PbPbcuts"       "2025 Pythia8+Hydjet Z->EE"      0    60     "true"
run_case "${ztoee_input}"   "MCZee30_100_2023PbPbcuts"     "2025 Pythia8+Hydjet Z->EE"     60   200     "true"

run_case "${jpsitoee_input}" "MCJpsiToEE0_100_2023PbPbcuts"  "2025 Pythia8+Hydjet J/psi->EE"  0   200     "false"
run_case "${jpsitoee_input}" "MCJpsiToEE0_30_2023PbPbcuts"   "2025 Pythia8+Hydjet J/psi->EE"  0    60     "false"
run_case "${jpsitoee_input}" "MCJpsiToEE30_100_2023PbPbcuts" "2025 Pythia8+Hydjet J/psi->EE" 60   200     "false"
