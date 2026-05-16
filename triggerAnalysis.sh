#!/usr/bin/env bash
set -euo pipefail

usage() {
  echo "Usage: $0 mc <year> <particle> <centMin> <centMax> <sample>"
  echo ""
  echo "  year     : 2025 or 2026"
  echo "  particle : ele, pho, or pho_over_ele"
  echo "  centMin  : centrality minimum (hiBin), e.g. 0"
  echo "  centMax  : centrality maximum (hiBin), e.g. 200"
  echo "  sample   : QCDPhoton, DY, ZToEE (2025 only), JpsiToEE"
  echo ""
  echo "  Supported combinations:"
  echo "    ele         : DY, ZToEE (2025 only), JpsiToEE"
  echo "    pho         : QCDPhoton"
  echo "    pho_over_ele: DY, ZToEE (2025 only), JpsiToEE"
  echo ""
  echo "Examples:"
  echo "  ./triggerAnalysis.sh mc 2025 ele          0 200 ZToEE"
  echo "  ./triggerAnalysis.sh mc 2026 ele          0 200 DY"
  echo "  ./triggerAnalysis.sh mc 2026 ele          0 200 JpsiToEE"
  echo "  ./triggerAnalysis.sh mc 2026 pho          0 200 QCDPhoton"
  echo "  ./triggerAnalysis.sh mc 2026 pho_over_ele 0 200 DY"
  exit 1
}

if [[ $# -ne 6 ]]; then
  usage
fi

MODE="$1"
YEAR="$2"
PARTICLE="$3"
CENT_MIN="$4"
CENT_MAX="$5"
SAMPLE="$6"

if [[ "${MODE}" != "mc" ]]; then
  echo "Error: mode must be 'mc', got '${MODE}'"
  usage
fi

REPO_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "${REPO_DIR}"

PLOT_DIR="plots${YEAR}"

# ──────────────────────────────────────────────────
#  Build the right executable
# ──────────────────────────────────────────────────

EXEC="Execute_triggerAnalysis_${PARTICLE}_mc"

make "${EXEC}"
mkdir -p "${PLOT_DIR}" output

# ──────────────────────────────────────────────────
#  MC input forest paths
# ──────────────────────────────────────────────────
#
#  Paths keyed by (year, sample):
#
#  ZToEE    (2025 only — no 2026 forest available)
#    2025: JpsiDielectron ... crab_Run3_PbPb_2025MC_ZToEE/251115_102359/0000
#  DY       (2026 only)
#    2026: DrellYan_HighMass_MadGraph ... crab_Run3_PbPb_2026MC_DY/260515_160610/0000
#  JpsiToEE
#    2025: JpsiDielectron ... crab_Run3_PbPb_2025MC_JpsiToEE/251113_150758/0000
#    2026: JpsiDielectron ... crab_Run3_PbPb_2026MC_JpsiToEE/260514_234022/0000
#  QCDPhoton
#    2025: PhotonQCD      ... crab_Run3_PbPb_2025MC_QCDPhoton/251023_035328/0000
#    2026: PhotonQCD      ... crab_Run3_PbPb_2026MC_ZToEE/260514_202916/0000

MC_BASE_2025="eos/cms/store/group/phys_heavyions/kdeverea/Run3_PbPb_2025MC"
MC_BASE_2026="eos/cms/store/group/phys_heavyions/kdeverea/Run3_PbPb_2026MC"

FOREST_ZToEE_2025="/${MC_BASE_2025}/JpsiDielectron_pTHatMin4_HydjetEmbedded_Pythia8_TuneCP5_1510pre6/crab_Run3_PbPb_2025MC_ZToEE/251115_102359/0000"

FOREST_DY_2026="/${MC_BASE_2026}/DrellYan_HighMass_MadGraph_HydjetEmbedded_1610pre3/crab_Run3_PbPb_2026MC_DY/260515_160610/0000"

FOREST_JpsiToEE_2025="/${MC_BASE_2025}/JpsiDielectron_pTHatMin4_HydjetEmbedded_Pythia8_TuneCP5_1510pre6/crab_Run3_PbPb_2025MC_JpsiToEE/251113_150758/0000"
FOREST_JpsiToEE_2026="/${MC_BASE_2026}/JpsiDielectron_pTHatMin4_HydjetEmbedded_Pythia8_TuneCP5_1610pre3/crab_Run3_PbPb_2026MC_JpsiToEE/260514_234022/0000"

FOREST_QCDPhoton_2025="/${MC_BASE_2025}/PhotonQCD_pTMin20_HydjetEmbedded_Pythia8_TuneCP5_1510pre6/crab_Run3_PbPb_2025MC_QCDPhoton/251023_035328/0000"
FOREST_QCDPhoton_2026="/${MC_BASE_2026}/PhotonQCD_pTHatMin20_HydjetEmbedded_1610pre3/crab_Run3_PbPb_2026MC_ZToEE/260514_202916/0000"

# ──────────────────────────────────────────────────
#  Resolve INPUT_FOREST and NAMETAG for MC samples
# ──────────────────────────────────────────────────

resolve_mc_forest() {
  local sample="$1" year="$2"
  case "${sample}:${year}" in
    ZToEE:2025)     echo "${FOREST_ZToEE_2025}" ;;
    DY:2026)        echo "${FOREST_DY_2026}" ;;
    JpsiToEE:2025)  echo "${FOREST_JpsiToEE_2025}" ;;
    JpsiToEE:2026)  echo "${FOREST_JpsiToEE_2026}" ;;
    QCDPhoton:2025) echo "${FOREST_QCDPhoton_2025}" ;;
    QCDPhoton:2026) echo "${FOREST_QCDPhoton_2026}" ;;
    ZToEE:2026)
      echo "Error: ZToEE is not available for 2026. Use DY or JpsiToEE." >&2
      exit 1 ;;
    *)
      echo "Error: no forest path defined for sample=${sample} year=${year}" >&2
      exit 1 ;;
  esac
}

resolve_mc_nametag() {
  local sample="$1" year="$2"
  case "${sample}:${year}" in
    ZToEE:2025)     echo "${year} Pythia8+Hydjet Z->EE" ;;
    DY:2026)        echo "${year} MadGraph+Hydjet DY" ;;
    JpsiToEE:2025)  echo "${year} Pythia8+Hydjet J/psi->EE" ;;
    JpsiToEE:2026)  echo "${year} Pythia8+Hydjet J/psi->EE" ;;
    QCDPhoton:2025) echo "${year} Pythia8+Hydjet QCD Photon" ;;
    QCDPhoton:2026) echo "${year} Pythia8+Hydjet QCD Photon" ;;
    *)
      echo "Error: no nametag defined for sample=${sample} year=${year}" >&2
      exit 1 ;;
  esac
}

# ──────────────────────────────────────────────────
#  Dispatch
# ──────────────────────────────────────────────────

if [[ "${MODE}" == "mc" && "${PARTICLE}" == "ele" ]]; then

  valid_ele=false
  [[ "${SAMPLE}" == "JpsiToEE" ]] && valid_ele=true
  [[ "${SAMPLE}" == "DY"       ]] && valid_ele=true
  [[ "${SAMPLE}" == "ZToEE" && "${YEAR}" == "2025" ]] && valid_ele=true
  if [[ "${valid_ele}" == "false" ]]; then
    echo "Error: sample '${SAMPLE}' not supported for mc+ele (year=${YEAR})."
    echo "  2025: ZToEE, JpsiToEE"
    echo "  2026: DY, JpsiToEE"
    exit 1
  fi

  INPUT_FOREST="$(resolve_mc_forest "${SAMPLE}" "${YEAR}")"
  NAMETAG="$(resolve_mc_nametag "${SAMPLE}" "${YEAR}")"
  OUTPUT_BASE="MC${SAMPLE}${YEAR}_${CENT_MIN}_$((CENT_MAX / 2))"

  ./"${EXEC}" \
    --inputForest "${INPUT_FOREST}" \
    --inputEmulation "" \
    --outputBase "${OUTPUT_BASE}" \
    --nametag "${NAMETAG}" \
    --matchToForest true \
    --matchToEmulation false \
    --nfiles -1 \
    --minHiBin "${CENT_MIN}" \
    --maxHiBin "${CENT_MAX}" \
    --year "${YEAR}" \
    --plotDir "${PLOT_DIR}"


elif [[ "${MODE}" == "mc" && "${PARTICLE}" == "pho" ]]; then

  if [[ "${SAMPLE}" != "QCDPhoton" ]]; then
    echo "Error: sample '${SAMPLE}' not supported for mc+pho. Valid: QCDPhoton"
    exit 1
  fi

  INPUT_FOREST="$(resolve_mc_forest "${SAMPLE}" "${YEAR}")"
  NAMETAG="$(resolve_mc_nametag "${SAMPLE}" "${YEAR}")"
  OUTPUT_BASE="MC${SAMPLE}${YEAR}_${CENT_MIN}_$((CENT_MAX / 2))"

  ./"${EXEC}" \
    --inputForest "${INPUT_FOREST}" \
    --inputEmulation "" \
    --outputBase "${OUTPUT_BASE}" \
    --nametag "${NAMETAG}" \
    --matchToForest true \
    --matchToEmulation false \
    --nfiles -1 \
    --minHiBin "${CENT_MIN}" \
    --maxHiBin "${CENT_MAX}" \
    --year "${YEAR}" \
    --plotDir "${PLOT_DIR}"


elif [[ "${MODE}" == "mc" && "${PARTICLE}" == "pho_over_ele" ]]; then

  valid_poe=false
  [[ "${SAMPLE}" == "JpsiToEE" ]] && valid_poe=true
  [[ "${SAMPLE}" == "DY"       ]] && valid_poe=true
  [[ "${SAMPLE}" == "ZToEE" && "${YEAR}" == "2025" ]] && valid_poe=true
  if [[ "${valid_poe}" == "false" ]]; then
    echo "Error: sample '${SAMPLE}' not supported for mc+pho_over_ele (year=${YEAR})."
    echo "  2025: ZToEE, JpsiToEE"
    echo "  2026: DY, JpsiToEE"
    exit 1
  fi

  INPUT_FOREST="$(resolve_mc_forest "${SAMPLE}" "${YEAR}")"
  NAMETAG="$(resolve_mc_nametag "${SAMPLE}" "${YEAR}")"
  OUTPUT_BASE="MC${SAMPLE}${YEAR}_${CENT_MIN}_$((CENT_MAX / 2))"

  ./"${EXEC}" \
    --inputForest "${INPUT_FOREST}" \
    --outputBase "${OUTPUT_BASE}" \
    --nametag "${NAMETAG}" \
    --nfiles -1 \
    --minHiBin "${CENT_MIN}" \
    --maxHiBin "${CENT_MAX}" \
    --year "${YEAR}" \
    --plotDir "${PLOT_DIR}"


else
  echo "Error: unsupported particle '${PARTICLE}' for mc mode"
  echo "Supported: ele, pho, pho_over_ele"
  exit 1
fi

echo "Done. Plots written to ${PLOT_DIR}/"
