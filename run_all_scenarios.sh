#!/usr/bin/env bash

set -euo pipefail

PYTHON=python3

ANALYSIS_SCRIPT="dgm_sel_rdf.py"
COMPARE_MC_DATA_SCRIPT="compare_mc_data_binned.py"
COMPARE_SCENARIOS_SCRIPT="compare_scenarios_pt.py"

MUON_TYPE="${1:-DSA}"
DATA_CAMPAIGN="${2:-2024I}"
MC_CAMPAIGN="${3:-}"
SELECTION_MODE="${4:-${SELECTION_MODE:-legacy}}"
TREE_NAME="Events"
MIN_PT="12.5"

# =========================================
# BINNING (DGL vs DSA)
# =========================================

if [[ "${MUON_TYPE}" == "DSA" ]]; then
    PT_BINS="20,30,40,50,65,85,120,200,1000"
    DZ_BINS="1,5,10,20,30,45,60,150,300"
    DXY_BINS="1,5,10,20,30,50,80,160,300"
else
    PT_BINS="20,30,40,50,65,85,120,200,400,1000"
    DZ_BINS="1,5,10,20,30,45,60,100,150"
    DXY_BINS="1,5,10,20,30,40,60,80"
fi

case "${DATA_CAMPAIGN}" in
    2022F)
        DATA_DIR='/eos/user/c/castaned/Cosmics/Cosmics/CosmicsPPreco-CosmicDTLocalReco_Run2022F_Ntuples_v3'
        DEFAULT_MC_CAMPAIGN="2022"
        ;;
    2023D)
        DATA_DIR='/eos/user/c/castaned/Cosmics/Cosmics/CosmicsPPreco-CosmicDTLocalReco_Run2023D_Ntuples_v3'
        DEFAULT_MC_CAMPAIGN="2023"
        ;;
    2024I)
        DATA_DIR='/eos/user/c/castaned/Cosmics/Cosmics/CosmicsPPreco-CosmicDTLocalReco_Run2024I_Ntuples_v3'
        DEFAULT_MC_CAMPAIGN="2024"
        ;;
    *)
        echo "Unsupported DATA campaign: ${DATA_CAMPAIGN}"
        echo "Supported campaigns: 2022F, 2023D, 2024I"
        exit 1
        ;;
esac

if [[ -z "${MC_CAMPAIGN}" ]]; then
    MC_CAMPAIGN="${DEFAULT_MC_CAMPAIGN}"
fi

case "${MC_CAMPAIGN}" in
    2022)
        MC_DIR='/eos/user/c/castaned/Cosmics/LooseMuCosmic_Bin-P-10to3000-T0-Minus50to0_cosmuogen/CosmicsMC_Run3_2022_Ntuples_v2'
        ;;
    2023)
        MC_DIR='/eos/user/c/castaned/Cosmics/LooseMuCosmic_Bin-P-10to3000-T0-Minus50to0_cosmuogen/CosmicsMC_Run3_2023_Ntuples_v2'
        ;;
    2024)
        MC_DIR='/eos/user/c/castaned/Cosmics/LooseMuCosmic_Bin-P-10to3000-T0-Minus50to0_cosmuogen/CosmicsMC_Run3_2024_Ntuples_v2'
        ;;
    *)
        echo "Unsupported MC campaign: ${MC_CAMPAIGN}"
        echo "Supported MC campaigns: 2022, 2023, 2024"
        exit 1
        ;;
esac

case "${SELECTION_MODE}" in
    legacy|symmetric)
        ;;
    *)
        echo "Unsupported selection mode: ${SELECTION_MODE}"
        echo "Supported modes: legacy, symmetric"
        exit 1
        ;;
esac

OUTDIR_BASE="results_full_scan_${MUON_TYPE}_${DATA_CAMPAIGN}_${SELECTION_MODE}"

# INPUTS
MC_INPUT="${MC_INPUT:-${MC_DIR}/*/*/*.root}"
DATA_INPUT="${DATA_INPUT:-${DATA_DIR}/*/*/*.root}"

TRIGGER_ARGS=""
REL_PT_ERR_CUT="0.2"

if [[ "${MUON_TYPE}" == "DSA" ]]; then
    TIGHT_REL_PT_ERR="${TIGHT_REL_PT_ERR:-0.2}"
    TIGHT_MAX_CHI2="${TIGHT_MAX_CHI2:-5.0}"
    TIGHT_MIN_PRIMARY_HITS="${TIGHT_MIN_PRIMARY_HITS:-31}"
    TIGHT_MIN_SECONDARY_HITS="${TIGHT_MIN_SECONDARY_HITS:--1}"
else
    TIGHT_REL_PT_ERR="${TIGHT_REL_PT_ERR:-0.3}"
    TIGHT_MAX_CHI2="${TIGHT_MAX_CHI2:-5.0}"
    TIGHT_MIN_PRIMARY_HITS="${TIGHT_MIN_PRIMARY_HITS:-13}"
    TIGHT_MIN_SECONDARY_HITS="${TIGHT_MIN_SECONDARY_HITS:-6}"
fi

# =========================================
# UTILS
# =========================================

run_cmd() {
    echo
    echo ">>> $*"
    eval "$*"
}

build_common_args() {
    local datatype="$1"
    local outdir="$2"

    echo "--datatype ${datatype} \
          --muon-type ${MUON_TYPE} \
          --tree-name ${TREE_NAME} \
          --min-pt ${MIN_PT} \
          --pt-bins \"${PT_BINS}\" \
          --dz-bins \"${DZ_BINS}\" \
          --dxy-bins \"${DXY_BINS}\" \
          ${TRIGGER_ARGS} \
          --outdir ${outdir}"
}

# =========================================
# PREP
# =========================================

echo "Running full pipeline for ${MUON_TYPE} with DATA campaign ${DATA_CAMPAIGN} and MC campaign ${MC_CAMPAIGN}"
echo "Selection mode: ${SELECTION_MODE}"
echo "MC input: ${MC_INPUT}"
echo "DATA input: ${DATA_INPUT}"
if [[ "${SELECTION_MODE}" == "symmetric" ]]; then
    echo "Tight category cuts: ptErr/pt < ${TIGHT_REL_PT_ERR}, normalizedChi2 < ${TIGHT_MAX_CHI2}, primaryHits >= ${TIGHT_MIN_PRIMARY_HITS}, secondaryHits >= ${TIGHT_MIN_SECONDARY_HITS}"
else
    echo "Legacy tag/probe study: baseline, tagOnly, tagProbe with ptErr/pt < ${REL_PT_ERR_CUT}"
fi

if [[ ! -f "${ANALYSIS_SCRIPT}" ]]; then
    echo "Missing ${ANALYSIS_SCRIPT}"
    exit 1
fi

mkdir -p "${OUTDIR_BASE}"

# =========================================
# SCENARIO DIRECTORIES
# =========================================

if [[ "${SELECTION_MODE}" == "symmetric" ]]; then
    MC_RELAXED_DIR="${OUTDIR_BASE}/MC_relaxed"
    MC_TIGHT_DIR="${OUTDIR_BASE}/MC_tight"
    DATA_RELAXED_DIR="${OUTDIR_BASE}/DATA_relaxed"
    DATA_TIGHT_DIR="${OUTDIR_BASE}/DATA_tight"
    mkdir -p "${MC_RELAXED_DIR}" "${MC_TIGHT_DIR}"
    mkdir -p "${DATA_RELAXED_DIR}" "${DATA_TIGHT_DIR}"
else
    MC_BASELINE_DIR="${OUTDIR_BASE}/MC_baseline"
    MC_TAGONLY_DIR="${OUTDIR_BASE}/MC_tagOnly"
    MC_TAGPROBE_DIR="${OUTDIR_BASE}/MC_tagProbe"
    DATA_BASELINE_DIR="${OUTDIR_BASE}/DATA_baseline"
    DATA_TAGONLY_DIR="${OUTDIR_BASE}/DATA_tagOnly"
    DATA_TAGPROBE_DIR="${OUTDIR_BASE}/DATA_tagProbe"
    mkdir -p "${MC_BASELINE_DIR}" "${MC_TAGONLY_DIR}" "${MC_TAGPROBE_DIR}"
    mkdir -p "${DATA_BASELINE_DIR}" "${DATA_TAGONLY_DIR}" "${DATA_TAGPROBE_DIR}"
fi

# =========================================
# COMMON ARGS
# =========================================

if [[ "${SELECTION_MODE}" == "symmetric" ]]; then
    MC_RELAXED_ARGS=$(build_common_args "MC" "${MC_RELAXED_DIR}")
    MC_TIGHT_ARGS=$(build_common_args "MC" "${MC_TIGHT_DIR}")
    DATA_RELAXED_ARGS=$(build_common_args "DATA" "${DATA_RELAXED_DIR}")
    DATA_TIGHT_ARGS=$(build_common_args "DATA" "${DATA_TIGHT_DIR}")
else
    MC_BASELINE_ARGS=$(build_common_args "MC" "${MC_BASELINE_DIR}")
    MC_TAGONLY_ARGS=$(build_common_args "MC" "${MC_TAGONLY_DIR}")
    MC_TAGPROBE_ARGS=$(build_common_args "MC" "${MC_TAGPROBE_DIR}")
    DATA_BASELINE_ARGS=$(build_common_args "DATA" "${DATA_BASELINE_DIR}")
    DATA_TAGONLY_ARGS=$(build_common_args "DATA" "${DATA_TAGONLY_DIR}")
    DATA_TAGPROBE_ARGS=$(build_common_args "DATA" "${DATA_TAGPROBE_DIR}")
fi

# =========================================
# 1) MC
# =========================================

if [[ "${SELECTION_MODE}" == "symmetric" ]]; then
    run_cmd "${PYTHON} ${ANALYSIS_SCRIPT} ${MC_RELAXED_ARGS} \
        --input ${MC_INPUT} \
        --output Cosmics_muons_MC_${MUON_TYPE}_relaxed.root"

    run_cmd "${PYTHON} ${ANALYSIS_SCRIPT} ${MC_TIGHT_ARGS} \
        --input ${MC_INPUT} \
        --max-rel-pt-err-tag ${TIGHT_REL_PT_ERR} \
        --max-rel-pt-err-probe ${TIGHT_REL_PT_ERR} \
        --max-normalized-chi2-tag ${TIGHT_MAX_CHI2} \
        --max-normalized-chi2-probe ${TIGHT_MAX_CHI2} \
        --min-primary-hit-count-tag ${TIGHT_MIN_PRIMARY_HITS} \
        --min-primary-hit-count-probe ${TIGHT_MIN_PRIMARY_HITS} \
        --min-secondary-hit-count-tag ${TIGHT_MIN_SECONDARY_HITS} \
        --min-secondary-hit-count-probe ${TIGHT_MIN_SECONDARY_HITS} \
        --output Cosmics_muons_MC_${MUON_TYPE}_tight.root"
else
    run_cmd "${PYTHON} ${ANALYSIS_SCRIPT} ${MC_BASELINE_ARGS} \
        --input ${MC_INPUT} \
        --output Cosmics_muons_MC_${MUON_TYPE}_baseline.root"

    run_cmd "${PYTHON} ${ANALYSIS_SCRIPT} ${MC_TAGONLY_ARGS} \
        --input ${MC_INPUT} \
        --max-rel-pt-err-tag ${REL_PT_ERR_CUT} \
        --output Cosmics_muons_MC_${MUON_TYPE}_tagOnly.root"

    run_cmd "${PYTHON} ${ANALYSIS_SCRIPT} ${MC_TAGPROBE_ARGS} \
        --input ${MC_INPUT} \
        --max-rel-pt-err-tag ${REL_PT_ERR_CUT} \
        --max-rel-pt-err-probe ${REL_PT_ERR_CUT} \
        --output Cosmics_muons_MC_${MUON_TYPE}_tagProbe.root"
fi

# =========================================
# 2) DATA
# =========================================

if [[ "${SELECTION_MODE}" == "symmetric" ]]; then
    run_cmd "${PYTHON} ${ANALYSIS_SCRIPT} ${DATA_RELAXED_ARGS} \
        --input ${DATA_INPUT} \
        --output Cosmics_muons_DATA_${MUON_TYPE}_relaxed.root"

    run_cmd "${PYTHON} ${ANALYSIS_SCRIPT} ${DATA_TIGHT_ARGS} \
        --input ${DATA_INPUT} \
        --max-rel-pt-err-tag ${TIGHT_REL_PT_ERR} \
        --max-rel-pt-err-probe ${TIGHT_REL_PT_ERR} \
        --max-normalized-chi2-tag ${TIGHT_MAX_CHI2} \
        --max-normalized-chi2-probe ${TIGHT_MAX_CHI2} \
        --min-primary-hit-count-tag ${TIGHT_MIN_PRIMARY_HITS} \
        --min-primary-hit-count-probe ${TIGHT_MIN_PRIMARY_HITS} \
        --min-secondary-hit-count-tag ${TIGHT_MIN_SECONDARY_HITS} \
        --min-secondary-hit-count-probe ${TIGHT_MIN_SECONDARY_HITS} \
        --output Cosmics_muons_DATA_${MUON_TYPE}_tight.root"
else
    run_cmd "${PYTHON} ${ANALYSIS_SCRIPT} ${DATA_BASELINE_ARGS} \
        --input ${DATA_INPUT} \
        --output Cosmics_muons_DATA_${MUON_TYPE}_baseline.root"

    run_cmd "${PYTHON} ${ANALYSIS_SCRIPT} ${DATA_TAGONLY_ARGS} \
        --input ${DATA_INPUT} \
        --max-rel-pt-err-tag ${REL_PT_ERR_CUT} \
        --output Cosmics_muons_DATA_${MUON_TYPE}_tagOnly.root"

    run_cmd "${PYTHON} ${ANALYSIS_SCRIPT} ${DATA_TAGPROBE_ARGS} \
        --input ${DATA_INPUT} \
        --max-rel-pt-err-tag ${REL_PT_ERR_CUT} \
        --max-rel-pt-err-probe ${REL_PT_ERR_CUT} \
        --output Cosmics_muons_DATA_${MUON_TYPE}_tagProbe.root"
fi

# =========================================
# 3) SCENARIO COMPARISON
# =========================================

if [[ "${SELECTION_MODE}" == "symmetric" ]]; then
    run_cmd "${PYTHON} ${COMPARE_SCENARIOS_SCRIPT} \
        --inputs \
        ${MC_RELAXED_DIR}/Cosmics_muons_MC_${MUON_TYPE}_relaxed.root \
        ${MC_TIGHT_DIR}/Cosmics_muons_MC_${MUON_TYPE}_tight.root \
        --labels \
        relaxed \
        tight \
        --names \
        relaxed \
        tight \
        --outdir ${OUTDIR_BASE}/compare_scenarios_MC_${MUON_TYPE}"

    run_cmd "${PYTHON} ${COMPARE_SCENARIOS_SCRIPT} \
        --inputs \
        ${DATA_RELAXED_DIR}/Cosmics_muons_DATA_${MUON_TYPE}_relaxed.root \
        ${DATA_TIGHT_DIR}/Cosmics_muons_DATA_${MUON_TYPE}_tight.root \
        --labels \
        relaxed \
        tight \
        --names \
        relaxed \
        tight \
        --outdir ${OUTDIR_BASE}/compare_scenarios_DATA_${MUON_TYPE}"
else
    run_cmd "${PYTHON} ${COMPARE_SCENARIOS_SCRIPT} \
        --baseline ${MC_BASELINE_DIR}/Cosmics_muons_MC_${MUON_TYPE}_baseline.root \
        --tagonly ${MC_TAGONLY_DIR}/Cosmics_muons_MC_${MUON_TYPE}_tagOnly.root \
        --tagprobe ${MC_TAGPROBE_DIR}/Cosmics_muons_MC_${MUON_TYPE}_tagProbe.root \
        --outdir ${OUTDIR_BASE}/compare_scenarios_MC_${MUON_TYPE}"

    run_cmd "${PYTHON} ${COMPARE_SCENARIOS_SCRIPT} \
        --baseline ${DATA_BASELINE_DIR}/Cosmics_muons_DATA_${MUON_TYPE}_baseline.root \
        --tagonly ${DATA_TAGONLY_DIR}/Cosmics_muons_DATA_${MUON_TYPE}_tagOnly.root \
        --tagprobe ${DATA_TAGPROBE_DIR}/Cosmics_muons_DATA_${MUON_TYPE}_tagProbe.root \
        --outdir ${OUTDIR_BASE}/compare_scenarios_DATA_${MUON_TYPE}"
fi

# =========================================
# 4) MC vs DATA
# =========================================

if [[ "${SELECTION_MODE}" == "symmetric" ]]; then
    SCENARIOS=(relaxed tight)
else
    SCENARIOS=(baseline tagOnly tagProbe)
fi

for scenario in "${SCENARIOS[@]}"; do

run_cmd "${PYTHON} ${COMPARE_MC_DATA_SCRIPT} \
    --mc ${OUTDIR_BASE}/MC_${scenario}/Cosmics_muons_MC_${MUON_TYPE}_${scenario}.root \
    --data ${OUTDIR_BASE}/DATA_${scenario}/Cosmics_muons_DATA_${MUON_TYPE}_${scenario}.root \
    --period ${DATA_CAMPAIGN} \
    --outdir ${OUTDIR_BASE}/compare_MC_DATA_${MUON_TYPE}_${scenario}"

done

# =========================================
# DONE
# =========================================

echo
echo "======================================="
echo "DONE for ${MUON_TYPE}"
echo "Output: ${OUTDIR_BASE}"
echo "======================================="
