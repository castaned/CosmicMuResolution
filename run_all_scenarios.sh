#!/usr/bin/env bash

set -euo pipefail

PYTHON=python3

ANALYSIS_SCRIPT="dgm_sel_rdf.py"
COMPARE_MC_DATA_SCRIPT="compare_mc_data_binned.py"
COMPARE_SCENARIOS_SCRIPT="compare_scenarios_pt.py"

MUON_TYPE="${1:-DSA}"
DATA_CAMPAIGN="${2:-2024I}"
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
        DATA_DIR='/eos/user/c/castaned/Cosmics/Cosmics/CosmicsPPreco-CosmicDTLocalReco_Run2022F_Ntuples_v2'
        ;;
    2023D)
        DATA_DIR='/eos/user/c/castaned/Cosmics/Cosmics/CosmicsPPreco-CosmicDTLocalReco_Run2023D_Ntuples_v2'
        ;;
    2024I)
        DATA_DIR='/eos/user/c/castaned/Cosmics/Cosmics/CosmicsPPreco-CosmicDTLocalReco_Run2024I_Ntuples_v4'
        ;;
    *)
        echo "Unsupported DATA campaign: ${DATA_CAMPAIGN}"
        echo "Supported campaigns: 2022F, 2023D, 2024I"
        exit 1
        ;;
esac

OUTDIR_BASE="results_full_scan_${MUON_TYPE}_${DATA_CAMPAIGN}"

# INPUTS
MC_INPUT="${MC_INPUT:-/eos/user/c/castaned/Cosmics/LooseMuCosmic_Bin-P-10to3000-T0-Minus50to0_cosmuogen/CosmicsMC_Run3_2024_Ntuples_v2/260310_205717/0000/*.root}"
DATA_INPUT="${DATA_INPUT:-${DATA_DIR}/*/*/*.root}"

REL_PT_ERR_CUT="0.2"
TRIGGER_ARGS=""

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

echo "Running full pipeline for ${MUON_TYPE} with DATA campaign ${DATA_CAMPAIGN}"
echo "MC input: ${MC_INPUT}"
echo "DATA input: ${DATA_INPUT}"

if [[ ! -f "${ANALYSIS_SCRIPT}" ]]; then
    echo "Missing ${ANALYSIS_SCRIPT}"
    exit 1
fi

mkdir -p "${OUTDIR_BASE}"

# =========================================
# SCENARIO DIRECTORIES
# =========================================

MC_BASELINE_DIR="${OUTDIR_BASE}/MC_baseline"
MC_TAGONLY_DIR="${OUTDIR_BASE}/MC_tagOnly"
MC_TAGPROBE_DIR="${OUTDIR_BASE}/MC_tagProbe"

DATA_BASELINE_DIR="${OUTDIR_BASE}/DATA_baseline"
DATA_TAGONLY_DIR="${OUTDIR_BASE}/DATA_tagOnly"
DATA_TAGPROBE_DIR="${OUTDIR_BASE}/DATA_tagProbe"

mkdir -p "${MC_BASELINE_DIR}" "${MC_TAGONLY_DIR}" "${MC_TAGPROBE_DIR}"
mkdir -p "${DATA_BASELINE_DIR}" "${DATA_TAGONLY_DIR}" "${DATA_TAGPROBE_DIR}"

# =========================================
# COMMON ARGS
# =========================================

MC_BASELINE_ARGS=$(build_common_args "MC" "${MC_BASELINE_DIR}")
MC_TAGONLY_ARGS=$(build_common_args "MC" "${MC_TAGONLY_DIR}")
MC_TAGPROBE_ARGS=$(build_common_args "MC" "${MC_TAGPROBE_DIR}")

DATA_BASELINE_ARGS=$(build_common_args "DATA" "${DATA_BASELINE_DIR}")
DATA_TAGONLY_ARGS=$(build_common_args "DATA" "${DATA_TAGONLY_DIR}")
DATA_TAGPROBE_ARGS=$(build_common_args "DATA" "${DATA_TAGPROBE_DIR}")

# =========================================
# 1) MC
# =========================================

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

# =========================================
# 2) DATA
# =========================================

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

# =========================================
# 3) SCENARIO COMPARISON
# =========================================

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

# =========================================
# 4) MC vs DATA
# =========================================

for scenario in baseline tagOnly tagProbe; do

run_cmd "${PYTHON} ${COMPARE_MC_DATA_SCRIPT} \
    --mc ${OUTDIR_BASE}/MC_${scenario}/Cosmics_muons_MC_${MUON_TYPE}_${scenario}.root \
    --data ${OUTDIR_BASE}/DATA_${scenario}/Cosmics_muons_DATA_${MUON_TYPE}_${scenario}.root \
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
