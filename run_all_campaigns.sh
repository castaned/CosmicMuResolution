#!/usr/bin/env bash

set -euo pipefail

SCRIPT_DIR="$(cd -- "$(dirname -- "${BASH_SOURCE[0]}")" && pwd)"
RUNNER="${SCRIPT_DIR}/run_all_scenarios.sh"

MUON_TYPES=(DGL DSA)
DATA_CAMPAIGNS=(2022F 2023D 2024I)
MC_CAMPAIGNS=(2022 2023 2024)

if [[ ! -x "${RUNNER}" ]]; then
    echo "Missing executable runner: ${RUNNER}" >&2
    exit 1
fi

cd "${SCRIPT_DIR}"

for muon_type in "${MUON_TYPES[@]}"; do
    for index in "${!DATA_CAMPAIGNS[@]}"; do
        data_campaign="${DATA_CAMPAIGNS[index]}"
        mc_campaign="${MC_CAMPAIGNS[index]}"

        echo
        echo "============================================================"
        echo "Running ${muon_type}: DATA ${data_campaign}, MC ${mc_campaign}"
        echo "============================================================"

        "${RUNNER}" \
            "${muon_type}" \
            "${data_campaign}" \
            "${mc_campaign}" \
            symmetric
    done
done

echo
echo "All DGL and DSA campaigns completed successfully."
