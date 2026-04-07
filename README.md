# CosmicMuResolution

Cosmic muon momentum-resolution study comparing `DATA` and `MC` samples for `DGL` and `DSA` muons.

## Current Workflow

The main implementation is [dgm_sel_rdf.py](./dgm_sel_rdf.py), a PyROOT `RDataFrame` pipeline that:

- reads ROOT inputs from files, directories, or glob patterns
- applies the tag-and-probe style event selection
- computes q/pT residual distributions
- fits those distributions in `pt`, `dz`, and `dxy` bins
- writes a ROOT output file plus PNG control/comparison plots

### Run The Analysis

Example:

```bash
python3 dgm_sel_rdf.py \
  --input "/path/to/input/*.root" \
  --datatype DATA \
  --muon-type DSA \
  --tree-name Events \
  --outdir results
```

Important options:

- `--datatype {MC,DATA}` selects sample type.
- `--muon-type {DSA,DGL}` selects the muon reconstruction.
- `--pt-bins`, `--dz-bins`, and `--dxy-bins` override the default binning.
- `--max-rel-pt-err-tag` and `--max-rel-pt-err-probe` enable tighter scenario studies.
- `--output` sets the ROOT file name written inside `--outdir`.

By default the script writes:

- `results/<output-root-file>`
- `results/<DATA-or-MC>/<DSA-or-DGL>/Control_plots/*.png`
- `results/<DATA-or-MC>/<DSA-or-DGL>/pt_Plots/*.png`
- `results/<DATA-or-MC>/<DSA-or-DGL>/dz_Plots/*.png`
- `results/<DATA-or-MC>/<DSA-or-DGL>/dxy_Plots/*.png`

## Batch Workflow

Use [run_all_scenarios.sh](./run_all_scenarios.sh) to run the full modern workflow for one muon type:

```bash
bash run_all_scenarios.sh DSA
```

This script:

- runs baseline, tag-only, and tag+probe scenarios for `MC`
- runs baseline, tag-only, and tag+probe scenarios for `DATA`
- compares scenarios with [compare_scenarios_pt.py](./compare_scenarios_pt.py)
- compares `MC` vs `DATA` with [compare_mc_data_binned.py](./compare_mc_data_binned.py)

The EOS input paths inside `run_all_scenarios.sh` are hard-coded and should be adjusted for your environment.

## Repository Layout

- [dgm_sel_rdf.py](./dgm_sel_rdf.py): current analysis implementation
- [run_all_scenarios.sh](./run_all_scenarios.sh): batch runner for the modern workflow
- [compare_mc_data_binned.py](./compare_mc_data_binned.py): overlays `MC` and `DATA` mean/sigma results for `pt`, `dz`, and `dxy`
- [compare_scenarios_pt.py](./compare_scenarios_pt.py): compares baseline, tag-only, and tag+probe scenarios
- `DATA/` and `MC/`: generated plot artifacts
- `old/legacy_workflow/`: archived macro-era scripts and top-level legacy plots

## Requirements

This repository assumes a local CERN ROOT / PyROOT environment. There is currently no packaged environment file or dependency manifest in the repo.
