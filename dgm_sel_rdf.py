#!/usr/bin/env python3

import ROOT
import os
import glob
import argparse
import array
import csv
import math

ROOT.gStyle.SetOptStat(0)

# ============================================================
# Defaults
# ============================================================
DEFAULT_PT_BINS = [20.0, 30.0, 40.0, 50.0, 65.0, 85.0, 120.0, 200.0, 400.0, 1000.0]
DEFAULT_DZ_BINS = [1.0,5.0, 10.0, 20.0, 30.0, 45.0, 60.0,100,150.0]
DEFAULT_DXY_BINS = [1.0,5.0, 10.0, 20.0, 30.0, 40.0,60.0, 80.0]
HYBRID_DOUBLE_MIN_ENTRIES = 25
HYBRID_MAX_REL_SIGMA_ERR = 0.5
HYBRID_MAX_SIGMA_ERR_RATIO_TO_SINGLE = 2.0
HYBRID_MIN_CORE_FRACTION = 0.10
HYBRID_MAX_CORE_FRACTION = 0.90
HYBRID_MAX_TAIL_TO_CORE_RATIO = 6.0


# ============================================================
# Helpers
# ============================================================
def parse_bin_list(text):
    vals = [float(x.strip()) for x in text.split(",") if x.strip()]
    if len(vals) < 2:
        raise argparse.ArgumentTypeError("Bin list must contain at least two edges.")
    return vals


def expand_inputs(inputs):
    files = []
    for item in inputs:
        if os.path.isdir(item):
            files.extend(sorted(glob.glob(os.path.join(item, "*.root"))))
        else:
            matched = glob.glob(item)
            if matched:
                files.extend(sorted(matched))
            elif item.endswith(".root") and os.path.isfile(item):
                files.append(item)
    return sorted(set(files))


def ensure_dir(path):
    os.makedirs(path, exist_ok=True)


def array_to_c(pylist):
    return array.array("d", pylist)


# ============================================================
# C++ helpers injected into ROOT
# ============================================================
ROOT.gInterpreter.Declare(r'''
#include <ROOT/RVec.hxx>
#include <cmath>
#include <vector>
#include <string>
#include <algorithm>

using namespace ROOT::VecOps;

struct AnalysisOutput {
    RVec<double> resolutions;
    RVec<double> symmetric_resolutions;
    RVec<double> tag_pt;
    RVec<double> tag_eta;
    RVec<double> tag_phi;
    RVec<double> tag_charge;
    RVec<double> tag_pt_err_over_pt;
    RVec<double> probe_pt;
    RVec<double> probe_eta;
    RVec<double> probe_phi;
    RVec<double> probe_charge;
    RVec<double> probe_pt_err_over_pt;
    RVec<double> all_muon_pt;
    RVec<double> tag_dz_abs;
    RVec<double> tag_dxy_abs;
};

template <typename T>
RVec<double> filter_by_bin(const RVec<double>& values, const RVec<T>& refs, double low, double high) {
    RVec<double> out;
    const auto n = std::min(values.size(), refs.size());
    out.reserve(n);
    for (size_t i = 0; i < n; ++i) {
        if (refs[i] >= low && refs[i] < high) out.push_back(values[i]);
    }
    return out;
}

AnalysisOutput analyze_event(
    int ndmu,
    bool hlt,
    const RVec<float>& pt,
    const RVec<float>& ptErr,
    const RVec<float>& normalizedChi2,
    const RVec<int>& primaryHitCount,
    const RVec<int>& secondaryHitCount,
    const RVec<float>& eta,
    const RVec<float>& phi,
    const RVec<float>& dz,
    const RVec<float>& dxy,
    const RVec<float>& charge,
    const RVec<char>& passTagID,
    const RVec<int>& probeID,
    const RVec<char>& hasProbe,
    const RVec<int>& isRecoType,
    bool requireTrigger,
    double min_pt,
    double max_rel_pt_err_tag,
    double max_rel_pt_err_probe,
    double max_normalized_chi2_tag,
    double max_normalized_chi2_probe,
    int min_primary_hit_count_tag,
    int min_primary_hit_count_probe,
    int min_secondary_hit_count_tag,
    int min_secondary_hit_count_probe
) {
    AnalysisOutput out;

    if (requireTrigger && !hlt) return out;
    if (ndmu < 2) return out;

    for (int i = 0; i < ndmu; ++i) {
        out.all_muon_pt.push_back(pt[i]);

        if (pt[i] < min_pt) continue;
        if (!passTagID[i]) continue;
        if (isRecoType[i] != 1) continue;
        if (!hasProbe[i]) continue;

        const int j = probeID[i];
        if (j < 0 || j >= ndmu) continue;
        if (pt[j] < min_pt) continue;
        if (pt[i] == 0.f || pt[j] == 0.f) continue;

        const double tag_rel_pt_err = ptErr[i] / pt[i];
        const double probe_rel_pt_err = ptErr[j] / pt[j];
        const double tag_normalized_chi2 = normalizedChi2[i];
        const double probe_normalized_chi2 = normalizedChi2[j];
        const int tag_primary_hit_count = primaryHitCount[i];
        const int probe_primary_hit_count = primaryHitCount[j];
        const int tag_secondary_hit_count = secondaryHitCount[i];
        const int probe_secondary_hit_count = secondaryHitCount[j];

        if (max_rel_pt_err_tag >= 0.0 && tag_rel_pt_err >= max_rel_pt_err_tag) continue;
        if (max_rel_pt_err_probe >= 0.0 && probe_rel_pt_err >= max_rel_pt_err_probe) continue;
        if (max_normalized_chi2_tag >= 0.0 && tag_normalized_chi2 >= max_normalized_chi2_tag) continue;
        if (max_normalized_chi2_probe >= 0.0 && probe_normalized_chi2 >= max_normalized_chi2_probe) continue;
        if (min_primary_hit_count_tag >= 0 && tag_primary_hit_count < min_primary_hit_count_tag) continue;
        if (min_primary_hit_count_probe >= 0 && probe_primary_hit_count < min_primary_hit_count_probe) continue;
        if (min_secondary_hit_count_tag >= 0 && tag_secondary_hit_count < min_secondary_hit_count_tag) continue;
        if (min_secondary_hit_count_probe >= 0 && probe_secondary_hit_count < min_secondary_hit_count_probe) continue;

        const double inv_up = std::abs(charge[j] / pt[j]);
        const double inv_down = std::abs(charge[i] / pt[i]);
        if (inv_down == 0.) continue;

        const double resolution = (inv_up - inv_down) / (inv_down);
        const double inv_avg = 0.5 * (inv_up + inv_down);
        if (inv_avg == 0.) continue;
        const double symmetric_resolution = (inv_up - inv_down) / inv_avg;

        out.resolutions.push_back(resolution);
        out.symmetric_resolutions.push_back(symmetric_resolution);

        out.tag_pt.push_back(pt[i]);
        out.tag_eta.push_back(eta[i]);
        out.tag_phi.push_back(phi[i]);
        out.tag_charge.push_back(charge[i]);
        out.tag_dz_abs.push_back(std::abs(dz[i]));
        out.tag_dxy_abs.push_back(std::abs(dxy[i]));
        out.tag_pt_err_over_pt.push_back(tag_rel_pt_err);

        out.probe_pt.push_back(pt[j]);
        out.probe_eta.push_back(eta[j]);
        out.probe_phi.push_back(phi[j]);
        out.probe_charge.push_back(charge[j]);
        out.probe_pt_err_over_pt.push_back(probe_rel_pt_err);
    }

    return out;
}
''')


def write_fit_csv(csv_path, rows):
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "bin_low", "bin_high", "entries",
            "mean", "mean_err",
            "sigma", "sigma_err",
            "chi2", "ndf", "chi2_ndf"
        ])
        writer.writerows(rows)


def write_double_fit_csv(csv_path, rows):
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "bin_low", "bin_high", "entries",
            "mean", "mean_err",
            "sigma_eff", "sigma_eff_err",
            "sigma_core", "sigma_core_err",
            "sigma_tail", "sigma_tail_err",
            "frac_core", "frac_core_err",
            "chi2", "ndf", "chi2_ndf",
            "fit_status"
        ])
        writer.writerows(rows)


def write_hybrid_fit_csv(csv_path, rows):
    with open(csv_path, "w", newline="") as f:
        writer = csv.writer(f)
        writer.writerow([
            "bin_low", "bin_high", "entries",
            "mean", "mean_err",
            "sigma", "sigma_err",
            "chi2", "ndf", "chi2_ndf",
            "fit_model", "selection_reason",
            "single_fit_status", "double_fit_status"
        ])
        writer.writerows(rows)


def fit_single_gaussian(hist, fit_name, fit_min, fit_max):
    fit = ROOT.TF1(fit_name, "gaus", fit_min, fit_max)

    entries = hist.GetEntries()
    result = None
    mean = sigma = mean_err = sigma_err = 0.0
    chi2 = chi2_ndf = 0.0
    ndf = 0
    fit_status = -1

    if entries > 5:
        result = hist.Fit(fit, "QSR")
        if result:
            fit_status = int(result.Status())
        mean = fit.GetParameter(1)
        sigma = abs(fit.GetParameter(2))
        mean_err = fit.GetParError(1)
        sigma_err = fit.GetParError(2)
        chi2 = fit.GetChisquare()
        ndf = fit.GetNDF()
        chi2_ndf = chi2 / ndf if ndf > 0 else 0.0

    return {
        "fit": fit,
        "result": result,
        "entries": entries,
        "fit_status": fit_status,
        "mean": mean,
        "mean_err": mean_err,
        "sigma": sigma,
        "sigma_err": sigma_err,
        "chi2": chi2,
        "ndf": ndf,
        "chi2_ndf": chi2_ndf,
    }


def fit_double_gaussian(hist, fit_name, fit_min, fit_max):
    entries = hist.GetEntries()
    fit = ROOT.TF1(fit_name, "gaus(0) + gaus(3)", fit_min, fit_max)
    fit.SetLineColor(ROOT.kBlue + 1)

    result = None
    fit_status = -1
    mean = mean_err = 0.0
    sigma_eff = sigma_eff_err = 0.0
    sigma_core = sigma_core_err = 0.0
    sigma_tail = sigma_tail_err = 0.0
    frac_core = frac_core_err = 0.0
    chi2 = chi2_ndf = 0.0
    ndf = 0

    if entries <= 10:
        return {
            "fit": fit,
            "result": result,
            "entries": entries,
            "fit_status": fit_status,
            "mean": mean,
            "mean_err": mean_err,
            "sigma_eff": sigma_eff,
            "sigma_eff_err": sigma_eff_err,
            "sigma_core": sigma_core,
            "sigma_core_err": sigma_core_err,
            "sigma_tail": sigma_tail,
            "sigma_tail_err": sigma_tail_err,
            "frac_core": frac_core,
            "frac_core_err": frac_core_err,
            "chi2": chi2,
            "ndf": ndf,
            "chi2_ndf": chi2_ndf,
        }

    hist_mean = hist.GetMean()
    hist_rms = max(hist.GetRMS(), 1e-4)
    hist_max = max(hist.GetMaximum(), 1.0)
    core_sigma_guess = max(0.5 * hist_rms, 1e-4)
    tail_sigma_guess = max(1.5 * hist_rms, core_sigma_guess * 1.2)

    fit.SetParNames("A1", "mean1", "sigma1", "A2", "mean2", "sigma2")
    fit.SetParameters(0.7 * hist_max, hist_mean, core_sigma_guess, 0.3 * hist_max, hist_mean, tail_sigma_guess)
    fit.SetParLimits(0, 0.0, max(10.0 * hist_max, 1.0))
    fit.SetParLimits(1, fit_min, fit_max)
    fit.SetParLimits(2, 1e-4, max(hist_rms * 5.0, (fit_max - fit_min)))
    fit.SetParLimits(3, 0.0, max(10.0 * hist_max, 1.0))
    fit.SetParLimits(4, fit_min, fit_max)
    fit.SetParLimits(5, 1e-4, max(hist_rms * 8.0, (fit_max - fit_min)))

    result = hist.Fit(fit, "QSR")
    if result:
        fit_status = int(result.Status())

    amp1 = fit.GetParameter(0)
    mean1 = fit.GetParameter(1)
    sigma1 = abs(fit.GetParameter(2))
    amp2 = fit.GetParameter(3)
    mean2 = fit.GetParameter(4)
    sigma2 = abs(fit.GetParameter(5))

    err_mean1 = fit.GetParError(1)
    err_sigma1 = fit.GetParError(2)
    err_mean2 = fit.GetParError(4)
    err_sigma2 = fit.GetParError(5)
    err_amp1 = fit.GetParError(0)
    err_amp2 = fit.GetParError(3)

    # Report the narrower component as the core to keep outputs stable across bins.
    if sigma1 <= sigma2:
        core_amp, core_mean, core_sigma = amp1, mean1, sigma1
        tail_amp, tail_mean, tail_sigma = amp2, mean2, sigma2
        core_mean_err, core_sigma_err = err_mean1, err_sigma1
        tail_mean_err, tail_sigma_err = err_mean2, err_sigma2
        core_amp_err, tail_amp_err = err_amp1, err_amp2
    else:
        core_amp, core_mean, core_sigma = amp2, mean2, sigma2
        tail_amp, tail_mean, tail_sigma = amp1, mean1, sigma1
        core_mean_err, core_sigma_err = err_mean2, err_sigma2
        tail_mean_err, tail_sigma_err = err_mean1, err_sigma1
        core_amp_err, tail_amp_err = err_amp2, err_amp1

    total_amp = core_amp + tail_amp
    if total_amp > 0.0:
        frac_core = core_amp / total_amp
        frac_tail = tail_amp / total_amp
    else:
        frac_core = 0.0
        frac_tail = 0.0

    mean = frac_core * core_mean + frac_tail * tail_mean
    second_moment = (
        frac_core * (core_sigma ** 2 + core_mean ** 2)
        + frac_tail * (tail_sigma ** 2 + tail_mean ** 2)
    )
    variance = max(second_moment - mean ** 2, 0.0)
    sigma_eff = variance ** 0.5

    if total_amp > 0.0:
        frac_core_err = ((tail_amp * core_amp_err) ** 2 + (core_amp * tail_amp_err) ** 2) ** 0.5 / (total_amp ** 2)
    else:
        frac_core_err = 0.0

    sigma_eff_err = (
        (frac_core * core_sigma_err) ** 2
        + (frac_tail * tail_sigma_err) ** 2
    ) ** 0.5
    mean_err = ((frac_core * core_mean_err) ** 2 + (frac_tail * tail_mean_err) ** 2) ** 0.5

    chi2 = fit.GetChisquare()
    ndf = fit.GetNDF()
    chi2_ndf = chi2 / ndf if ndf > 0 else 0.0

    return {
        "fit": fit,
        "result": result,
        "entries": entries,
        "fit_status": fit_status,
        "mean": mean,
        "mean_err": mean_err,
        "sigma_eff": sigma_eff,
        "sigma_eff_err": sigma_eff_err,
        "sigma_core": core_sigma,
        "sigma_core_err": core_sigma_err,
        "sigma_tail": tail_sigma,
        "sigma_tail_err": tail_sigma_err,
        "frac_core": frac_core,
        "frac_core_err": frac_core_err,
        "chi2": chi2,
        "ndf": ndf,
        "chi2_ndf": chi2_ndf,
    }


def relative_error(value, error):
    if not math.isfinite(value) or value == 0.0:
        return float("inf")
    if not math.isfinite(error) or error < 0.0:
        return float("inf")
    return abs(error / value)


def choose_hybrid_result(single_result, double_result):
    reason = "double"
    use_double = True

    single_rel_sigma_err = relative_error(single_result["sigma"], single_result["sigma_err"])
    double_rel_sigma_err = relative_error(double_result["sigma_eff"], double_result["sigma_eff_err"])
    sigma_err_ratio_to_single = float("inf")
    if math.isfinite(single_result["sigma_err"]) and single_result["sigma_err"] > 0.0:
        sigma_err_ratio_to_single = double_result["sigma_eff_err"] / single_result["sigma_err"]

    if double_result["entries"] < HYBRID_DOUBLE_MIN_ENTRIES:
        use_double = False
        reason = "low_entries"
    elif double_result["fit_status"] != 0:
        use_double = False
        reason = "fit_status"
    elif double_result["ndf"] <= 0:
        use_double = False
        reason = "ndf"
    elif not math.isfinite(double_result["sigma_eff"]) or double_result["sigma_eff"] <= 0.0:
        use_double = False
        reason = "sigma_eff"
    elif double_rel_sigma_err > HYBRID_MAX_REL_SIGMA_ERR:
        use_double = False
        reason = "sigma_eff_err"
    elif not math.isfinite(double_result["frac_core"]) or not (0.0 < double_result["frac_core"] < 1.0):
        use_double = False
        reason = "frac_core"
    elif not (HYBRID_MIN_CORE_FRACTION <= double_result["frac_core"] <= HYBRID_MAX_CORE_FRACTION):
        use_double = False
        reason = "core_fraction_extreme"
    elif (
        not math.isfinite(double_result["sigma_core"]) or double_result["sigma_core"] <= 0.0
        or not math.isfinite(double_result["sigma_tail"]) or double_result["sigma_tail"] <= 0.0
    ):
        use_double = False
        reason = "component_sigma"
    elif double_result["sigma_tail"] < double_result["sigma_core"]:
        use_double = False
        reason = "component_order"
    elif double_result["sigma_tail"] / double_result["sigma_core"] > HYBRID_MAX_TAIL_TO_CORE_RATIO:
        use_double = False
        reason = "tail_too_broad"
    elif (
        math.isfinite(sigma_err_ratio_to_single)
        and sigma_err_ratio_to_single > HYBRID_MAX_SIGMA_ERR_RATIO_TO_SINGLE
    ):
        use_double = False
        reason = "abs_err_vs_single"

    if use_double:
        return {
            "fit": double_result["fit"],
            "entries": double_result["entries"],
            "mean": double_result["mean"],
            "mean_err": double_result["mean_err"],
            "sigma": double_result["sigma_eff"],
            "sigma_err": double_result["sigma_eff_err"],
            "chi2": double_result["chi2"],
            "ndf": double_result["ndf"],
            "chi2_ndf": double_result["chi2_ndf"],
            "fit_model": "double",
            "selection_reason": reason,
        }

    return {
        "fit": single_result["fit"],
        "entries": single_result["entries"],
        "mean": single_result["mean"],
        "mean_err": single_result["mean_err"],
        "sigma": single_result["sigma"],
        "sigma_err": single_result["sigma_err"],
        "chi2": single_result["chi2"],
        "ndf": single_result["ndf"],
        "chi2_ndf": single_result["chi2_ndf"],
        "fit_model": "single",
        "selection_reason": reason,
    }


def draw_graph_and_save(graph, canvas_name, out_file, out_path, logx):
    canvas = ROOT.TCanvas(canvas_name, "", 1000, 750)
    canvas.SetGrid()
    if logx:
        canvas.SetLogx()
    canvas.SetLeftMargin(0.18)
    graph.Draw("AP")
    out_file.cd()
    graph.Write()
    canvas.SaveAs(out_path)


# ============================================================
# Plot and fit
# ============================================================
def fit_and_draw(
    hlist,
    bin_type,
    x_title,
    base_dir,
    muon_type,
    res_min,
    res_max,
    out_file,
    logx=True,
    name_suffix="",
    residual_axis_label="q/p_{T} residual",
    residual_summary_label="q/p_{T} relative residual",
):
    means, mean_errs = [], []
    sigmas, sigma_errs = [], []
    chi2ndf_vals = []
    hybrid_means, hybrid_mean_errs = [], []
    hybrid_sigmas, hybrid_sigma_errs = [], []
    hybrid_chi2ndf_vals = []
    dg_means, dg_mean_errs = [], []
    dg_sigma_eff, dg_sigma_eff_errs = [], []
    dg_sigma_core, dg_sigma_core_errs = [], []
    dg_sigma_tail, dg_sigma_tail_errs = [], []
    dg_frac_core, dg_frac_core_errs = [], []
    dg_chi2ndf_vals = []
    centers, halfwidths = [], []
    csv_rows = []
    hybrid_csv_rows = []
    double_csv_rows = []

    canvas_all = ROOT.TCanvas(f"All_plots_{bin_type}{name_suffix}", f"{bin_type} Histograms{name_suffix}", 1500, 1000)
    canvas_all_hybrid = ROOT.TCanvas(f"All_plots_{bin_type}{name_suffix}_hybrid", f"{bin_type} Histograms Hybrid{name_suffix}", 1500, 1000)
    canvas_all_double = ROOT.TCanvas(f"All_plots_{bin_type}{name_suffix}_double", f"{bin_type} Histograms Double{name_suffix}", 1500, 1000)
    n = len(hlist)
    nx = 3
    ny = (n + nx - 1) // nx
    canvas_all.Divide(nx, ny)
    canvas_all_hybrid.Divide(nx, ny)
    canvas_all_double.Divide(nx, ny)

    for i, (low, high, hptr) in enumerate(hlist):
        h = hptr.GetValue()
        centers.append(0.5 * (low + high))
        halfwidths.append(0.5 * (high - low))

        c = ROOT.TCanvas(f"c_{bin_type}_{i}{name_suffix}", f"{bin_type}_{low}_{high}{name_suffix}", 800, 600)
        c_hybrid = ROOT.TCanvas(f"c_{bin_type}_{i}{name_suffix}_hybrid", f"{bin_type}_{low}_{high}{name_suffix}_hybrid", 800, 600)
        c_double = ROOT.TCanvas(f"c_{bin_type}_{i}{name_suffix}_double", f"{bin_type}_{low}_{high}{name_suffix}_double", 800, 600)
        c.SetGrid()
        c_hybrid.SetGrid()
        c_double.SetGrid()

        single_result = fit_single_gaussian(h, f"gauss_{bin_type}_{i}", res_min, res_max)
        double_result = fit_double_gaussian(h, f"double_gauss_{bin_type}_{i}", res_min, res_max)
        hybrid_result = choose_hybrid_result(single_result, double_result)

        entries = single_result["entries"]
        fit = single_result["fit"]
        mean = single_result["mean"]
        sigma = single_result["sigma"]
        mean_err = single_result["mean_err"]
        sigma_err = single_result["sigma_err"]
        chi2_ndf = single_result["chi2_ndf"]
        hybrid_fit = hybrid_result["fit"]
        dg_fit = double_result["fit"]

        means.append(mean)
        mean_errs.append(mean_err)
        sigmas.append(sigma)
        sigma_errs.append(sigma_err)
        chi2ndf_vals.append(chi2_ndf)

        hybrid_means.append(hybrid_result["mean"])
        hybrid_mean_errs.append(hybrid_result["mean_err"])
        hybrid_sigmas.append(hybrid_result["sigma"])
        hybrid_sigma_errs.append(hybrid_result["sigma_err"])
        hybrid_chi2ndf_vals.append(hybrid_result["chi2_ndf"])

        dg_means.append(double_result["mean"])
        dg_mean_errs.append(double_result["mean_err"])
        dg_sigma_eff.append(double_result["sigma_eff"])
        dg_sigma_eff_errs.append(double_result["sigma_eff_err"])
        dg_sigma_core.append(double_result["sigma_core"])
        dg_sigma_core_errs.append(double_result["sigma_core_err"])
        dg_sigma_tail.append(double_result["sigma_tail"])
        dg_sigma_tail_errs.append(double_result["sigma_tail_err"])
        dg_frac_core.append(double_result["frac_core"])
        dg_frac_core_errs.append(double_result["frac_core_err"])
        dg_chi2ndf_vals.append(double_result["chi2_ndf"])

        csv_rows.append([
            low, high, int(entries),
            mean, mean_err,
            sigma, sigma_err,
            single_result["chi2"], single_result["ndf"], chi2_ndf
        ])

        hybrid_csv_rows.append([
            low, high, int(entries),
            hybrid_result["mean"], hybrid_result["mean_err"],
            hybrid_result["sigma"], hybrid_result["sigma_err"],
            hybrid_result["chi2"], hybrid_result["ndf"], hybrid_result["chi2_ndf"],
            hybrid_result["fit_model"], hybrid_result["selection_reason"],
            single_result["fit_status"], double_result["fit_status"]
        ])

        double_csv_rows.append([
            low, high, int(entries),
            double_result["mean"], double_result["mean_err"],
            double_result["sigma_eff"], double_result["sigma_eff_err"],
            double_result["sigma_core"], double_result["sigma_core_err"],
            double_result["sigma_tail"], double_result["sigma_tail_err"],
            double_result["frac_core"], double_result["frac_core_err"],
            double_result["chi2"], double_result["ndf"], double_result["chi2_ndf"],
            double_result["fit_status"]
        ])

        h.SetStats(0)
        h.Draw()
        h.GetYaxis().SetTitle("Events")
        h.GetXaxis().SetTitle(residual_axis_label)

        if entries > 5:
            fit.Draw("same")
            leg = ROOT.TLegend(0.50, 0.46, 0.88, 0.84)
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            leg.AddEntry(0, f"Entries = {int(entries)}", "")
            leg.AddEntry(0, f"Mean = {mean:.4f} #pm {mean_err:.4f}", "")
            leg.AddEntry(0, f"#sigma = {sigma:.4f} #pm {sigma_err:.4f}", "")
            leg.AddEntry(0, f"#chi^2 = {single_result['chi2']:.4f}", "")
            leg.AddEntry(0, f"NDF = {int(single_result['ndf'])}", "")
            leg.AddEntry(0, f"#chi^2/NDF = {chi2_ndf:.4f}", "")
            leg.Draw()

        out_file.cd()
        h.Write()
        fit.Write(f"fit_{bin_type}_{int(low)}_{int(high)}{name_suffix}")
        c.SaveAs(os.path.join(base_dir, f"{bin_type}_Plots", f"{bin_type}_{muon_type}_{int(low)}_{int(high)}{name_suffix}.png"))

        canvas_all.cd(i + 1)
        h.Draw()
        if entries > 5:
            fit.Draw("same")

        c_hybrid.cd()
        h.Draw()
        h.GetYaxis().SetTitle("Events")
        h.GetXaxis().SetTitle(residual_axis_label)

        if hybrid_result["fit_model"] == "double" and entries > 10:
            hybrid_fit.Draw("same")
        elif entries > 5:
            hybrid_fit.Draw("same")

        leg_hybrid = ROOT.TLegend(0.42, 0.42, 0.88, 0.84)
        leg_hybrid.SetFillStyle(0)
        leg_hybrid.SetBorderSize(0)
        leg_hybrid.AddEntry(0, f"Entries = {int(entries)}", "")
        leg_hybrid.AddEntry(0, f"Model = {hybrid_result['fit_model']}", "")
        leg_hybrid.AddEntry(0, f"Reason = {hybrid_result['selection_reason']}", "")
        leg_hybrid.AddEntry(0, f"Mean = {hybrid_result['mean']:.4f} #pm {hybrid_result['mean_err']:.4f}", "")
        leg_hybrid.AddEntry(0, f"#sigma = {hybrid_result['sigma']:.4f} #pm {hybrid_result['sigma_err']:.4f}", "")
        leg_hybrid.AddEntry(0, f"#chi^2/NDF = {hybrid_result['chi2_ndf']:.4f}", "")
        leg_hybrid.Draw()

        hybrid_fit.Write(f"fit_hybrid_{bin_type}_{int(low)}_{int(high)}{name_suffix}")
        c_hybrid.SaveAs(os.path.join(base_dir, f"{bin_type}_Plots", f"{bin_type}_{muon_type}_{int(low)}_{int(high)}{name_suffix}_hybrid.png"))

        canvas_all_hybrid.cd(i + 1)
        h.Draw()
        if hybrid_result["fit_model"] == "double" and entries > 10:
            hybrid_fit.Draw("same")
        elif entries > 5:
            hybrid_fit.Draw("same")

        c_double.cd()
        h.Draw()
        h.GetYaxis().SetTitle("Events")
        h.GetXaxis().SetTitle(residual_axis_label)

        if entries > 10:
            dg_fit.Draw("same")
            leg_double = ROOT.TLegend(0.42, 0.42, 0.88, 0.84)
            leg_double.SetFillStyle(0)
            leg_double.SetBorderSize(0)
            leg_double.AddEntry(0, f"Entries = {int(entries)}", "")
            leg_double.AddEntry(0, f"Mean = {double_result['mean']:.4f} #pm {double_result['mean_err']:.4f}", "")
            leg_double.AddEntry(0, f"#sigma_{{eff}} = {double_result['sigma_eff']:.4f} #pm {double_result['sigma_eff_err']:.4f}", "")
            leg_double.AddEntry(0, f"#sigma_{{core}} = {double_result['sigma_core']:.4f} #pm {double_result['sigma_core_err']:.4f}", "")
            leg_double.AddEntry(0, f"#sigma_{{tail}} = {double_result['sigma_tail']:.4f} #pm {double_result['sigma_tail_err']:.4f}", "")
            leg_double.AddEntry(0, f"f_{{core}} = {double_result['frac_core']:.4f} #pm {double_result['frac_core_err']:.4f}", "")
            leg_double.AddEntry(0, f"#chi^2/NDF = {double_result['chi2_ndf']:.4f}", "")
            leg_double.Draw()

        dg_fit.Write(f"fit_double_{bin_type}_{int(low)}_{int(high)}{name_suffix}")
        c_double.SaveAs(os.path.join(base_dir, f"{bin_type}_Plots", f"{bin_type}_{muon_type}_{int(low)}_{int(high)}{name_suffix}_double.png"))

        canvas_all_double.cd(i + 1)
        h.Draw()
        if entries > 10:
            dg_fit.Draw("same")

    canvas_all.SaveAs(os.path.join(base_dir, f"{bin_type}_Plots", f"{muon_type}_{bin_type}_all_fits{name_suffix}.png"))
    canvas_all_hybrid.SaveAs(os.path.join(base_dir, f"{bin_type}_Plots", f"{muon_type}_{bin_type}_all_fits{name_suffix}_hybrid.png"))
    canvas_all_double.SaveAs(os.path.join(base_dir, f"{bin_type}_Plots", f"{muon_type}_{bin_type}_all_fits{name_suffix}_double.png"))

    g_mean = ROOT.TGraphErrors(
        len(hlist),
        array_to_c(centers),
        array_to_c(means),
        array_to_c(halfwidths),
        array_to_c(mean_errs),
    )
    g_mean.SetName(f"{bin_type}_mean_total{name_suffix}")
    g_mean.SetTitle(f"{bin_type} mean;{x_title};Mean of {residual_summary_label}")
    g_mean.SetMarkerStyle(8)
    g_mean.SetLineWidth(2)
    draw_graph_and_save(
        g_mean,
        f"c_mean_{bin_type}",
        out_file,
        os.path.join(base_dir, f"{bin_type}_Plots", f"{bin_type}_mean_total{name_suffix}.png"),
        logx
    )

    g_sigma = ROOT.TGraphErrors(
        len(hlist),
        array_to_c(centers),
        array_to_c(sigmas),
        array_to_c(halfwidths),
        array_to_c(sigma_errs),
    )
    g_sigma.SetName(f"{bin_type}_sigma_total{name_suffix}")
    g_sigma.SetTitle(f"{bin_type} sigma;{x_title};#sigma of {residual_summary_label}")
    g_sigma.SetMarkerStyle(8)
    g_sigma.SetLineWidth(2)
    draw_graph_and_save(
        g_sigma,
        f"c_sigma_{bin_type}",
        out_file,
        os.path.join(base_dir, f"{bin_type}_Plots", f"{bin_type}_sigma_total{name_suffix}.png"),
        logx
    )

    g_chi2ndf = ROOT.TGraphErrors(
        len(hlist),
        array_to_c(centers),
        array_to_c(chi2ndf_vals),
        array_to_c(halfwidths),
        array_to_c([0.0] * len(hlist)),
    )
    g_chi2ndf.SetName(f"{bin_type}_chi2ndf_total{name_suffix}")
    g_chi2ndf.SetTitle(f"{bin_type} #chi^{{2}}/NDF;{x_title};#chi^{{2}}/NDF")
    g_chi2ndf.SetMarkerStyle(8)
    g_chi2ndf.SetLineWidth(2)
    draw_graph_and_save(
        g_chi2ndf,
        f"c_chi2ndf_{bin_type}",
        out_file,
        os.path.join(base_dir, f"{bin_type}_Plots", f"{bin_type}_chi2ndf_total{name_suffix}.png"),
        logx
    )

    g_hybrid_mean = ROOT.TGraphErrors(
        len(hlist),
        array_to_c(centers),
        array_to_c(hybrid_means),
        array_to_c(halfwidths),
        array_to_c(hybrid_mean_errs),
    )
    g_hybrid_mean.SetName(f"{bin_type}_mean_total{name_suffix}_hybrid")
    g_hybrid_mean.SetTitle(f"{bin_type} hybrid mean;{x_title};Mean of {residual_summary_label}")
    g_hybrid_mean.SetMarkerStyle(29)
    g_hybrid_mean.SetMarkerColor(ROOT.kOrange + 7)
    g_hybrid_mean.SetLineColor(ROOT.kOrange + 7)
    g_hybrid_mean.SetLineWidth(2)
    draw_graph_and_save(
        g_hybrid_mean,
        f"c_mean_{bin_type}_hybrid",
        out_file,
        os.path.join(base_dir, f"{bin_type}_Plots", f"{bin_type}_mean_total{name_suffix}_hybrid.png"),
        logx
    )

    g_hybrid_sigma = ROOT.TGraphErrors(
        len(hlist),
        array_to_c(centers),
        array_to_c(hybrid_sigmas),
        array_to_c(halfwidths),
        array_to_c(hybrid_sigma_errs),
    )
    g_hybrid_sigma.SetName(f"{bin_type}_sigma_total{name_suffix}_hybrid")
    g_hybrid_sigma.SetTitle(f"{bin_type} hybrid sigma;{x_title};#sigma of {residual_summary_label}")
    g_hybrid_sigma.SetMarkerStyle(29)
    g_hybrid_sigma.SetMarkerColor(ROOT.kOrange + 7)
    g_hybrid_sigma.SetLineColor(ROOT.kOrange + 7)
    g_hybrid_sigma.SetLineWidth(2)
    draw_graph_and_save(
        g_hybrid_sigma,
        f"c_sigma_{bin_type}_hybrid",
        out_file,
        os.path.join(base_dir, f"{bin_type}_Plots", f"{bin_type}_sigma_total{name_suffix}_hybrid.png"),
        logx
    )

    g_hybrid_chi2ndf = ROOT.TGraphErrors(
        len(hlist),
        array_to_c(centers),
        array_to_c(hybrid_chi2ndf_vals),
        array_to_c(halfwidths),
        array_to_c([0.0] * len(hlist)),
    )
    g_hybrid_chi2ndf.SetName(f"{bin_type}_chi2ndf_total{name_suffix}_hybrid")
    g_hybrid_chi2ndf.SetTitle(f"{bin_type} hybrid #chi^{{2}}/NDF;{x_title};#chi^{{2}}/NDF")
    g_hybrid_chi2ndf.SetMarkerStyle(29)
    g_hybrid_chi2ndf.SetMarkerColor(ROOT.kOrange + 7)
    g_hybrid_chi2ndf.SetLineColor(ROOT.kOrange + 7)
    g_hybrid_chi2ndf.SetLineWidth(2)
    draw_graph_and_save(
        g_hybrid_chi2ndf,
        f"c_chi2ndf_{bin_type}_hybrid",
        out_file,
        os.path.join(base_dir, f"{bin_type}_Plots", f"{bin_type}_chi2ndf_total{name_suffix}_hybrid.png"),
        logx
    )

    g_double_mean = ROOT.TGraphErrors(
        len(hlist),
        array_to_c(centers),
        array_to_c(dg_means),
        array_to_c(halfwidths),
        array_to_c(dg_mean_errs),
    )
    g_double_mean.SetName(f"{bin_type}_mean_total{name_suffix}_double")
    g_double_mean.SetTitle(f"{bin_type} double-gaussian mean;{x_title};Mean of {residual_summary_label}")
    g_double_mean.SetMarkerStyle(22)
    g_double_mean.SetMarkerColor(ROOT.kBlue + 1)
    g_double_mean.SetLineColor(ROOT.kBlue + 1)
    g_double_mean.SetLineWidth(2)
    draw_graph_and_save(
        g_double_mean,
        f"c_mean_{bin_type}_double",
        out_file,
        os.path.join(base_dir, f"{bin_type}_Plots", f"{bin_type}_mean_total{name_suffix}_double.png"),
        logx
    )

    g_double_sigma_eff = ROOT.TGraphErrors(
        len(hlist),
        array_to_c(centers),
        array_to_c(dg_sigma_eff),
        array_to_c(halfwidths),
        array_to_c(dg_sigma_eff_errs),
    )
    g_double_sigma_eff.SetName(f"{bin_type}_sigma_eff_total{name_suffix}_double")
    g_double_sigma_eff.SetTitle(f"{bin_type} double-gaussian effective sigma;{x_title};#sigma_{{eff}} of {residual_summary_label}")
    g_double_sigma_eff.SetMarkerStyle(22)
    g_double_sigma_eff.SetMarkerColor(ROOT.kBlue + 1)
    g_double_sigma_eff.SetLineColor(ROOT.kBlue + 1)
    g_double_sigma_eff.SetLineWidth(2)
    draw_graph_and_save(
        g_double_sigma_eff,
        f"c_sigmaeff_{bin_type}_double",
        out_file,
        os.path.join(base_dir, f"{bin_type}_Plots", f"{bin_type}_sigma_eff_total{name_suffix}_double.png"),
        logx
    )

    g_double_sigma_core = ROOT.TGraphErrors(
        len(hlist),
        array_to_c(centers),
        array_to_c(dg_sigma_core),
        array_to_c(halfwidths),
        array_to_c(dg_sigma_core_errs),
    )
    g_double_sigma_core.SetName(f"{bin_type}_sigma_core_total{name_suffix}_double")
    g_double_sigma_core.SetTitle(f"{bin_type} double-gaussian core sigma;{x_title};#sigma_{{core}} of {residual_summary_label}")
    g_double_sigma_core.SetMarkerStyle(23)
    g_double_sigma_core.SetMarkerColor(ROOT.kGreen + 2)
    g_double_sigma_core.SetLineColor(ROOT.kGreen + 2)
    g_double_sigma_core.SetLineWidth(2)
    draw_graph_and_save(
        g_double_sigma_core,
        f"c_sigmacore_{bin_type}_double",
        out_file,
        os.path.join(base_dir, f"{bin_type}_Plots", f"{bin_type}_sigma_core_total{name_suffix}_double.png"),
        logx
    )

    g_double_sigma_tail = ROOT.TGraphErrors(
        len(hlist),
        array_to_c(centers),
        array_to_c(dg_sigma_tail),
        array_to_c(halfwidths),
        array_to_c(dg_sigma_tail_errs),
    )
    g_double_sigma_tail.SetName(f"{bin_type}_sigma_tail_total{name_suffix}_double")
    g_double_sigma_tail.SetTitle(f"{bin_type} double-gaussian tail sigma;{x_title};#sigma_{{tail}} of {residual_summary_label}")
    g_double_sigma_tail.SetMarkerStyle(21)
    g_double_sigma_tail.SetMarkerColor(ROOT.kRed + 1)
    g_double_sigma_tail.SetLineColor(ROOT.kRed + 1)
    g_double_sigma_tail.SetLineWidth(2)
    draw_graph_and_save(
        g_double_sigma_tail,
        f"c_sigmatail_{bin_type}_double",
        out_file,
        os.path.join(base_dir, f"{bin_type}_Plots", f"{bin_type}_sigma_tail_total{name_suffix}_double.png"),
        logx
    )

    g_double_frac_core = ROOT.TGraphErrors(
        len(hlist),
        array_to_c(centers),
        array_to_c(dg_frac_core),
        array_to_c(halfwidths),
        array_to_c(dg_frac_core_errs),
    )
    g_double_frac_core.SetName(f"{bin_type}_frac_core_total{name_suffix}_double")
    g_double_frac_core.SetTitle(f"{bin_type} double-gaussian core fraction;{x_title};f_{{core}}")
    g_double_frac_core.SetMarkerStyle(20)
    g_double_frac_core.SetMarkerColor(ROOT.kMagenta + 2)
    g_double_frac_core.SetLineColor(ROOT.kMagenta + 2)
    g_double_frac_core.SetLineWidth(2)
    draw_graph_and_save(
        g_double_frac_core,
        f"c_fraccore_{bin_type}_double",
        out_file,
        os.path.join(base_dir, f"{bin_type}_Plots", f"{bin_type}_frac_core_total{name_suffix}_double.png"),
        logx
    )

    g_double_chi2ndf = ROOT.TGraphErrors(
        len(hlist),
        array_to_c(centers),
        array_to_c(dg_chi2ndf_vals),
        array_to_c(halfwidths),
        array_to_c([0.0] * len(hlist)),
    )
    g_double_chi2ndf.SetName(f"{bin_type}_chi2ndf_total{name_suffix}_double")
    g_double_chi2ndf.SetTitle(f"{bin_type} double-gaussian #chi^{{2}}/NDF;{x_title};#chi^{{2}}/NDF")
    g_double_chi2ndf.SetMarkerStyle(22)
    g_double_chi2ndf.SetMarkerColor(ROOT.kBlue + 1)
    g_double_chi2ndf.SetLineColor(ROOT.kBlue + 1)
    g_double_chi2ndf.SetLineWidth(2)
    draw_graph_and_save(
        g_double_chi2ndf,
        f"c_chi2ndf_{bin_type}_double",
        out_file,
        os.path.join(base_dir, f"{bin_type}_Plots", f"{bin_type}_chi2ndf_total{name_suffix}_double.png"),
        logx
    )

    csv_path = os.path.join(base_dir, f"{bin_type}_Plots", f"{bin_type}_fit_summary{name_suffix}.csv")
    write_fit_csv(csv_path, csv_rows)
    hybrid_csv_path = os.path.join(base_dir, f"{bin_type}_Plots", f"{bin_type}_fit_summary{name_suffix}_hybrid.csv")
    write_hybrid_fit_csv(hybrid_csv_path, hybrid_csv_rows)
    double_csv_path = os.path.join(base_dir, f"{bin_type}_Plots", f"{bin_type}_fit_summary{name_suffix}_double.csv")
    write_double_fit_csv(double_csv_path, double_csv_rows)


def make_overlay_plot(h_tag, h_probe, out_file, out_png, name="overlay"):
    c = ROOT.TCanvas(f"c_{name}", "", 800, 600)
    c.SetGrid()

    h_tag = h_tag.Clone(f"{name}_tag_clone")
    h_probe = h_probe.Clone(f"{name}_probe_clone")

    h_tag.SetLineColor(ROOT.kRed + 1)
    h_probe.SetLineColor(ROOT.kBlue + 1)
    h_tag.SetLineWidth(2)
    h_probe.SetLineWidth(2)
    h_tag.SetStats(0)
    h_probe.SetStats(0)

    max_y = max(h_tag.GetMaximum(), h_probe.GetMaximum())
    h_tag.SetMaximum(1.2 * max_y if max_y > 0 else 1.0)

    h_tag.Draw("hist")
    h_probe.Draw("hist same")

    leg = ROOT.TLegend(0.65, 0.75, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(h_tag, "Tag", "l")
    leg.AddEntry(h_probe, "Probe", "l")
    leg.Draw()

    out_file.cd()
    c.Write(f"c_{name}")
    c.SaveAs(out_png)


# ============================================================
# Main
# ============================================================
def main():
    parser = argparse.ArgumentParser(description="Tag-and-probe cosmic muon analysis with PyROOT RDataFrame")

    parser.add_argument("-i", "--input", nargs="+", required=True, help="Input ROOT files, directories, or glob patterns")
    parser.add_argument("--tree-name", default="Events", help="TTree name")
    parser.add_argument("--datatype", choices=["MC", "DATA"], default="MC", help="Sample type")
    parser.add_argument("--muon-type", choices=["DSA", "DGL"], default="DSA", help="Muon type")
    parser.add_argument("--output", default=None, help="Output ROOT file name")
    parser.add_argument("--outdir", default=".", help="Base output directory")
    parser.add_argument("--pt-bins", type=parse_bin_list, default=DEFAULT_PT_BINS, help="Comma-separated pT bin edges")
    parser.add_argument("--dz-bins", type=parse_bin_list, default=DEFAULT_DZ_BINS, help="Comma-separated |dz| bin edges")
    parser.add_argument("--dxy-bins", type=parse_bin_list, default=DEFAULT_DXY_BINS, help="Comma-separated |dxy| bin edges")
    parser.add_argument("--min-pt", type=float, default=12.5, help="Minimum pT cut")
    parser.add_argument("--threads", type=int, default=0, help="Number of threads for implicit MT. 0 = ROOT default")
    parser.add_argument("--no-trigger", action="store_true", help="Ignore HLT requirement for both DATA and MC")
    parser.add_argument("--require-trigger-for-mc", action="store_true", help="Apply the same HLT requirement to MC samples")
    parser.add_argument("--res-range", type=parse_bin_list, default=None, help='Residual plot range as "xmin,xmax"')
    parser.add_argument("--max-rel-pt-err-tag", type=float, default=None, help="Maximum allowed ptError/pt for tag muon")
    parser.add_argument("--max-rel-pt-err-probe", type=float, default=None, help="Maximum allowed ptError/pt for probe muon")
    parser.add_argument("--max-normalized-chi2-tag", type=float, default=None, help="Maximum allowed normalized chi2 for tag muon")
    parser.add_argument("--max-normalized-chi2-probe", type=float, default=None, help="Maximum allowed normalized chi2 for probe muon")
    parser.add_argument("--min-primary-hit-count-tag", type=int, default=None, help="Minimum allowed primary hit count for tag muon")
    parser.add_argument("--min-primary-hit-count-probe", type=int, default=None, help="Minimum allowed primary hit count for probe muon")
    parser.add_argument("--min-secondary-hit-count-tag", type=int, default=None, help="Minimum allowed secondary hit count for tag muon")
    parser.add_argument("--min-secondary-hit-count-probe", type=int, default=None, help="Minimum allowed secondary hit count for probe muon")

    args = parser.parse_args()

    if args.threads > 0:
        ROOT.ROOT.EnableImplicitMT(args.threads)
    else:
        ROOT.ROOT.EnableImplicitMT()

    files = expand_inputs(args.input)
    if not files:
        raise RuntimeError("No ROOT files found from --input")

    output_name = args.output or f"Cosmics_muons_{args.datatype}_{args.muon_type}_RDF.root"

    base_dir = os.path.join(args.outdir, args.datatype, args.muon_type)
    ensure_dir(os.path.join(base_dir, "pt_Plots"))
    ensure_dir(os.path.join(base_dir, "dz_Plots"))
    ensure_dir(os.path.join(base_dir, "dxy_Plots"))
    ensure_dir(os.path.join(base_dir, "Control_plots"))

    print(f"Found {len(files)} ROOT files")
    print(f"Running with datatype={args.datatype}, muon_type={args.muon_type}")

    max_rel_pt_err_tag = args.max_rel_pt_err_tag if args.max_rel_pt_err_tag is not None else -1.0
    max_rel_pt_err_probe = args.max_rel_pt_err_probe if args.max_rel_pt_err_probe is not None else -1.0
    max_normalized_chi2_tag = args.max_normalized_chi2_tag if args.max_normalized_chi2_tag is not None else -1.0
    max_normalized_chi2_probe = args.max_normalized_chi2_probe if args.max_normalized_chi2_probe is not None else -1.0
    min_primary_hit_count_tag = args.min_primary_hit_count_tag if args.min_primary_hit_count_tag is not None else -1
    min_primary_hit_count_probe = args.min_primary_hit_count_probe if args.min_primary_hit_count_probe is not None else -1
    min_secondary_hit_count_tag = args.min_secondary_hit_count_tag if args.min_secondary_hit_count_tag is not None else -1
    min_secondary_hit_count_probe = args.min_secondary_hit_count_probe if args.min_secondary_hit_count_probe is not None else -1

    if max_rel_pt_err_tag >= 0:
        print(f"Applying tag ptError/pt cut   < {max_rel_pt_err_tag}")
    if max_rel_pt_err_probe >= 0:
        print(f"Applying probe ptError/pt cut < {max_rel_pt_err_probe}")
    if max_normalized_chi2_tag >= 0:
        print(f"Applying tag normalized chi2 cut   < {max_normalized_chi2_tag}")
    if max_normalized_chi2_probe >= 0:
        print(f"Applying probe normalized chi2 cut < {max_normalized_chi2_probe}")
    if min_primary_hit_count_tag >= 0:
        print(f"Applying tag primary-hit cut   >= {min_primary_hit_count_tag}")
    if min_primary_hit_count_probe >= 0:
        print(f"Applying probe primary-hit cut >= {min_primary_hit_count_probe}")
    if min_secondary_hit_count_tag >= 0:
        print(f"Applying tag secondary-hit cut   >= {min_secondary_hit_count_tag}")
    if min_secondary_hit_count_probe >= 0:
        print(f"Applying probe secondary-hit cut >= {min_secondary_hit_count_probe}")

    file_vec = ROOT.std.vector("string")()
    for f in files:
        file_vec.push_back(f)

    df = ROOT.RDataFrame(args.tree_name, file_vec)

    require_trigger = (not args.no_trigger) and (args.datatype == "DATA" or args.require_trigger_for_mc)
    require_trigger_cpp = "true" if require_trigger else "false"

    if require_trigger:
        print("Applying HLT_L2Mu10_NoVertex_NoBPTX3BX trigger requirement")
    else:
        print("Trigger requirement disabled for this sample")

    if args.muon_type == "DSA":
        pt_branch = "dmu_dsa_pt"
        pterr_branch = "dmu_dsa_ptError"
        chi2_branch = "dmu_dsa_normalizedChi2"
        primary_hit_count_branch = "dmu_dsa_nValidMuonDTHits"
        secondary_hit_count_branch = "dmu_dsa_nValidStripHits"
        eta_branch = "dmu_dsa_eta"
        phi_branch = "dmu_dsa_phi"
        dz_branch = "dmu_dsa_dz"
        dxy_branch = "dmu_dsa_dxy"
        charge_branch = "dmu_dsa_charge"
        pass_tag_branch = "dmu_dsa_passTagID"
        probe_id_branch = "dmu_dsa_probeID"
        has_probe_branch = "dmu_dsa_hasProbe"
        reco_type_branch = "dmu_isDSA"
        default_res_min, default_res_max = -5.0, 5.0
        all_muon_label = "DSA"
    else:
        pt_branch = "dmu_dgl_pt"
        pterr_branch = "dmu_dgl_ptError"
        chi2_branch = "dmu_dgl_normalizedChi2"
        primary_hit_count_branch = "dmu_dgl_nMuonHits"
        secondary_hit_count_branch = "dmu_dgl_nValidStripHits"
        eta_branch = "dmu_dgl_eta"
        phi_branch = "dmu_dgl_phi"
        dz_branch = "dmu_dgl_dz"
        dxy_branch = "dmu_dgl_dxy"
        charge_branch = "dmu_dgl_charge"
        pass_tag_branch = "dmu_dgl_passTagID"
        probe_id_branch = "dmu_dgl_probeID"
        has_probe_branch = "dmu_dgl_hasProbe"
        reco_type_branch = "dmu_isDGL"
        default_res_min, default_res_max = -0.3, 0.3
        all_muon_label = "DGL"

    df1 = df.Define(
        "ana",
        f"""
        analyze_event(
            ndmu,
            HLT_L2Mu10_NoVertex_NoBPTX3BX,
            {pt_branch},
            {pterr_branch},
            {chi2_branch},
            {primary_hit_count_branch},
            {secondary_hit_count_branch},
            {eta_branch},
            {phi_branch},
            {dz_branch},
            {dxy_branch},
            {charge_branch},
            {pass_tag_branch},
            {probe_id_branch},
            {has_probe_branch},
            {reco_type_branch},
            {require_trigger_cpp},
            {args.min_pt},
            {max_rel_pt_err_tag},
            {max_rel_pt_err_probe},
            {max_normalized_chi2_tag},
            {max_normalized_chi2_probe},
            {min_primary_hit_count_tag},
            {min_primary_hit_count_probe},
            {min_secondary_hit_count_tag},
            {min_secondary_hit_count_probe}
        )
        """
    )

    if args.res_range is not None:
        if len(args.res_range) != 2:
            raise ValueError("--res-range must contain exactly two values")
        res_min, res_max = args.res_range
    else:
        res_min, res_max = default_res_min, default_res_max

    df2 = (
        df1.Define("resolutions", "ana.resolutions")
           .Define("symmetric_resolutions", "ana.symmetric_resolutions")
           .Define("tag_pt", "ana.tag_pt")
           .Define("tag_eta", "ana.tag_eta")
           .Define("tag_phi", "ana.tag_phi")
           .Define("tag_charge", "ana.tag_charge")
           .Define("tag_pt_err_over_pt", "ana.tag_pt_err_over_pt")
           .Define("probe_pt", "ana.probe_pt")
           .Define("probe_eta", "ana.probe_eta")
           .Define("probe_phi", "ana.probe_phi")
           .Define("probe_charge", "ana.probe_charge")
           .Define("probe_pt_err_over_pt", "ana.probe_pt_err_over_pt")
           .Define("all_muon_pt", "ana.all_muon_pt")
           .Define("tag_dz_abs", "ana.tag_dz_abs")
           .Define("tag_dxy_abs", "ana.tag_dxy_abs")
    )

    h_total = df2.Histo1D(
        ("Tot_pthist", "q/pT Residual distribution;q/p_{T} residual;Events", 80, res_min, res_max),
        "resolutions"
    )
    h_total_sym = df2.Histo1D(
        ("Tot_pthist_sym", "Symmetric q/pT Residual distribution;symmetric q/p_{T} residual;Events", 80, res_min, res_max),
        "symmetric_resolutions"
    )

    h_pt_tag = df2.Histo1D(("hist_pt_tag", "Tag muon p_{T};p_{T} [GeV];Events", 350, 0, 300), "tag_pt")
    h_eta_tag = df2.Histo1D(("hist_eta_tag", "Tag muon #eta;#eta;Events", 60, -3, 3), "tag_eta")
    h_phi_tag = df2.Histo1D(("hist_phi_tag", "Tag muon #phi;#phi;Events", 64, -3.2, 3.2), "tag_phi")
    h_charge_tag = df2.Histo1D(("hist_charge_tag", "Tag muon charge;charge;Events", 10, -4, 4), "tag_charge")

    h_pt_probe = df2.Histo1D(("hist_pt_probe", "Probe muon p_{T};p_{T} [GeV];Events", 350, 0, 300), "probe_pt")
    h_eta_probe = df2.Histo1D(("hist_eta_probe", "Probe muon #eta;#eta;Events", 60, -3, 3), "probe_eta")
    h_phi_probe = df2.Histo1D(("hist_phi_probe", "Probe muon #phi;#phi;Events", 64, -3.2, 3.2), "probe_phi")
    h_charge_probe = df2.Histo1D(("hist_charge_probe", "Probe muon charge;charge;Events", 10, -4, 4), "probe_charge")

    h_tag_dz_abs = df2.Histo1D(
        ("hist_tag_dz_abs", "Tag |dz|;|dz|;Events", 100, 0.0, 150.0),
        "tag_dz_abs"
    )

    h_tag_dxy_abs = df2.Histo1D(
        ("hist_tag_dxy_abs", "Tag |dxy|;|dxy|;Events", 100, 0.0, 80.0),
        "tag_dxy_abs"
    )

    h_tag_pt_err_over_pt = df2.Histo1D(
        ("hist_tag_ptErrOverPt", "Tag p_{T}^{error}/p_{T};p_{T}^{error}/p_{T};Events", 100, 0.0, 1.0),
        "tag_pt_err_over_pt"
    )
    h_probe_pt_err_over_pt = df2.Histo1D(
        ("hist_probe_ptErrOverPt", "Probe p_{T}^{error}/p_{T};p_{T}^{error}/p_{T};Events", 100, 0.0, 1.0),
        "probe_pt_err_over_pt"
    )

    h_all_muon_pt = df2.Histo1D(
        (f"hist_pt_{all_muon_label}", f"{all_muon_label} muon p_{{T}};p_{{T}} [GeV];Events", 350, 0, 300),
        "all_muon_pt"
    )

    df_bins = df2
    pt_histos = []
    pt_histos_sym = []
    dz_histos = []
    dz_histos_sym = []
    dxy_histos = []
    dxy_histos_sym = []

    for i in range(len(args.pt_bins) - 1):
        low, high = args.pt_bins[i], args.pt_bins[i + 1]
        cname = f"res_pt_bin_{i}"
        cname_sym = f"res_pt_bin_{i}_sym"
        df_bins = df_bins.Define(cname, f"filter_by_bin(resolutions, tag_pt, {low}, {high})")
        df_bins = df_bins.Define(cname_sym, f"filter_by_bin(symmetric_resolutions, tag_pt, {low}, {high})")
        h = df_bins.Histo1D(
            (f"pt_{int(low)}_{int(high)}", f"pT {low}-{high};q/p_{{T}} residual;Events", 100, res_min, res_max),
            cname
        )
        h_sym = df_bins.Histo1D(
            (f"pt_{int(low)}_{int(high)}_sym", f"pT {low}-{high};symmetric q/p_{{T}} residual;Events", 100, res_min, res_max),
            cname_sym
        )
        pt_histos.append((low, high, h))
        pt_histos_sym.append((low, high, h_sym))

    for i in range(len(args.dz_bins) - 1):
        low, high = args.dz_bins[i], args.dz_bins[i + 1]
        cname = f"res_dz_bin_{i}"
        cname_sym = f"res_dz_bin_{i}_sym"
        df_bins = df_bins.Define(cname, f"filter_by_bin(resolutions, tag_dz_abs, {low}, {high})")
        df_bins = df_bins.Define(cname_sym, f"filter_by_bin(symmetric_resolutions, tag_dz_abs, {low}, {high})")
        h = df_bins.Histo1D(
            (f"dz_{int(low)}_{int(high)}", f"|dz| {low}-{high};q/p_{{T}} residual;Events", 100, res_min, res_max),
            cname
        )
        h_sym = df_bins.Histo1D(
            (f"dz_{int(low)}_{int(high)}_sym", f"|dz| {low}-{high};symmetric q/p_{{T}} residual;Events", 100, res_min, res_max),
            cname_sym
        )
        dz_histos.append((low, high, h))
        dz_histos_sym.append((low, high, h_sym))

    for i in range(len(args.dxy_bins) - 1):
        low, high = args.dxy_bins[i], args.dxy_bins[i + 1]
        cname = f"res_dxy_bin_{i}"
        cname_sym = f"res_dxy_bin_{i}_sym"
        df_bins = df_bins.Define(cname, f"filter_by_bin(resolutions, tag_dxy_abs, {low}, {high})")
        df_bins = df_bins.Define(cname_sym, f"filter_by_bin(symmetric_resolutions, tag_dxy_abs, {low}, {high})")
        h = df_bins.Histo1D(
            (f"dxy_{int(low)}_{int(high)}", f"|dxy| {low}-{high};q/p_{{T}} residual;Events", 100, res_min, res_max),
            cname
        )
        h_sym = df_bins.Histo1D(
            (f"dxy_{int(low)}_{int(high)}_sym", f"|dxy| {low}-{high};symmetric q/p_{{T}} residual;Events", 100, res_min, res_max),
            cname_sym
        )
        dxy_histos.append((low, high, h))
        dxy_histos_sym.append((low, high, h_sym))

    out_path = os.path.join(args.outdir, output_name)
    out = ROOT.TFile(out_path, "RECREATE")

    main_hists = [
        h_total_sym,
        h_total,
        h_pt_tag, h_eta_tag, h_phi_tag, h_charge_tag,
        h_pt_probe, h_eta_probe, h_phi_probe, h_charge_probe,
        h_tag_dz_abs, h_tag_dxy_abs,
        h_tag_pt_err_over_pt, h_probe_pt_err_over_pt,
        h_all_muon_pt
    ]

    for hptr in main_hists:
        h = hptr.GetValue()
        out.cd()
        h.Write()

    control_specs = [
        (h_total_sym.GetValue(), "Tot_pthist_sym.png"),
        (h_pt_tag.GetValue(), "hist_pt_tag.png"),
        (h_eta_tag.GetValue(), "hist_eta_tag.png"),
        (h_phi_tag.GetValue(), "hist_phi_tag.png"),
        (h_charge_tag.GetValue(), "hist_charge_tag.png"),
        (h_pt_probe.GetValue(), "hist_pt_probe.png"),
        (h_eta_probe.GetValue(), "hist_eta_probe.png"),
        (h_phi_probe.GetValue(), "hist_phi_probe.png"),
        (h_charge_probe.GetValue(), "hist_charge_probe.png"),
        (h_tag_dz_abs.GetValue(), "hist_tag_dz_abs.png"),
        (h_tag_dxy_abs.GetValue(), "hist_tag_dxy_abs.png"),
        (h_tag_pt_err_over_pt.GetValue(), "hist_tag_ptErrOverPt.png"),
        (h_probe_pt_err_over_pt.GetValue(), "hist_probe_ptErrOverPt.png"),
        (h_all_muon_pt.GetValue(), f"hist_pT_{all_muon_label}.png"),
    ]

    for hist, png_name in control_specs:
        c = ROOT.TCanvas(f"c_{png_name}", "", 800, 600)
        hist.SetStats(0)
        hist.Draw()
        c.SaveAs(os.path.join(base_dir, "Control_plots", png_name))

    make_overlay_plot(
        h_tag_pt_err_over_pt.GetValue(),
        h_probe_pt_err_over_pt.GetValue(),
        out,
        os.path.join(base_dir, "Control_plots", "ptErr_over_pt_comparison.png"),
        name="ptErrOverPt_compare"
    )

    make_overlay_plot(
        h_pt_tag.GetValue(),
        h_pt_probe.GetValue(),
        out,
        os.path.join(base_dir, "Control_plots", "pt_tag_probe_comparison.png"),
        name="pt_tag_probe_compare"
    )

    fit_and_draw(pt_histos, "pt", "p_{T}^{tag} [GeV]", base_dir, args.muon_type, res_min, res_max, out, logx=True)
    fit_and_draw(dz_histos, "dz", "|dz|", base_dir, args.muon_type, res_min, res_max, out, logx=True)
    fit_and_draw(dxy_histos, "dxy", "|dxy|", base_dir, args.muon_type, res_min, res_max, out, logx=True)
    fit_and_draw(
        pt_histos_sym,
        "pt",
        "p_{T}^{tag} [GeV]",
        base_dir,
        args.muon_type,
        res_min,
        res_max,
        out,
        logx=True,
        name_suffix="_sym",
        residual_axis_label="symmetric q/p_{T} residual",
        residual_summary_label="symmetric q/p_{T} relative residual",
    )
    fit_and_draw(
        dz_histos_sym,
        "dz",
        "|dz|",
        base_dir,
        args.muon_type,
        res_min,
        res_max,
        out,
        logx=True,
        name_suffix="_sym",
        residual_axis_label="symmetric q/p_{T} residual",
        residual_summary_label="symmetric q/p_{T} relative residual",
    )
    fit_and_draw(
        dxy_histos_sym,
        "dxy",
        "|dxy|",
        base_dir,
        args.muon_type,
        res_min,
        res_max,
        out,
        logx=True,
        name_suffix="_sym",
        residual_axis_label="symmetric q/p_{T} residual",
        residual_summary_label="symmetric q/p_{T} relative residual",
    )

    out.Close()
    print(f"Done. Output written to: {out_path}")


if __name__ == "__main__":
    main()
