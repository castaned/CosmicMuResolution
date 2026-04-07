#!/usr/bin/env python3

import ROOT
import os
import glob
import argparse
import array
import csv

ROOT.gStyle.SetOptStat(0)

# ============================================================
# Defaults
# ============================================================
DEFAULT_PT_BINS = [20.0, 30.0, 40.0, 50.0, 65.0, 85.0, 120.0, 200.0, 400.0, 1000.0]
DEFAULT_DZ_BINS = [1.0,5.0, 10.0, 20.0, 30.0, 45.0, 60.0,100,150.0]
DEFAULT_DXY_BINS = [1.0,5.0, 10.0, 20.0, 30.0, 40.0,60.0, 80.0]


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
    const RVec<float>& eta,
    const RVec<float>& phi,
    const RVec<float>& dz,
    const RVec<float>& dxy,
    const RVec<float>& charge,
    const RVec<char>& passTagID,
    const RVec<int>& probeID,
    const RVec<char>& hasProbe,
    const RVec<int>& isRecoType,
    bool isMC,
    double min_pt,
    double max_rel_pt_err_tag,
    double max_rel_pt_err_probe
) {
    AnalysisOutput out;

    if (!isMC && !hlt) return out;
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

        if (max_rel_pt_err_tag >= 0.0 && tag_rel_pt_err >= max_rel_pt_err_tag) continue;
        if (max_rel_pt_err_probe >= 0.0 && probe_rel_pt_err >= max_rel_pt_err_probe) continue;

        const double inv_up = std::abs(charge[j] / pt[j]);
        const double inv_down = std::abs(charge[i] / pt[i]);
        if (inv_down == 0.) continue;

        const double resolution = (inv_up - inv_down) / (inv_down);

        out.resolutions.push_back(resolution);

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


# ============================================================
# Plot and fit
# ============================================================
def fit_and_draw(hlist, bin_type, x_title, base_dir, muon_type, res_min, res_max, out_file, logx=True):
    means, mean_errs = [], []
    sigmas, sigma_errs = [], []
    chi2_vals, ndf_vals, chi2ndf_vals = [], [], []
    centers, halfwidths = [], []
    csv_rows = []

    canvas_all = ROOT.TCanvas(f"All_plots_{bin_type}", f"{bin_type} Histograms", 1500, 1000)
    n = len(hlist)
    nx = 3
    ny = (n + nx - 1) // nx
    canvas_all.Divide(nx, ny)

    for i, (low, high, hptr) in enumerate(hlist):
        h = hptr.GetValue()
        centers.append(0.5 * (low + high))
        halfwidths.append(0.5 * (high - low))

        c = ROOT.TCanvas(f"c_{bin_type}_{i}", f"{bin_type}_{low}_{high}", 800, 600)
        c.SetGrid()

        fit = ROOT.TF1(f"gauss_{bin_type}_{i}", "gaus", res_min, res_max)

        entries = h.GetEntries()
        mean = sigma = mean_err = sigma_err = 0.0
        chi2 = chi2_ndf = 0.0
        ndf = 0

        if entries > 5:
            h.Fit(fit, "QSR")
            mean = fit.GetParameter(1)
            sigma = fit.GetParameter(2)
            mean_err = fit.GetParError(1)
            sigma_err = fit.GetParError(2)
            chi2 = fit.GetChisquare()
            ndf = fit.GetNDF()
            chi2_ndf = chi2 / ndf if ndf > 0 else 0.0

        means.append(mean)
        mean_errs.append(mean_err)
        sigmas.append(sigma)
        sigma_errs.append(sigma_err)
        chi2_vals.append(chi2)
        ndf_vals.append(float(ndf))
        chi2ndf_vals.append(chi2_ndf)

        csv_rows.append([
            low, high, int(entries),
            mean, mean_err,
            sigma, sigma_err,
            chi2, ndf, chi2_ndf
        ])

        h.SetStats(0)
        h.Draw()
        h.GetYaxis().SetTitle("Events")
        h.GetXaxis().SetTitle("q/p_{T} residual")

        if entries > 5:
            fit.Draw("same")
            leg = ROOT.TLegend(0.50, 0.46, 0.88, 0.84)
            leg.SetFillStyle(0)
            leg.SetBorderSize(0)
            leg.AddEntry(0, f"Entries = {int(entries)}", "")
            leg.AddEntry(0, f"Mean = {mean:.4f} #pm {mean_err:.4f}", "")
            leg.AddEntry(0, f"#sigma = {sigma:.4f} #pm {sigma_err:.4f}", "")
            leg.AddEntry(0, f"#chi^2 = {chi2:.4f}", "")
            leg.AddEntry(0, f"NDF = {int(ndf)}", "")
            leg.AddEntry(0, f"#chi^2/NDF = {chi2_ndf:.4f}", "")
            leg.Draw()

        out_file.cd()
        h.Write()
        fit.Write(f"fit_{bin_type}_{int(low)}_{int(high)}")
        c.SaveAs(os.path.join(base_dir, f"{bin_type}_Plots", f"{bin_type}_{muon_type}_{int(low)}_{int(high)}.png"))

        canvas_all.cd(i + 1)
        h.Draw()
        if entries > 5:
            fit.Draw("same")

    canvas_all.SaveAs(os.path.join(base_dir, f"{bin_type}_Plots", f"{muon_type}_{bin_type}_all_fits.png"))

    g_mean = ROOT.TGraphErrors(
        len(hlist),
        array_to_c(centers),
        array_to_c(means),
        array_to_c(halfwidths),
        array_to_c(mean_errs),
    )
    g_mean.SetName(f"{bin_type}_mean_total")
    g_mean.SetTitle(f"{bin_type} mean;{x_title};Mean of q/p_{{T}} relative residual")
    g_mean.SetMarkerStyle(8)
    g_mean.SetLineWidth(2)

    c_mean = ROOT.TCanvas(f"c_mean_{bin_type}", "", 1000, 750)
    c_mean.SetGrid()
    if logx:
        c_mean.SetLogx()
    c_mean.SetLeftMargin(0.18)
    g_mean.Draw("AP")
    out_file.cd()
    g_mean.Write()
    c_mean.SaveAs(os.path.join(base_dir, f"{bin_type}_Plots", f"{bin_type}_mean_total.png"))

    g_sigma = ROOT.TGraphErrors(
        len(hlist),
        array_to_c(centers),
        array_to_c(sigmas),
        array_to_c(halfwidths),
        array_to_c(sigma_errs),
    )
    g_sigma.SetName(f"{bin_type}_sigma_total")
    g_sigma.SetTitle(f"{bin_type} sigma;{x_title};#sigma of q/p_{{T}} relative residual")
    g_sigma.SetMarkerStyle(8)
    g_sigma.SetLineWidth(2)

    c_sigma = ROOT.TCanvas(f"c_sigma_{bin_type}", "", 1000, 750)
    c_sigma.SetGrid()
    if logx:
        c_sigma.SetLogx()
    c_sigma.SetLeftMargin(0.18)
    g_sigma.Draw("AP")
    out_file.cd()
    g_sigma.Write()
    c_sigma.SaveAs(os.path.join(base_dir, f"{bin_type}_Plots", f"{bin_type}_sigma_total.png"))

    g_chi2ndf = ROOT.TGraphErrors(
        len(hlist),
        array_to_c(centers),
        array_to_c(chi2ndf_vals),
        array_to_c(halfwidths),
        array_to_c([0.0] * len(hlist)),
    )
    g_chi2ndf.SetName(f"{bin_type}_chi2ndf_total")
    g_chi2ndf.SetTitle(f"{bin_type} #chi^{{2}}/NDF;{x_title};#chi^{{2}}/NDF")
    g_chi2ndf.SetMarkerStyle(8)
    g_chi2ndf.SetLineWidth(2)

    c_chi2 = ROOT.TCanvas(f"c_chi2ndf_{bin_type}", "", 1000, 750)
    c_chi2.SetGrid()
    if logx:
        c_chi2.SetLogx()
    c_chi2.SetLeftMargin(0.18)
    g_chi2ndf.Draw("AP")
    out_file.cd()
    g_chi2ndf.Write()
    c_chi2.SaveAs(os.path.join(base_dir, f"{bin_type}_Plots", f"{bin_type}_chi2ndf_total.png"))

    csv_path = os.path.join(base_dir, f"{bin_type}_Plots", f"{bin_type}_fit_summary.csv")
    write_fit_csv(csv_path, csv_rows)


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
    parser.add_argument("--no-trigger", action="store_true", help="Ignore HLT requirement even for DATA")
    parser.add_argument("--res-range", type=parse_bin_list, default=None, help='Residual plot range as "xmin,xmax"')
    parser.add_argument("--max-rel-pt-err-tag", type=float, default=None, help="Maximum allowed ptError/pt for tag muon")
    parser.add_argument("--max-rel-pt-err-probe", type=float, default=None, help="Maximum allowed ptError/pt for probe muon")

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

    if max_rel_pt_err_tag >= 0:
        print(f"Applying tag ptError/pt cut   < {max_rel_pt_err_tag}")
    if max_rel_pt_err_probe >= 0:
        print(f"Applying probe ptError/pt cut < {max_rel_pt_err_probe}")

    file_vec = ROOT.std.vector("string")()
    for f in files:
        file_vec.push_back(f)

    df = ROOT.RDataFrame(args.tree_name, file_vec)

    is_mc_cpp = "true" if args.datatype == "MC" else "false"
    use_trigger_cpp = "false" if args.no_trigger else "true"

    if args.muon_type == "DSA":
        pt_branch = "dmu_dsa_pt"
        pterr_branch = "dmu_dsa_ptError"
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
            ({use_trigger_cpp} ? HLT_L2Mu10_NoVertex_NoBPTX3BX : true),
            {pt_branch},
            {pterr_branch},
            {eta_branch},
            {phi_branch},
            {dz_branch},
            {dxy_branch},
            {charge_branch},
            {pass_tag_branch},
            {probe_id_branch},
            {has_probe_branch},
            {reco_type_branch},
            {is_mc_cpp},
            {args.min_pt},
            {max_rel_pt_err_tag},
            {max_rel_pt_err_probe}
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
    dz_histos = []
    dxy_histos = []

    for i in range(len(args.pt_bins) - 1):
        low, high = args.pt_bins[i], args.pt_bins[i + 1]
        cname = f"res_pt_bin_{i}"
        df_bins = df_bins.Define(cname, f"filter_by_bin(resolutions, tag_pt, {low}, {high})")
        h = df_bins.Histo1D(
            (f"pt_{int(low)}_{int(high)}", f"pT {low}-{high};q/p_{{T}} residual;Events", 100, res_min, res_max),
            cname
        )
        pt_histos.append((low, high, h))

    for i in range(len(args.dz_bins) - 1):
        low, high = args.dz_bins[i], args.dz_bins[i + 1]
        cname = f"res_dz_bin_{i}"
        df_bins = df_bins.Define(cname, f"filter_by_bin(resolutions, tag_dz_abs, {low}, {high})")
        h = df_bins.Histo1D(
            (f"dz_{int(low)}_{int(high)}", f"|dz| {low}-{high};q/p_{{T}} residual;Events", 100, res_min, res_max),
            cname
        )
        dz_histos.append((low, high, h))

    for i in range(len(args.dxy_bins) - 1):
        low, high = args.dxy_bins[i], args.dxy_bins[i + 1]
        cname = f"res_dxy_bin_{i}"
        df_bins = df_bins.Define(cname, f"filter_by_bin(resolutions, tag_dxy_abs, {low}, {high})")
        h = df_bins.Histo1D(
            (f"dxy_{int(low)}_{int(high)}", f"|dxy| {low}-{high};q/p_{{T}} residual;Events", 100, res_min, res_max),
            cname
        )
        dxy_histos.append((low, high, h))

    out_path = os.path.join(args.outdir, output_name)
    out = ROOT.TFile(out_path, "RECREATE")

    main_hists = [
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

    out.Close()
    print(f"Done. Output written to: {out_path}")


if __name__ == "__main__":
    main()


