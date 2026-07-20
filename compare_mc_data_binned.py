#!/usr/bin/env python3

import ROOT
import os
import argparse

ROOT.gStyle.SetOptStat(0)


def style_graph(gr, color, marker):
    gr.SetLineColor(color)
    gr.SetMarkerColor(color)
    gr.SetMarkerStyle(marker)
    gr.SetLineWidth(2)


def get_graph_y_range(graphs, mode="robust"):
    vals = []

    for g in graphs:
        n = g.GetN()
        y = g.GetY()
        ey = g.GetEY()

        for i in range(n):
            yi = y[i]
            eyi = ey[i]

            if mode == "robust":
                if yi > 0 and eyi > 2.0 * abs(yi):
                    continue
                if eyi > 5.0:
                    continue

            low = yi - eyi
            high = yi + eyi
            vals.append((low, high))

    if not vals:
        return 0.0, 1.0

    ymin = min(v[0] for v in vals)
    ymax = max(v[1] for v in vals)

    if ymin == ymax:
        if ymin == 0:
            return -1.0, 1.0
        return ymin * 0.8, ymax * 1.2

    pad = 0.15 * (ymax - ymin)
    return ymin - pad, ymax + pad


def choose_fixed_range(out_png, ymin_auto, ymax_auto):
    name = out_png.lower()

    is_dsa = "dsa" in name
    is_dgl = "dgl" in name
    is_sym = "_sym_" in name or "sym_mc_vs_data" in name
    is_hybrid = "hybrid" in name
    is_double = "double" in name

    is_pt = "pt_" in name
    is_dz = "dz_" in name
    is_dxy = "dxy_" in name

    is_standard_dsa_pt = is_dsa and is_pt and not (is_sym or is_hybrid or is_double)

    if is_standard_dsa_pt:
        if "sigma" in name:
            return 0.0, 2.0
        if "mean" in name:
            return -0.35, 0.40

    if is_dsa and is_sym:
        if "sigma" in name and is_pt:
            return 0.2, 0.8
        if "mean" in name and is_pt:
            return -0.4, 0.65

    if "sigma" in name:
        if is_pt:
            return (0.0, 3.6) if is_dsa else (0.0, 0.12)
        if is_dz:
            return (0.0, 0.8) if is_dsa else (0.0, 0.08)
        if is_dxy:
            return (0.0, 0.6) if is_dsa else (0.0, 0.08)

    if "mean" in name:
        if is_pt:
            if is_dsa:
                return -1.0, 0.6
            if is_dgl:
                return -0.008, 0.004
            return -0.2, 0.2

        if is_dz:
            return (-0.4, 0.05) if is_dsa else (-0.025, 0.01)

        if is_dxy:
            return (-0.4, 0.05) if is_dsa else (-0.01, 0.005)

    return ymin_auto, ymax_auto


def get_x_axis(bin_type):
    if bin_type == "pt":
        return 15.0, 1200.0, "p_{T}^{tag} [GeV]", True
    if bin_type == "dz":
        return 0.8, 200.0, "|dz|", True
    if bin_type == "dxy":
        return 0.8, 100.0, "|dxy|", True
    return 0.0, 1.0, bin_type, False


def make_frame(xmin, xmax, ymin, ymax, xtitle, ytitle, name):
    h = ROOT.TH1F(name, "", 100, xmin, xmax)
    h.SetMinimum(ymin)
    h.SetMaximum(ymax)
    h.GetXaxis().SetTitle(xtitle)
    h.GetYaxis().SetTitle(ytitle)
    h.SetStats(0)
    return h


def draw_header(canvas, left_label, right_label):
    latex = ROOT.TLatex()
    latex.SetNDC()
    left_margin = canvas.GetLeftMargin()
    top_y = 1.0 - canvas.GetTopMargin() + 0.034

    cms_label = "CMS"
    extra_label = left_label.replace("CMS", "", 1).strip() if left_label.startswith("CMS") else left_label

    latex.SetTextAlign(13)
    latex.SetTextFont(61)
    latex.SetTextSize(0.050)
    latex.DrawLatex(left_margin, top_y, cms_label)

    if extra_label:
        latex.SetTextFont(52)
        latex.SetTextSize(0.038)
        latex.DrawLatex(left_margin + 0.105, top_y, extra_label)

    if right_label:
        latex.SetTextAlign(33)
        latex.SetTextFont(42)
        latex.SetTextSize(0.040)
        latex.DrawLatex(0.88, top_y, right_label)


def draw_overlay(g_mc, g_data, title, ytitle, out_png, bin_type, left_label, right_label):
    c = ROOT.TCanvas(f"c_{os.path.basename(out_png)}", "", 900, 700)
    c.SetGrid()

    xmin, xmax, xtitle, logx = get_x_axis(bin_type)
    if logx:
        c.SetLogx()
    c.SetLeftMargin(0.16)
    c.SetTopMargin(0.09)

    ymin_auto, ymax_auto = get_graph_y_range([g_mc, g_data], mode="robust")
    ymin, ymax = choose_fixed_range(out_png, ymin_auto, ymax_auto)

    frame = make_frame(
        xmin,
        xmax,
        ymin,
        ymax,
        xtitle,
        ytitle,
        f"frame_{os.path.basename(out_png).replace('.', '_')}"
    )
    frame.Draw()

    g_mc.Draw("P SAME")
    g_data.Draw("P SAME")

    leg = ROOT.TLegend(0.65, 0.75, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    leg.AddEntry(g_mc, "MC", "pl")
    leg.AddEntry(g_data, "Data", "pl")
    leg.Draw()
    draw_header(c, left_label, right_label)

    c.SaveAs(out_png)


def load_comparison_graphs(f_mc, f_data, source_name, mc_name, data_name):
    g_mc = f_mc.Get(source_name)
    g_data = f_data.Get(source_name)

    if not g_mc or not g_data:
        raise RuntimeError(f"Could not find {source_name} in one or both files")

    return g_mc.Clone(mc_name), g_data.Clone(data_name)


def compare_family(f_mc, f_data, outdir, bin_type, left_label, right_label, source_suffix="", output_suffix="", ytitle_prefix=""):
    mean_name = f"{bin_type}_mean_total{source_suffix}"
    sigma_name = f"{bin_type}_sigma_total{source_suffix}"
    hybrid_mean_name = f"{bin_type}_mean_total{source_suffix}_hybrid"
    hybrid_sigma_name = f"{bin_type}_sigma_total{source_suffix}_hybrid"
    double_mean_name = f"{bin_type}_mean_total{source_suffix}_double"
    sigma_eff_name = f"{bin_type}_sigma_eff_total{source_suffix}_double"

    g_mean_mc, g_mean_data = load_comparison_graphs(
        f_mc, f_data, mean_name, f"{bin_type}_mean_mc", f"{bin_type}_mean_data"
    )
    g_sigma_mc, g_sigma_data = load_comparison_graphs(
        f_mc, f_data, sigma_name, f"{bin_type}_sigma_mc", f"{bin_type}_sigma_data"
    )
    g_hybrid_mean_mc, g_hybrid_mean_data = load_comparison_graphs(
        f_mc, f_data, hybrid_mean_name, f"{bin_type}_mean_hybrid_mc", f"{bin_type}_mean_hybrid_data"
    )
    g_hybrid_sigma_mc, g_hybrid_sigma_data = load_comparison_graphs(
        f_mc, f_data, hybrid_sigma_name, f"{bin_type}_sigma_hybrid_mc", f"{bin_type}_sigma_hybrid_data"
    )
    g_double_mean_mc, g_double_mean_data = load_comparison_graphs(
        f_mc, f_data, double_mean_name, f"{bin_type}_mean_double_mc", f"{bin_type}_mean_double_data"
    )
    g_sigma_eff_mc, g_sigma_eff_data = load_comparison_graphs(
        f_mc, f_data, sigma_eff_name, f"{bin_type}_sigma_eff_double_mc", f"{bin_type}_sigma_eff_double_data"
    )

    style_graph(g_mean_mc, ROOT.kRed + 1, 20)
    style_graph(g_mean_data, ROOT.kBlue + 1, 21)
    style_graph(g_sigma_mc, ROOT.kRed + 1, 20)
    style_graph(g_sigma_data, ROOT.kBlue + 1, 21)
    style_graph(g_hybrid_mean_mc, ROOT.kRed + 1, 20)
    style_graph(g_hybrid_mean_data, ROOT.kBlue + 1, 21)
    style_graph(g_hybrid_sigma_mc, ROOT.kRed + 1, 20)
    style_graph(g_hybrid_sigma_data, ROOT.kBlue + 1, 21)
    style_graph(g_double_mean_mc, ROOT.kRed + 1, 20)
    style_graph(g_double_mean_data, ROOT.kBlue + 1, 21)
    style_graph(g_sigma_eff_mc, ROOT.kRed + 1, 20)
    style_graph(g_sigma_eff_data, ROOT.kBlue + 1, 21)

    draw_overlay(
        g_mean_mc,
        g_mean_data,
        f"{bin_type} mean comparison",
        f"Mean of {ytitle_prefix}q/p_{{T}} relative residual",
        os.path.join(outdir, f"{bin_type}_mean{output_suffix}_mc_vs_data.png"),
        bin_type,
        left_label,
        right_label
    )

    draw_overlay(
        g_sigma_mc,
        g_sigma_data,
        f"{bin_type} sigma comparison",
        f"#sigma of {ytitle_prefix}q/p_{{T}} relative residual",
        os.path.join(outdir, f"{bin_type}_sigma{output_suffix}_mc_vs_data.png"),
        bin_type,
        left_label,
        right_label
    )

    draw_overlay(
        g_hybrid_mean_mc,
        g_hybrid_mean_data,
        f"{bin_type} hybrid mean comparison",
        f"Mean of {ytitle_prefix}q/p_{{T}} relative residual",
        os.path.join(outdir, f"{bin_type}_mean{output_suffix}_hybrid_mc_vs_data.png"),
        bin_type,
        left_label,
        right_label
    )

    draw_overlay(
        g_hybrid_sigma_mc,
        g_hybrid_sigma_data,
        f"{bin_type} hybrid sigma comparison",
        f"#sigma of {ytitle_prefix}q/p_{{T}} relative residual",
        os.path.join(outdir, f"{bin_type}_sigma{output_suffix}_hybrid_mc_vs_data.png"),
        bin_type,
        left_label,
        right_label
    )

    draw_overlay(
        g_double_mean_mc,
        g_double_mean_data,
        f"{bin_type} double-gaussian mean comparison",
        f"Mean of {ytitle_prefix}q/p_{{T}} relative residual",
        os.path.join(outdir, f"{bin_type}_mean{output_suffix}_double_mc_vs_data.png"),
        bin_type,
        left_label,
        right_label
    )

    draw_overlay(
        g_sigma_eff_mc,
        g_sigma_eff_data,
        f"{bin_type} double-gaussian effective sigma comparison",
        f"#sigma_{{eff}} of {ytitle_prefix}q/p_{{T}} relative residual",
        os.path.join(outdir, f"{bin_type}_sigma_eff{output_suffix}_double_mc_vs_data.png"),
        bin_type,
        left_label,
        right_label
    )

    return {
        "mean_mc": g_mean_mc,
        "mean_data": g_mean_data,
        "sigma_mc": g_sigma_mc,
        "sigma_data": g_sigma_data,
        "mean_hybrid_mc": g_hybrid_mean_mc,
        "mean_hybrid_data": g_hybrid_mean_data,
        "sigma_hybrid_mc": g_hybrid_sigma_mc,
        "sigma_hybrid_data": g_hybrid_sigma_data,
        "mean_double_mc": g_double_mean_mc,
        "mean_double_data": g_double_mean_data,
        "sigma_eff_double_mc": g_sigma_eff_mc,
        "sigma_eff_double_data": g_sigma_eff_data,
    }


def has_graph_pair(f_mc, f_data, source_name):
    return bool(f_mc.Get(source_name) and f_data.Get(source_name))


def main():
    parser = argparse.ArgumentParser(description="Compare mean/sigma between MC and DATA for pt, dz, dxy")
    parser.add_argument("--mc", required=True, help="MC ROOT file")
    parser.add_argument("--data", required=True, help="DATA ROOT file")
    parser.add_argument("--outdir", default="comparison_plots", help="Output directory")
    parser.add_argument("--period", default="", help="Period label shown at the top right, e.g. 2023D")
    parser.add_argument("--left-label", default="CMS Cosmics Preliminary", help="Header label shown at the top left")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    f_mc = ROOT.TFile.Open(args.mc)
    f_data = ROOT.TFile.Open(args.data)

    if not f_mc or f_mc.IsZombie():
        raise RuntimeError(f"Could not open MC file: {args.mc}")
    if not f_data or f_data.IsZombie():
        raise RuntimeError(f"Could not open DATA file: {args.data}")

    out_root = ROOT.TFile(os.path.join(args.outdir, "mc_data_comparison.root"), "RECREATE")

    for bin_type in ["pt", "dz", "dxy"]:
        graphs = compare_family(f_mc, f_data, args.outdir, bin_type, args.left_label, args.period)
        out_root.cd()
        graphs["mean_mc"].Write(f"{bin_type}_mean_mc")
        graphs["mean_data"].Write(f"{bin_type}_mean_data")
        graphs["sigma_mc"].Write(f"{bin_type}_sigma_mc")
        graphs["sigma_data"].Write(f"{bin_type}_sigma_data")
        graphs["mean_hybrid_mc"].Write(f"{bin_type}_mean_hybrid_mc")
        graphs["mean_hybrid_data"].Write(f"{bin_type}_mean_hybrid_data")
        graphs["sigma_hybrid_mc"].Write(f"{bin_type}_sigma_hybrid_mc")
        graphs["sigma_hybrid_data"].Write(f"{bin_type}_sigma_hybrid_data")
        graphs["mean_double_mc"].Write(f"{bin_type}_mean_double_mc")
        graphs["mean_double_data"].Write(f"{bin_type}_mean_double_data")
        graphs["sigma_eff_double_mc"].Write(f"{bin_type}_sigma_eff_double_mc")
        graphs["sigma_eff_double_data"].Write(f"{bin_type}_sigma_eff_double_data")
        if has_graph_pair(f_mc, f_data, f"{bin_type}_mean_total_sym"):
            sym_graphs = compare_family(
                f_mc,
                f_data,
                args.outdir,
                bin_type,
                args.left_label,
                args.period,
                source_suffix="_sym",
                output_suffix="_sym",
                ytitle_prefix="symmetric ",
            )
            out_root.cd()
            sym_graphs["mean_mc"].Write(f"{bin_type}_mean_mc_sym")
            sym_graphs["mean_data"].Write(f"{bin_type}_mean_data_sym")
            sym_graphs["sigma_mc"].Write(f"{bin_type}_sigma_mc_sym")
            sym_graphs["sigma_data"].Write(f"{bin_type}_sigma_data_sym")
            sym_graphs["mean_hybrid_mc"].Write(f"{bin_type}_mean_hybrid_mc_sym")
            sym_graphs["mean_hybrid_data"].Write(f"{bin_type}_mean_hybrid_data_sym")
            sym_graphs["sigma_hybrid_mc"].Write(f"{bin_type}_sigma_hybrid_mc_sym")
            sym_graphs["sigma_hybrid_data"].Write(f"{bin_type}_sigma_hybrid_data_sym")
            sym_graphs["mean_double_mc"].Write(f"{bin_type}_mean_double_mc_sym")
            sym_graphs["mean_double_data"].Write(f"{bin_type}_mean_double_data_sym")
            sym_graphs["sigma_eff_double_mc"].Write(f"{bin_type}_sigma_eff_double_mc_sym")
            sym_graphs["sigma_eff_double_data"].Write(f"{bin_type}_sigma_eff_double_data_sym")

    out_root.Close()
    f_mc.Close()
    f_data.Close()

    print(f"Saved comparison plots in: {args.outdir}")


if __name__ == "__main__":
    main()
