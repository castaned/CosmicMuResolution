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

    is_pt = "pt_" in name
    is_dz = "dz_" in name
    is_dxy = "dxy_" in name

    if "sigma" in name:
        if is_pt:
            return (0.0, 3.6) if is_dsa else (0.0, 0.15)
        if is_dz:
            return (0.0, 0.8) if is_dsa else (0.0, 0.08)
        if is_dxy:
            return (0.0, 0.6) if is_dsa else (0.0, 0.08)

    if "mean" in name:
        if is_pt:
            if is_dsa:
                return -1.0, 0.6
            if is_dgl:
                return -0.01, 0.02
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


def draw_overlay(g_mc, g_data, title, ytitle, out_png, bin_type):
    c = ROOT.TCanvas(f"c_{os.path.basename(out_png)}", "", 900, 700)
    c.SetGrid()

    xmin, xmax, xtitle, logx = get_x_axis(bin_type)
    if logx:
        c.SetLogx()
    c.SetLeftMargin(0.16)

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

    c.SaveAs(out_png)


def compare_family(f_mc, f_data, outdir, bin_type):
    mean_name = f"{bin_type}_mean_total"
    sigma_name = f"{bin_type}_sigma_total"

    g_mean_mc = f_mc.Get(mean_name)
    g_mean_data = f_data.Get(mean_name)
    g_sigma_mc = f_mc.Get(sigma_name)
    g_sigma_data = f_data.Get(sigma_name)

    if not g_mean_mc or not g_mean_data:
        raise RuntimeError(f"Could not find {mean_name} in one or both files")
    if not g_sigma_mc or not g_sigma_data:
        raise RuntimeError(f"Could not find {sigma_name} in one or both files")

    g_mean_mc = g_mean_mc.Clone(f"{bin_type}_mean_mc")
    g_mean_data = g_mean_data.Clone(f"{bin_type}_mean_data")
    g_sigma_mc = g_sigma_mc.Clone(f"{bin_type}_sigma_mc")
    g_sigma_data = g_sigma_data.Clone(f"{bin_type}_sigma_data")

    style_graph(g_mean_mc, ROOT.kRed + 1, 20)
    style_graph(g_mean_data, ROOT.kBlue + 1, 21)
    style_graph(g_sigma_mc, ROOT.kRed + 1, 20)
    style_graph(g_sigma_data, ROOT.kBlue + 1, 21)

    draw_overlay(
        g_mean_mc,
        g_mean_data,
        f"{bin_type} mean comparison",
        "Mean of q/p_{T} relative residual",
        os.path.join(outdir, f"{bin_type}_mean_mc_vs_data.png"),
        bin_type
    )

    draw_overlay(
        g_sigma_mc,
        g_sigma_data,
        f"{bin_type} sigma comparison",
        "#sigma of q/p_{T} relative residual",
        os.path.join(outdir, f"{bin_type}_sigma_mc_vs_data.png"),
        bin_type
    )

    return g_mean_mc, g_mean_data, g_sigma_mc, g_sigma_data


def main():
    parser = argparse.ArgumentParser(description="Compare mean/sigma between MC and DATA for pt, dz, dxy")
    parser.add_argument("--mc", required=True, help="MC ROOT file")
    parser.add_argument("--data", required=True, help="DATA ROOT file")
    parser.add_argument("--outdir", default="comparison_plots", help="Output directory")
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
        g_mean_mc, g_mean_data, g_sigma_mc, g_sigma_data = compare_family(f_mc, f_data, args.outdir, bin_type)
        out_root.cd()
        g_mean_mc.Write(f"{bin_type}_mean_mc")
        g_mean_data.Write(f"{bin_type}_mean_data")
        g_sigma_mc.Write(f"{bin_type}_sigma_mc")
        g_sigma_data.Write(f"{bin_type}_sigma_data")

    out_root.Close()
    f_mc.Close()
    f_data.Close()

    print(f"Saved comparison plots in: {args.outdir}")


if __name__ == "__main__":
    main()
