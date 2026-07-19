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
            return 0.0, 1.2
        if is_dz:
            return (0.0, 0.5) if is_dsa else (0.0, 0.3)
        if is_dxy:
            return (0.0, 0.3) if is_dsa else (0.0, 0.1)

    if "mean" in name:
        if is_pt:
            if is_dsa:
                return -0.5, 2.0
            if is_dgl:
                return -0.1, 0.1
            return -0.2, 0.2

        if is_dz:
            return (-0.2, 0.2) if is_dsa else (-0.05, 0.05)

        if is_dxy:
            return (-0.1, 0.1) if is_dsa else (-0.02, 0.02)

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


def draw_overlay(graphs, labels, title, ytitle, out_png, bin_type):
    c = ROOT.TCanvas(f"c_{os.path.basename(out_png)}", "", 900, 700)
    c.SetGrid()

    xmin, xmax, xtitle, logx = get_x_axis(bin_type)
    if logx:
        c.SetLogx()
    c.SetLeftMargin(0.16)

    ymin_auto, ymax_auto = get_graph_y_range(graphs, mode="robust")
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

    for g in graphs:
        g.Draw("P SAME")

    leg = ROOT.TLegend(0.58, 0.72, 0.88, 0.88)
    leg.SetBorderSize(0)
    leg.SetFillStyle(0)
    for g, label in zip(graphs, labels):
        leg.AddEntry(g, label, "pl")
    leg.Draw()

    c.SaveAs(out_png)


def load_scenario_graphs(files, scenario_names, source_name, bin_type, suffix):
    graphs = []
    for scenario in scenario_names:
        graph = files[scenario].Get(source_name)
        if not graph:
            raise RuntimeError(f"Could not find {source_name} in one or more files")
        graphs.append(graph.Clone(f"{bin_type}_{suffix}_{scenario}"))
    return graphs


def compare_family(files, outdir, bin_type, scenario_names, scenario_labels):
    mean_name = f"{bin_type}_mean_total"
    sigma_name = f"{bin_type}_sigma_total"
    hybrid_mean_name = f"{bin_type}_mean_total_hybrid"
    hybrid_sigma_name = f"{bin_type}_sigma_total_hybrid"

    mean_graphs = load_scenario_graphs(files, scenario_names, mean_name, bin_type, "mean")
    sigma_graphs = load_scenario_graphs(files, scenario_names, sigma_name, bin_type, "sigma")
    hybrid_mean_graphs = load_scenario_graphs(files, scenario_names, hybrid_mean_name, bin_type, "mean_hybrid")
    hybrid_sigma_graphs = load_scenario_graphs(files, scenario_names, hybrid_sigma_name, bin_type, "sigma_hybrid")

    colors = [ROOT.kBlack, ROOT.kRed + 1, ROOT.kBlue + 1]
    markers = [20, 21, 22]

    for i, graph in enumerate(mean_graphs):
        style_graph(graph, colors[i], markers[i])
    for i, graph in enumerate(sigma_graphs):
        style_graph(graph, colors[i], markers[i])
    for i, graph in enumerate(hybrid_mean_graphs):
        style_graph(graph, colors[i], markers[i])
    for i, graph in enumerate(hybrid_sigma_graphs):
        style_graph(graph, colors[i], markers[i])

    draw_overlay(
        mean_graphs,
        scenario_labels,
        f"{bin_type} mean scenario comparison",
        "Mean of q/p_{T} relative residual",
        os.path.join(outdir, f"{bin_type}_mean_scenarios.png"),
        bin_type
    )

    draw_overlay(
        sigma_graphs,
        scenario_labels,
        f"{bin_type} sigma scenario comparison",
        "#sigma of q/p_{T} relative residual",
        os.path.join(outdir, f"{bin_type}_sigma_scenarios.png"),
        bin_type
    )

    draw_overlay(
        hybrid_mean_graphs,
        scenario_labels,
        f"{bin_type} hybrid mean scenario comparison",
        "Mean of q/p_{T} relative residual",
        os.path.join(outdir, f"{bin_type}_mean_hybrid_scenarios.png"),
        bin_type
    )

    draw_overlay(
        hybrid_sigma_graphs,
        scenario_labels,
        f"{bin_type} hybrid sigma scenario comparison",
        "#sigma of q/p_{T} relative residual",
        os.path.join(outdir, f"{bin_type}_sigma_hybrid_scenarios.png"),
        bin_type
    )

    graphs = {}
    for scenario_name, graph in zip(scenario_names, mean_graphs):
        graphs[f"mean_{scenario_name}"] = graph
    for scenario_name, graph in zip(scenario_names, sigma_graphs):
        graphs[f"sigma_{scenario_name}"] = graph
    for scenario_name, graph in zip(scenario_names, hybrid_mean_graphs):
        graphs[f"mean_hybrid_{scenario_name}"] = graph
    for scenario_name, graph in zip(scenario_names, hybrid_sigma_graphs):
        graphs[f"sigma_hybrid_{scenario_name}"] = graph
    return graphs


def resolve_scenarios(args):
    if args.inputs or args.labels or args.names:
        if not (args.inputs and args.labels and args.names):
            raise ValueError("--inputs, --labels, and --names must be provided together")
        if not (2 <= len(args.inputs) <= 3):
            raise ValueError("--inputs must contain two or three ROOT files")
        if len(args.inputs) != len(args.labels) or len(args.inputs) != len(args.names):
            raise ValueError("--inputs, --labels, and --names must have the same length")
        return args.names, args.labels, args.inputs

    legacy_args = [args.baseline, args.tagonly, args.tagprobe]
    if all(legacy_args):
        return (
            ["baseline", "tagonly", "tagprobe"],
            [
                "Baseline",
                "Tag p_{T}^{err}/p_{T}<0.2",
                "Tag & Probe p_{T}^{err}/p_{T}<0.2",
            ],
            legacy_args,
        )

    raise ValueError("Provide either legacy --baseline/--tagonly/--tagprobe inputs or the new --inputs/--labels/--names triplet")


def main():
    parser = argparse.ArgumentParser(description="Compare baseline, tag-only, and tag+probe for pt, dz, dxy")
    parser.add_argument("--baseline", help="Baseline ROOT file")
    parser.add_argument("--tagonly", help="Tag-only cut ROOT file")
    parser.add_argument("--tagprobe", help="Tag+probe cut ROOT file")
    parser.add_argument("--inputs", nargs="+", help="Scenario ROOT files to compare")
    parser.add_argument("--labels", nargs="+", help="Legend labels matching --inputs")
    parser.add_argument("--names", nargs="+", help="Short scenario names matching --inputs, used in ROOT object names")
    parser.add_argument("--outdir", default="scenario_comparison_plots", help="Output directory")
    args = parser.parse_args()

    os.makedirs(args.outdir, exist_ok=True)

    scenario_names, scenario_labels, input_paths = resolve_scenarios(args)
    files = {name: ROOT.TFile.Open(path) for name, path in zip(scenario_names, input_paths)}

    for name, f in files.items():
        if not f or f.IsZombie():
            raise RuntimeError(f"Could not open {name} file")

    out_root = ROOT.TFile(os.path.join(args.outdir, "scenario_comparison.root"), "RECREATE")

    for bin_type in ["pt", "dz", "dxy"]:
        graphs = compare_family(files, args.outdir, bin_type, scenario_names, scenario_labels)
        out_root.cd()
        for key, graph in graphs.items():
            graph.Write(f"{bin_type}_{key}")

    out_root.Close()

    for f in files.values():
        f.Close()

    print(f"Saved scenario comparison plots in: {args.outdir}")


if __name__ == "__main__":
    main()
