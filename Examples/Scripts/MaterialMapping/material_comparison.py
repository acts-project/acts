import ROOT
import argparse
import math

from ROOT import (
    TCanvas,
    TFile,
    TTree,
    gDirectory,
    gStyle,
    TH1D,
    TH2F,
    TProfile,
    TRatioPlot,
)


def TH1D_from_TProf(tprof):
    h = TH1D(
        tprof.GetName() + "_th1",
        tprof.GetTitle(),
        tprof.GetNbinsX(),
        tprof.GetXaxis().GetXmin(),
        tprof.GetXaxis().GetXmax(),
    )
    for i in range(tprof.GetNbinsX()):
        if tprof.GetBinContent(i + 1) == 0.0:
            h.SetBinContent(i + 1, 0.0)
            h.SetBinError(i + 1, 1000.0)
            continue
        h.SetBinContent(i + 1, tprof.GetBinContent(i + 1))
        h.SetBinError(i + 1, tprof.GetBinError(i + 1))
    return h


if "__main__" == __name__:
    p = argparse.ArgumentParser()

    p.add_argument(
        "-e", "--entries", type=int, default=100000, help="Number of events to process"
    )
    p.add_argument(
        "-i",
        "--input",
        type=str,
        nargs="+",
        default="",
        help="Input files with material tracks",
    )
    p.add_argument(
        "-o",
        "--output",
        type=str,
        default="",
        help="Output file with produced material overview plots",
    )
    p.add_argument(
        "-l",
        "--labels",
        type=str,
        nargs="+",
        default="",
        help="The labels for the input files",
    )
    p.add_argument(
        "-v",
        "--variables",
        type=str,
        nargs="+",
        default=["t_X0", "t_L0"],
        help="The variables to be plotted",
    )
    p.add_argument(
        "-m",
        "--max",
        type=float,
        nargs="+",
        default=[4.5, 1.5],
        help="The variables to be plotted",
    )
    p.add_argument(
        "-a",
        "--axes",
        type=str,
        nargs="+",
        default=["v_eta", "v_phi"],
        help="The axes versus which the variables to be plotted",
    )
    p.add_argument(
        "--eta",
        type=float,
        nargs=2,
        default=[-3.0, 3.0],
        help="Eta range for the plotting",
    )
    p.add_argument("--eta-bins", type=int, default=60, help="Eta bins for the plotting")
    p.add_argument(
        "--phi",
        type=float,
        nargs=2,
        default=[-math.pi, math.pi],
        help="Phi range for the plotting",
    )
    p.add_argument("--phi-bins", type=int, default=72, help="Phi bins for the plotting")
    p.add_argument(
        "-d",
        "--detray-input",
        type=str,
        default="",
        help="Optional detray input csv file",
    )

    args = p.parse_args()

    ttrees = []

    if len(args.input) != len(args.labels):
        print("** ERROR ** The number of input files and labels must match")
        exit(1)

    # The histogram map
    h_dict = {}
    tfiles = []

    ofile = ROOT.TFile.Open(args.output, "RECREATE") if args.output else None

    input_files = args.input

    # Detray csv file
    if args.detray_input:
        input_files.append(args.detray_input)

    # Loop over the files and create the comparison histograms
    for fi, ifile in enumerate(input_files):
        # Special treament for detray input
        if ifile == args.detray_input:
            print("Reading detray input from file: ", args.detray_input)
            ttree = TTree("t", "material read from detray")
            ttree.ReadFile(args.detray_input, "v_eta/D:v_phi:t_X0:t_L0:t_X0:t_L0")
            lbl = "detray"
        else:
            # Get the file
            tfile = ROOT.TFile.Open(ifile)
            tfiles.append(tfile)
            ttree = tfile.Get("material-tracks")
            ttrees.append(ttree)
            # The label
            lbl = args.labels[fi]

        # Loop over the variables and axes
        for iv, var in enumerate(args.variables):
            for ax in args.axes:
                if ax == "v_eta":
                    bins = args.eta_bins
                    low = args.eta[0]
                    high = args.eta[1]
                    cut = f"v_phi > {args.phi[0]} && v_phi < {args.phi[1]}"
                elif ax == "v_phi":
                    bins = args.phi_bins
                    low = args.phi[0]
                    high = args.phi[1]
                    cut = f"v_eta > {args.eta[0]} && v_eta < {args.eta[1]}"
                # Some naming magic
                hcmmd = f"{var}:{ax}>>"
                hname = f"{var}_vs_{ax}"
                hfname = f"{lbl}_{var}_vs_{ax}"
                hrange = f"({bins},{low},{high},100,0,{args.max[iv]})"
                htitle = f"{var} vs {ax}"
                ttree.Draw(hcmmd + hfname + hrange, cut, "")
                h = ROOT.gDirectory.Get(hfname)
                # Fill into comparison histogram
                if h_dict.get(hname):
                    h_dict[hname].append(h)
                else:
                    h_dict[hname] = [h]

                # Write to file
                if ofile:
                    ofile.cd()
                    h.Write()

    colors = [
        ROOT.kBlack,
        ROOT.kRed + 2,
        ROOT.kBlue - 4,
        ROOT.kGreen + 1,
        ROOT.kYellow - 2,
    ]
    markers = [
        ROOT.kFullCircle,
        ROOT.kFullCircle,
        ROOT.kFullTriangleUp,
        ROOT.kFullTriangleDown,
        ROOT.kFullDiamond,
    ]

    # Now create the comparison histograms
    c = ROOT.TCanvas("Comparison", "Comparison", 1200, 800)
    c.Divide(2, 2)

    # Remove the stat box
    gStyle.SetOptStat(0)

    # Memory garbage collection, thanks ROOT
    hist_memory_pool = []

    ic = 0
    for hname in h_dict:
        ic += 1
        c.cd(ic)
        h_list = h_dict[hname]
        h_ref = None
        h_prof_ref = None
        for ih, h in enumerate(h_list):
            h_prof = TH1D_from_TProf(h.ProfileX())
            hist_memory_pool.append(h_prof)
            h_prof.SetObjectStat(0)
            h_prof.SetLineColor(colors[ih])
            h_prof.SetMarkerColor(colors[ih])
            h_prof.SetMarkerSize(0.5)
            h_prof.SetMarkerStyle(markers[ih])
            h_prof.GetYaxis().SetRangeUser(0.0, 1.3 * h_prof.GetMaximum())
            if ih == 0:
                h_ref = h
                h_prof_ref = h_prof
            else:
                h_ratio = ROOT.TRatioPlot(h_prof_ref, h_prof)
                h_ratio.SetGraphDrawOpt("pe")
                h_ratio.SetSeparationMargin(0.005)
                drawOption = "e,same" if ih > 1 else "e"
                h_ratio.Draw(drawOption)
                h_ratio.GetLowerRefGraph().SetLineColor(colors[ih])
                h_ratio.GetLowerRefGraph().SetMarkerColor(colors[ih])
                h_ratio.GetLowerRefGraph().SetMarkerSize(0.5)
                h_ratio.GetLowerRefGraph().SetMarkerStyle(markers[ih])
                c.Update()
                hist_memory_pool.append(h_ratio)
        h_prof.Draw(("same" if ih > 0 else ""))
        c.SaveAs(f"{hname}.png")
