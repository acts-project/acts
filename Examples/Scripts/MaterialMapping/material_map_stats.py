import ROOT


def make_bin_counts_dist(filename="material-maps.root"):

    f = ROOT.TFile.Open(filename, "READ")
    if not f or f.IsZombie():
        return

    for key in f.GetListOfKeys():

        obj = key.ReadObj()

        if not obj.InheritsFrom("TDirectory"):
            continue

        directory = obj

        h2 = directory.Get("binCounts")

        if not h2:
            print("Histogram not found - continue")
            continue

        print(f"\nProcessing: {directory.GetName()}")

        hdist = ROOT.TH1F(
            "binCountsDist",
            "BinCounts distribution;count;bins",
            20,
            0,
            h2.GetMaximum() + 1,
        )

        n_empty_bins = 0

        for ix in range(1, h2.GetNbinsX() + 1):
            for iy in range(1, h2.GetNbinsY() + 1):

                val = h2.GetBinContent(ix, iy)

                if val == 0:
                    n_empty_bins += 1

                hdist.Fill(val)

        print(
            f"Mean = {hdist.GetMean():.3f}, "
            f"StdDev = {hdist.GetStdDev():.3f}, "
            f"Empty bins counted = {n_empty_bins}"
        )

        c = ROOT.TCanvas(f"{directory.GetName()}binCounts", "", 800, 600)

        hdist.Draw()
        n_total_bins = h2.GetNbinsX() * h2.GetNbinsY()
        empty_fraction = n_empty_bins / n_total_bins
        pt = ROOT.TPaveText(0.60, 0.70, 0.90, 0.88, "NDC")
        pt.AddText(f"Empty bins: {n_empty_bins}")
        pt.AddText(f"Empty fraction: {100.0 * empty_fraction:.2f}%")
        pt.AddText(f"Mean: {hdist.GetMean():.2f}")
        pt.AddText(f"StdDev: {hdist.GetStdDev():.2f}")
        pt.Draw()

        c.SaveAs(f"{directory.GetName()}_binCounts.png")


if __name__ == "__main__":
    make_bin_counts_dist("material-maps.root")
