import ROOT

# Build a TCanvas
canvas = ROOT.TCanvas("canvas", "Scatter Plot Comparison", 800, 600)

# Build WO-Time TGraph
time_stddev = [0, 2, 4, 6, 8, 10]

te_wotime = [0.944756, 0.944756, 0.944756, 0.944756, 0.944756, 0.944756]
n_points = len(time_stddev)
graph1 = ROOT.TGraph(n_points)

for i, (x, y) in enumerate(zip(time_stddev, te_wotime)):
    graph1.SetPoint(i, x, y)

graph1.SetMarkerStyle(ROOT.kFullCircle)
graph1.SetMarkerSize(1.0)
graph1.SetMarkerColor(ROOT.kBlue)
graph1.GetYaxis().SetRangeUser(0.7, 1.2)
graph1.SetLineStyle(1)
graph1.SetLineColor(ROOT.kBlue)
graph1.SetTitle("Tracking Efficiency Compare about Time Distribution")
graph1.Draw("APL")  # "A" as Axis，"P" as Plot

# Build W-Time TGraph
te_wtime = [0.943171, 1, 1, 1, 1, 1]
graph2 = ROOT.TGraph(n_points)

for i, (x, y) in enumerate(zip(time_stddev, te_wtime)):
    graph2.SetPoint(i, x, y)

graph2.SetMarkerStyle(ROOT.kFullTriangleUp)
graph2.SetMarkerSize(1.0)
graph2.SetMarkerColor(ROOT.kRed)
graph2.GetYaxis().SetRangeUser(0.7, 1.2)
graph2.SetLineStyle(1)
graph2.SetLineColor(ROOT.kRed)
graph2.Draw("PL same")  # "same"

# Axis Title
graph1.GetXaxis().SetTitle("StdDev_T")
graph1.GetYaxis().SetTitle("Tracking Efficiency")

# Legend
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(graph1, "WO-Time", "p")
legend.AddEntry(graph2, "W-Time", "p")
legend.Draw()

# Title
canvas.SetTitle("Tracking Efficiency Compare about Time Distribution")

canvas.SaveAs("TE-stddev-time.pdf")

