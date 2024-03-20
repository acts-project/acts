import ROOT

# Build a TCanvas
canvas = ROOT.TCanvas("canvas", "Scatter Plot Comparison", 800, 600)

# Build WO-Time TGraph
pos_stddev = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200]

dr_wotime = [0.928412, 0.925492, 0.919233, 0.908279, 0.893536, 0.877115, 0.869841, 0.849324, 0.831738, 0.810732, 0.757979]
n_points = len(pos_stddev)
graph1 = ROOT.TGraph(n_points)

for i, (x, y) in enumerate(zip(pos_stddev, dr_wotime)):
    graph1.SetPoint(i, x, y)

graph1.SetMarkerStyle(ROOT.kFullCircle)
graph1.SetMarkerSize(1.0)
graph1.SetMarkerColor(ROOT.kBlue)
graph1.GetYaxis().SetRangeUser(0, 1.2)
graph1.SetLineStyle(1)
graph1.SetLineColor(ROOT.kBlue)
graph1.SetTitle("Duplicate Rate Compare about Position Distribution")
graph1.Draw("APL")  # "A" as Axis，"P" as Plot

# Build W-Time TGraph
dr_wtime = [0.340659, 0.34434, 0.125, 0.115702, 0.0892857, 0.0673077, 0.065, 0.0276243, 0.0285714, 0.0120482, 0.00641026]
graph2 = ROOT.TGraph(n_points)

for i, (x, y) in enumerate(zip(pos_stddev, dr_wtime)):
    graph2.SetPoint(i, x, y)

graph2.SetMarkerStyle(ROOT.kFullTriangleUp)
graph2.SetMarkerSize(1.0)
graph2.SetMarkerColor(ROOT.kRed)
graph2.GetYaxis().SetRangeUser(0, 1.2)
graph2.SetLineStyle(1)
graph2.SetLineColor(ROOT.kRed)
graph2.Draw("PL same")  # "same"

# Axis Title
graph1.GetXaxis().SetTitle("StdDev_Pos")
graph1.GetYaxis().SetTitle("Duplicate Rate")

# Legend
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(graph1, "WO-Time", "p")
legend.AddEntry(graph2, "W-Time", "p")
legend.Draw()

# Title
canvas.SetTitle("Duplicate Rate Compare about Position Distribution")

canvas.SaveAs("DR-stddev-pos.pdf")

