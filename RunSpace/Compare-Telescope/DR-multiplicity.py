import ROOT

# Build a TCanvas
canvas = ROOT.TCanvas("canvas", "Scatter Plot Comparison", 800, 600)

# Build WO-Time TGraph
multiplicity = [1, 2, 3, 4, 5]

dr_wotime = [0, 0.895397, 0.928412, 0.890093, 0.811839, ]
n_points = len(multiplicity)
graph1 = ROOT.TGraph(n_points)

for i, (x, y) in enumerate(zip(multiplicity, dr_wotime)):
    graph1.SetPoint(i, x, y)

graph1.SetMarkerStyle(ROOT.kFullCircle)
graph1.SetMarkerSize(1.0)
graph1.SetMarkerColor(ROOT.kBlue)
graph1.GetYaxis().SetRangeUser(0, 1.1)
graph1.SetLineStyle(1)
graph1.SetLineColor(ROOT.kBlue)
graph1.SetTitle("Duplicate Rate Compare about Multiplicity")
graph1.Draw("APL")  # "A" as Axis，"P" as Plot

# Build W-Time TGraph
dr_wtime = [0, 0.884925, 0.925721, 0.897075, 0.810726]
graph2 = ROOT.TGraph(n_points)

for j, (x, y) in enumerate(zip(multiplicity, dr_wtime)):
    graph2.SetPoint(j, x, y)

graph2.SetMarkerStyle(ROOT.kFullTriangleUp)
graph2.SetMarkerSize(1.0)
graph2.SetMarkerColor(ROOT.kRed)
graph2.GetYaxis().SetRangeUser(0, 1.1)
graph2.SetLineStyle(1)
graph2.SetLineColor(ROOT.kRed)
graph2.Draw("PL same")  # "same"

# Axis Title
graph1.GetXaxis().SetTitle("Multiplicity")
graph1.GetYaxis().SetTitle("Duplicate Rate")

# Legend
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(graph1, "WO-Time", "p")
legend.AddEntry(graph2, "W-Time", "p")
legend.Draw()

# Title
canvas.SetTitle("Duplicate Rate Compare about Multiiplicity")

canvas.SaveAs("DR-Multiplicity.pdf")

