import ROOT

# 文件路径列表
file_paths_wot = [
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wot_multiplicity_1/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wot_multiplicity_2/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wot_multiplicity_3/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wot_multiplicity_4/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wot_multiplicity_5/performance_ckf.root"
]

file_paths_wt = [
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wt_multiplicity_1/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wt_multiplicity_2/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wt_multiplicity_3/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wt_multiplicity_4/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wt_multiplicity_5/performance_ckf.root"
]

# Build a TCanvas
canvas = ROOT.TCanvas("canvas", "Scatter Plot Comparison", 800, 600)

multiplicity = [1, 2, 3, 4, 5]
te_wotime = []
te_wtime = []
TEErr_WoT = []
TEErr_WT = []

# Build WO-Time TGraph

for file_path in file_paths_wot:
    file = ROOT.TFile(file_path, "READ")
    tefficiency_name = "perfSummary"
    tefficiency = file.Get(tefficiency_name)
    value = tefficiency.GetEfficiency(1)
    print(value)
    teerr_wot = (tefficiency.GetEfficiencyErrorLow(1)+tefficiency.GetEfficiencyErrorUp(1))/2
    te_wotime.append(value)
    TEErr_WoT.append(teerr_wot)
    file.Close()

n_points = len(multiplicity)

graph1 = ROOT.TGraphErrors(n_points)

for i, (x, y) in enumerate(zip(multiplicity, te_wotime)):
    graph1.SetPoint(i, x, y)
for i, (e) in enumerate(zip(TEErr_WoT)):
    graph1.SetPointError(i, 0, TEErr_WoT[i-1])

graph1.SetMarkerStyle(ROOT.kFullCircle)
graph1.SetMarkerSize(0.6)
graph1.SetMarkerColor(ROOT.kBlue)
graph1.GetXaxis().SetRangeUser(0, 6)
graph1.GetYaxis().SetRangeUser(0.7, 1.2)
graph1.SetLineStyle(1)
graph1.SetLineColor(ROOT.kBlue)
graph1.SetTitle("Tracking Efficiency Compare about Multiplicity")
graph1.Draw("APL")  # "A" as Axis，"P" as Plot

# Build W-Time TGraph

for file_path in file_paths_wt:
    file = ROOT.TFile(file_path, "READ")
    tefficiency_name = "perfSummary"
    tefficiency = file.Get(tefficiency_name)
    value = tefficiency.GetEfficiency(1)
    print(value)
    teerr_wt = (tefficiency.GetEfficiencyErrorLow(1)+tefficiency.GetEfficiencyErrorUp(1))/2
    te_wtime.append(value)
    TEErr_WT.append(teerr_wt)
    file.Close()
n_points = len(multiplicity)

graph2 = ROOT.TGraphErrors(n_points)

for i, (x, y) in enumerate(zip(multiplicity, te_wtime)):
    graph2.SetPoint(i, x, y)
for i, (e) in enumerate(zip(TEErr_WT)):
    graph2.SetPointError(i, 0, TEErr_WT[i-1])

graph2.SetMarkerStyle(ROOT.kFullSquare)
graph2.SetMarkerSize(0.6)
graph2.SetMarkerColor(ROOT.kRed)
graph2.GetXaxis().SetRangeUser(0, 6)
graph2.GetYaxis().SetRangeUser(0.7, 1.2)
graph2.SetLineStyle(1)
graph2.SetLineColor(ROOT.kRed)
graph2.Draw("PL same") 

# Axis Title
graph1.GetXaxis().SetTitle("Multiplicity")
graph1.GetYaxis().SetTitle("Tracking Efficiency")

# Legend
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(graph1, "WO-Time", "p")
legend.AddEntry(graph2, "W-Time", "p")
legend.Draw()

canvas.Update()
canvas.Draw()
# Title
canvas.SetTitle("Tracking Efficiency Compare about Multiplicity")

canvas.SaveAs("TE-Multiplicity.pdf")

