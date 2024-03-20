import ROOT

# 文件路径列表
file_paths_wot = [
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wot_pos_stddev_0/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wot_pos_stddev_20/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wot_pos_stddev_40/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wot_pos_stddev_60/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wot_pos_stddev_80/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wot_pos_stddev_100/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wot_pos_stddev_120/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wot_pos_stddev_140/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wot_pos_stddev_160/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wot_pos_stddev_180/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wot_pos_stddev_200/performance_ckf.root",
]

file_paths_wt = [
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wt_pos_stddev_0/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wt_pos_stddev_20/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wt_pos_stddev_40/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wt_pos_stddev_60/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wt_pos_stddev_80/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wt_pos_stddev_100/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wt_pos_stddev_120/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wt_pos_stddev_140/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wt_pos_stddev_160/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wt_pos_stddev_180/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wt_pos_stddev_200/performance_ckf.root",
]

# Build a TCanvas
canvas = ROOT.TCanvas("canvas", "Scatter Plot Comparison", 800, 600)

p_stddev = [0, 20, 40, 60, 80, 100, 120, 140, 160, 180, 200]
fr_wotime = []
fr_wtime = []
FRErr_WoT = []
FRErr_WT = []

# Build WO-Time TGraph

for file_path in file_paths_wot:
    file = ROOT.TFile(file_path, "READ")
    tefficiency_name = "perfSummary"
    tefficiency = file.Get(tefficiency_name)
    value = tefficiency.GetEfficiency(2)
    print(value)
    frerr_wot = (tefficiency.GetEfficiencyErrorLow(2)+tefficiency.GetEfficiencyErrorUp(2))/2
    print(frerr_wot)
    fr_wotime.append(value)
    FRErr_WoT.append(frerr_wot)
    file.Close()

n_points = len(p_stddev)

graph1 = ROOT.TGraphErrors(n_points)

for i, (x, y) in enumerate(zip(p_stddev, fr_wotime)):
    graph1.SetPoint(i, x, y)
for i, (e) in enumerate(zip(FRErr_WoT)):
    graph1.SetPointError(i, 0, FRErr_WoT[i-1])

graph1.SetMarkerStyle(ROOT.kFullCircle)
graph1.SetMarkerSize(0.6)
graph1.SetMarkerColor(ROOT.kBlue)
graph1.GetXaxis().SetRangeUser(0, 200)
graph1.GetYaxis().SetRangeUser(0, 0.3)
graph1.SetLineStyle(1)
graph1.SetLineColor(ROOT.kBlue)
graph1.SetTitle("Fake Rate Compare about Position StdDev")
graph1.Draw("APL")  # "A" as Axis，"P" as Plot

# Build W-Time TGraph
for file_path in file_paths_wt:
    file = ROOT.TFile(file_path, "READ")
    tefficiency_name = "perfSummary"
    tefficiency = file.Get(tefficiency_name)
    value = tefficiency.GetEfficiency(2)
    print(value)
    frerr_wt = (tefficiency.GetEfficiencyErrorLow(2)+tefficiency.GetEfficiencyErrorUp(2))/2
    print(frerr_wt)
    fr_wtime.append(value)
    FRErr_WT.append(frerr_wt)
    file.Close()
n_points = len(p_stddev)

graph2 = ROOT.TGraphErrors(n_points)

for i, (x, y) in enumerate(zip(p_stddev, fr_wtime)):
    graph2.SetPoint(i, x, y)
for i, (e) in enumerate(zip(FRErr_WT)):
    graph2.SetPointError(i, 0, FRErr_WT[i-1])

graph2.SetMarkerStyle(ROOT.kFullSquare)
graph2.SetMarkerSize(0.6)
graph2.SetMarkerColor(ROOT.kRed)
graph2.GetXaxis().SetRangeUser(0, 200)
graph2.GetYaxis().SetRangeUser(0, 0.3)
graph2.SetLineStyle(1)
graph2.SetLineColor(ROOT.kRed)
graph2.Draw("PL same") 

# Axis Title
graph1.GetXaxis().SetTitle("Position StdDev")
graph1.GetYaxis().SetTitle("Fake Rate")

# Legend
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(graph1, "WO-Time", "p")
legend.AddEntry(graph2, "W-Time", "p")
legend.Draw()

canvas.Update()
canvas.Draw()
# Title
canvas.SetTitle("Fake Rate Compare about Position StdDev")

canvas.SaveAs("FR-Pos-StdDev.pdf")

