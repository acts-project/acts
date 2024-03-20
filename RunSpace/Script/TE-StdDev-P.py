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
    print(teerr_wot)
    te_wotime.append(value)
    TEErr_WoT.append(teerr_wot)
    file.Close()

n_points = len(p_stddev)

graph1 = ROOT.TGraphErrors(n_points)

for i, (x, y) in enumerate(zip(p_stddev, te_wotime)):
    graph1.SetPoint(i, x, y)
for i, (e) in enumerate(zip(TEErr_WoT)):
    graph1.SetPointError(i, 0, TEErr_WoT[i-1])

graph1.SetMarkerStyle(ROOT.kFullCircle)
graph1.SetMarkerSize(0.6)
graph1.SetMarkerColor(ROOT.kBlue)
graph1.GetXaxis().SetRangeUser(0, 200)
graph1.GetYaxis().SetRangeUser(0.7, 1.2)
graph1.SetLineStyle(1)
graph1.SetLineColor(ROOT.kBlue)
graph1.SetTitle("Tracking Efficiency Compare about Position StdDev")
graph1.Draw("APL")  # "A" as Axis，"P" as Plot

# Build W-Time TGraph
for file_path in file_paths_wt:
    file = ROOT.TFile(file_path, "READ")
    tefficiency_name = "perfSummary"
    tefficiency = file.Get(tefficiency_name)
    value = tefficiency.GetEfficiency(1)
    print(value)
    teerr_wt = (tefficiency.GetEfficiencyErrorLow(1)+tefficiency.GetEfficiencyErrorUp(1))/2
    print(teerr_wt)
    te_wtime.append(value)
    TEErr_WT.append(teerr_wt)
    file.Close()
n_points = len(p_stddev)

graph2 = ROOT.TGraphErrors(n_points)

for i, (x, y) in enumerate(zip(p_stddev, te_wtime)):
    graph2.SetPoint(i, x, y)
for i, (e) in enumerate(zip(TEErr_WT)):
    graph2.SetPointError(i, 0, TEErr_WT[i-1])

graph2.SetMarkerStyle(ROOT.kFullSquare)
graph2.SetMarkerSize(0.6)
graph2.SetMarkerColor(ROOT.kRed)
graph2.GetXaxis().SetRangeUser(0, 200)
graph2.GetYaxis().SetRangeUser(0.7, 1.2)
graph2.SetLineStyle(1)
graph2.SetLineColor(ROOT.kRed)
graph2.Draw("PL same") 

# Axis Title
graph1.GetXaxis().SetTitle("Position StdDev")
graph1.GetYaxis().SetTitle("Tracking Efficiency")

# Legend
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(graph1, "WO-Time", "p")
legend.AddEntry(graph2, "W-Time", "p")
legend.Draw()

canvas.Update()
canvas.Draw()
# Title
canvas.SetTitle("Tracking Efficiency Compare about Position StdDev")

canvas.SaveAs("TE-Pos-StdDev.pdf")

