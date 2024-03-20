import ROOT

# 文件路径列表
file_paths = [
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wot_multiplicity_1/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wot_multiplicity_2/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wot_multiplicity_3/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wot_multiplicity_4/performance_ckf.root",
    "/home/ZhangHaolin/ACTS/acts-telescope/RunSpace/wot_multiplicity_5/performance_ckf.root"
]

# Build a TCanvas
canvas = ROOT.TCanvas("canvas", "Scatter Plot Comparison", 800, 600)

# Build WO-Time TGraph
multiplicity = [1, 2, 3, 4, 5]
te_wotime = []
# te_wotime = [1, 1, 0.944756, 0.89689, 0.81487]

# 遍历每个文件
for file_path in file_paths:
    # 打开 ROOT 文件
    file = ROOT.TFile(file_path, "READ")
    # 获取 TEfficiency
    tefficiency_name = "perfSummary"  # 替换为你实际的 TEfficiency 名称
    tefficiency = file.Get(tefficiency_name)
    # 获取第一个 bin 的值
    value = tefficiency.GetEfficiency(1)  # 假设你想读取第一个 bin 的值
    te_wotime.append(value)
    file.Close()

n_points = len(multiplicity)
graph1 = ROOT.TGraph(n_points)

for i, (x, y) in enumerate(zip(multiplicity, te_wotime)):
    graph1.SetPoint(i, x, y)

graph1.SetMarkerStyle(ROOT.kFullCircle)
graph1.SetMarkerSize(1.0)
graph1.SetMarkerColor(ROOT.kBlue)
graph1.GetYaxis().SetRangeUser(0.7, 1.2)
graph1.SetLineStyle(1)
graph1.SetLineColor(ROOT.kBlue)
graph1.SetTitle("Tracking Efficiency Compare about Multiplicity")
graph1.Draw("APL")  # "A" as Axis，"P" as Plot

# Build W-Time TGraph
te_wtime = [1, 1, 0.943171, 0.904585, 0.813911]
graph2 = ROOT.TGraph(n_points)

for i, (x, y) in enumerate(zip(multiplicity, te_wtime)):
    graph2.SetPoint(i, x, y)

graph2.SetMarkerStyle(ROOT.kFullTriangleUp)
graph2.SetMarkerSize(1.0)
graph2.SetMarkerColor(ROOT.kRed)
graph2.GetYaxis().SetRangeUser(0.7, 1.2)
graph2.SetLineStyle(1)
graph2.SetLineColor(ROOT.kRed)
graph2.Draw("PL same")  # "same"

# Axis Title
graph1.GetXaxis().SetTitle("Multiplicity")
graph1.GetYaxis().SetTitle("Tracking Efficiency")

# Legend
legend = ROOT.TLegend(0.7, 0.7, 0.9, 0.9)
legend.AddEntry(graph1, "WO-Time", "p")
legend.AddEntry(graph2, "W-Time", "p")
legend.Draw()

# Title
canvas.SetTitle("Tracking Efficiency Compare about Multiiplicity")

canvas.SaveAs("TE-Multiplicity.pdf")

