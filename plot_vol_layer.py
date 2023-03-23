#code that open vol_layer.csv 
#manually add Dmitry info to the csv before 
#makes histogram 
#defines dmitry colour scale 
import pandas as pd 
import ROOT
from ROOT import TFile, TTree, gStyle , TColor 
import numpy as np 
import csv 


#open csv 

df = pd.read_csv('vol_layer.csv')
array = df.to_numpy()

acts_vol = array[:,0]
acts_lay = array[:,1]
z = array[:,2]
r =array[:,3]
FTF_id= array[:,4]


#make histogram 
c1 = ROOT.TCanvas("c1","c1",800,600)
h1 = ROOT.TH2D("h1","vol_layer;z;r",10000,-3000,3000,10000,-3000,3000) 


for i in range(0,len(acts_vol)): 
    if FTF_id[i] != 0: 
        h1.Fill(z[i],r[i],FTF_id[i])
    # else: 
    #     h1.Fill(z[i],r[i])



stops = [0.00, 0.34, 0.61, 0.84, 1.00]
red   = [0.00, 0.00, 0.87, 1.00, 0.51]
green = [0.00, 0.81, 1.00, 0.20, 0.00]
blue  = [0.51, 1.00, 0.12, 0.00, 0.00]

TColor.CreateGradientColorTable(7, stops, red, green, blue, 999)
gStyle.SetNumberContours(999)



#gStyle.SetPalette(80)

c1.cd()
h1.SetStats(0)
h1.Draw("COLZ") 
c1.SaveAs("vol_layer_better.pdf")
