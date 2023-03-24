#code that open vol_layer.csv 
#manually add Dmitry info to the csv before 
#makes histogram 
#defines dmitry colour scale 
import pandas as pd 
import ROOT
from ROOT import TFile, TTree, gStyle , TColor 
import numpy as np 
import csv 
import array as arr 


#open csv 

df = pd.read_csv('vol_layer_morestats.csv', skiprows=1)
array = df.to_numpy()

acts_vol = array[:,0]
acts_lay = array[:,1]
z = array[:,2]
r =array[:,3]
FTF_id= array[:,4]


#make histogram 
c1 = ROOT.TCanvas("c1","c1",800,600)
h1 = ROOT.TH2D("h1","ACTS volume & layers in terms of FTF volume id;z;r",100,-3000,3000,100,0,1000) 

# #ACTS: 
# for i in range(0,len(acts_vol)): 

#     h1.Fill(z[i],r[i],acts_vol[i])

#FTF
for i in range(0,len(acts_vol)): 
    if FTF_id[i] != 0: 
        h1.Fill(z[i],r[i],FTF_id[i])
        

#zcontours = [70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,90,91,92,93,94,95,96,97,98]
zcontours = [70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,85,86,87,88,89,90,91,92,93,94,95,96,97,98]
zcountour = arr.array('d', zcontours)

h1.SetContour(len(zcontours),zcountour)

h1.SetMinimum(70)
h1.SetMaximum(98) 

        

gStyle.SetPalette(91)

c1.cd()
h1.SetStats(0)
h1.Draw("COLZ") #text if want numbers 
c1.SaveAs("vol_layer_FTF.pdf")
