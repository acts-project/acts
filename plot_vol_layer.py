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
h1 = ROOT.TH2D("h1","Pixel geometry in terms of FTF layer id;z;r",100,-3000,3000,100,0,350) 
# h1 = ROOT.TH2D("h1","ITK geometry in terms of ACTS volume id;z;r",100,-3000,3000,100,0,1000) 


# # #ACTS: 
# for i in range(0,len(acts_vol)): 

#     h1.Fill(z[i],r[i],acts_vol[i])

# FTF
for i in range(0,len(acts_vol)): 
    if FTF_id[i] != 1: 
        h1.Fill(z[i],r[i],FTF_id[i])
    if FTF_id[i] == 1: 
        print("unidentified ", "volume:" , acts_vol[i] ,"layer: ", acts_lay[i])

    
        

# #zcontours = [70,71,72,73,74,75,76,77,78,79,80,81,82,83,84,90,91,92,93,94,95,96,97,98]
zcontours = [0,70.5,71.5,72.5,73.5,74.5,75.5,76.5,77.5,78.5,79.5,80.5,81.5,82.5,83.5,84.5,85.5,86.5,87.5,88.5,89.5,90.5,91.5,92.5,93.5,94.5,95.5,96.5,97.5,98.5]
zcountour = arr.array('d', zcontours)

h1.SetContour(len(zcontours),zcountour)

h1.SetMinimum(70)
h1.SetMaximum(98) 

        

gStyle.SetPalette(62)
#gStyle.SetPalette(91) 

c1.cd()
h1.SetStats(0)
h1.Draw("COLZ") #text if want numbers 
c1.SaveAs("vol_layer_FTF.pdf")
