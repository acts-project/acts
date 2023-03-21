import ROOT
from ROOT import TFile, TTree 
import numpy as np 


inFile = ROOT.TFile.Open("itk_output/trackstates_ckf.root","READ")

tree = inFile.Get("trackstates")

c1 = ROOT.TCanvas("c1","c1",800,600)
h1 = ROOT.TH2D("h1","test",10000,-3000,3000,10000,0,1200) 


#looping over entries, each have arrays of x,y and z
for entry in range(0,tree.GetEntries()): 
    tree.GetEntry(entry) 
    vol_id = getattr(tree, "volume_id")
    #print(vol_id)
    #this is not a single value its a list of them 
    x = getattr(tree,"g_x_hit")
    y = getattr(tree,"g_y_hit")
    z = getattr(tree,"g_z_hit")
    #looping over the hits in each entry 
    for i in range(0,len(x)): 
        r = np.sqrt(x[i]**2 + y[i]**2) 
        h1.Fill(z[i],r)        
        
    #each entry has different volume numbers       
        # for j in range(0,25): 
        #     if vol_id[i] == j: 
        # #if len(x) != 0:  
        #         print("vol Id= ", j , "r = " ,r , "z= ", z[i] )



# inFile.Close() 

# outHistFile = TFile('output.root','RECREATE')
# # outHistFile.cd ()
# #outHistFile["h1"] = h1 
# h1.Write() 

# outHistFile.Close() 

#save as image instead, used this in the dimuon code? 

c1.cd()
h1.Draw() 
c1.SaveAs("Python_output.jpg")