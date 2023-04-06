import ROOT
from ROOT import TFile, TTree, gStyle 
import numpy as np 
import csv 

# making histogram: 
# at top 
c1 = ROOT.TCanvas("c1","c1",800,600)
h2 = ROOT.TH2D("h2","CKF hits in Itk in terms of ACTS volume id",8000,-3000,3000,8000,0,1200) 

#open tree 
inFile = ROOT.TFile.Open("itk_output/trackstates_ckf.root","READ")
tree = inFile.Get("trackstates")
# f = open('vol_layer_morestats.csv','w')
# writer = csv.writer(f) 


geometry = {}
#looping over entries, each have arrays of x,y and z
for entry in range(0,10000):#tree.GetEntries()): 
    tree.GetEntry(entry) 
    vol_id = getattr(tree, "volume_id")
    lay_id = getattr(tree, "layer_id")
    #this is not a single value its a list of them 
    x = getattr(tree,"g_x_hit")
    y = getattr(tree,"g_y_hit")
    z = getattr(tree,"g_z_hit")
    #looping over the hits in each entry 

    #make an empty dictionary  
     

    for i in range(0,len(x)): 
        r = np.sqrt(x[i]**2 + y[i]**2) 
        #make single value of r for each hit, just use as r as currently in hit loop not event loop 
        #h1.Fill(z[i],r)        
        
    # each entry has different volume numbers       
        for j in range(2,25): 
            # use i to access vol and layer number of the hit 
            if vol_id[i] == j: 
                for k in range(2,60): 
                    #60 layer numbers (not all of them)
                    if lay_id[i] == k: 
                        #if len(x) != 0:  #check are hits, not sure need this, may speed it up 
                        #check if currently in dictionary 
                        key = str(j)+ "_" + str(k) #current vol and layer, turn into string 

                        if key in geometry: 
                            #need to add to z,r,n 
                            geometry[key][0] += z[i]  #alter z 
                            geometry[key][1] +=  r #alter r 
                            geometry[key][2] += 1 #increment n 

                        else:
                            #fill dictionary normally 
                            n=1 
                            geometry[key] = [z[i],r,n]


                        
#average values 
for key in geometry: 
    geometry[key][0] /= geometry[key][2] #average z 
    geometry[key][1] /= geometry[key][2] #average r 
    
print(geometry)
print(len(geometry)) 

#writing to csv 

header = ["vol_id","layer_id","z","r"]
writer.writerow(header)

for i in geometry: 
    
    #retrive vol and layer integers for csv 
    split = i.split("_") 
    vol = int(split[0])
    layer = int(split[1]) 
    #add to csv here 
    row = [vol,layer,geometry[i][0],geometry[i][1] ]
    writer.writerow(row)

f.close() 


#making histogram: 
#at top 
# c1 = ROOT.TCanvas("c1","c1",800,600)
# h1 = ROOT.TH2D("h1","test",10000,-3000,3000,10000,0,1200) 

# inFile.Close() 

# c1.cd()
# h1.Draw() 
# c1.SaveAs("Python_output.jpg")
