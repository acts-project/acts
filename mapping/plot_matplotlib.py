
import pandas as pd 
import numpy as np 
import csv 
import array as arr 
import matplotlib.pyplot as plt 
import matplotlib.colors 


#open csv 

df = pd.read_csv('FTF_SP_output.csv', skiprows=0)
array = df.to_numpy()

FTF_id= array[:,0]
z = array[:,1]
r =array[:,2]
eta_mod =array[:,3]



# cm = plt.cm.get_cmap('gist_rainbow') #eta mod 
# cm = plt.cm.get_cmap('spring') #ftf 

cm = matplotlib.colors.LinearSegmentedColormap.from_list("", ["mediumvioletred","palevioletred","thistle","lavender","lightsteelblue","cornflowerblue"]) 

sc = plt.scatter(z,r,c=FTF_id, cmap=cm, marker = '.', s = 2) 
plt.colorbar(sc)
plt.xlabel('z')
plt.ylabel('r')
plt.title("Space points in terms of FTF Layer id")
plt.show() 
plt.savefig('./mapping/FTF_spacepoints_plt.pdf')


# #modules 
# df = pd.read_csv('mapping/Module_plot.csv', skiprows=1)
# array = df.to_numpy()

# ACTS_vol= array[:,0]
# z = array[:,3]
# r =array[:,4]
# FTF_id =array[:,5]


# # cm = plt.cm.get_cmap('cool') #acts 
# cm = matplotlib.colors.LinearSegmentedColormap.from_list("", ["mediumvioletred","palevioletred","thistle","lavender","lightsteelblue","cornflowerblue"]) 
# # cm = plt.cm.get_cmap('spring') #ftf 
# sc = plt.scatter(z,r,c=FTF_id, cmap=cm, marker = '.', s = 2) 
# plt.colorbar(sc)
# plt.clim(70,98)
# plt.xlabel('z')
# plt.ylabel('r')
# plt.xlim(-3000,3000)
# plt.ylim(0,350)
# plt.title("Pixel geometry in terms of FTF Layer id")
# plt.show() 
# plt.savefig('./mapping/modules_FTF_plt.pdf')