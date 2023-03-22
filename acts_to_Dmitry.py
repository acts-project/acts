import csv 

f = open('acts_to_Dmitry.csv','w')
writer = csv.writer(f) 



header = ["ACTS_vol_id","ACTS_layer_id","Dmitry"]
writer.writerow(header)

for j in range(2,20,2): 
    row = [15,j,78] 
    writer.writerow(row)

for j in range(20,38,2): 
    row = [15,j,77] 
    writer.writerow(row)

#doesnt work had to edit 
for i in range(80,81): #dmirty number 
    for j in range(2,6,2): #layer number 
        #vol nuumber is 16 
        row = [9,j,i] 
        writer.writerow(row)

for i in range(82,85,1): #dmirty number 
    for j in range(2,6,2): #layer number 
        #vol nuumber is 16 
        row = [16,j,i] 
        writer.writerow(row)

for j in range(2,20,2): 
    row = [20,j,97] 
    writer.writerow(row)

for j in range(20,38,2): 
    row = [20,j,98] 
    writer.writerow(row)






# row = [ ]
# writer.writerow(row)

f.close() 
 
#open file to read csv then find right row? 
# loop over rows, if row [0] = and row [1] equal return row[3]  
#load file of ACTS volumes and layers- 