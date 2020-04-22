import pandas as pd
import numpy as np
import func

#MAIN
data = pd.read_csv('GSE64881_segmentation_at_30000bp.passqc.multibam.txt', sep='\t', header = None)
array = np.array(data)

feature_data = pd.read_csv('gam_feature_community.csv')

names = data.iloc[0,3:].values
nums = data.iloc[1:,0:3].values
alldata = data.iloc[1:,3:].values

colcount = (alldata == 1).sum(axis = 0)

rowcount = (alldata == 1).sum(axis = 1)
row100 = np.argwhere(np.array(rowcount).flatten() < 100).flatten()
rowcount = rowcount[row100]

histrange = data.iloc[69715:69796,3:].values
histrangenums = data.iloc[69715:69796,1:3].values

while (1):
    print("\n\n\n\n\n\nMain Menu:\nCo-segregation(1)\nNetwork Centrality(2)\nCommunity Detection(3)\nQuit(q)")
    val = input("")

    if val == '1' or val == '2' or val == '3': #General Statistics
        func.co_segregation(histrange,histrangenums, val, feature_data)

    elif val == 'q' or val == 'Q': #Quit
        break

    else:
        print("Invalid Input")


    input("\nPress Enter to continue...")