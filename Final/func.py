import random
import math
import pandas as pd
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
from datetime import datetime
from functools import reduce

def co_segregation (histrange,histrangenums, val, feature_data):
    total = len(histrange[1])
    histrowcount = (histrange == 1).sum(axis = 1)
    
    detection = []
    for x in histrowcount:
        detection.append(round((x/total),2))

    Fab = []
    for x in histrange:
        row = []
        for y in histrange:
            temp1 = np.argwhere(x == 1)
            temp2 = np.argwhere(y == 1)
            intersect = np.intersect1d(temp1,temp2)
            row.append(round((len(intersect)/total),2))
        Fab.append(row)

    linkage = []
    for x in range(len(detection)):
        row = []
        for y in range(len(detection)):
            num = Fab[x][y] - (detection[x] * detection[y])
            row.append(num)
        linkage.append(row)

    normlinkage = []
    for x in range(len(detection)):
        row = []
        for y in range(len(detection)):
            if linkage[x][y] > 0:
                Dmax = min((detection[y]*(1-detection[x]),detection[x]*(1-detection[y])))
                if Dmax == 0:
                    row.append(0)
                else: 
                    row.append(linkage[x][y]/Dmax)

            elif linkage[x][y] < 0:
                Dmax = min(Fab[x][y],(1 - detection[x]) * (1 - detection[y]))
                if Dmax == 0:
                    row.append(0)
                else: 
                    row.append(linkage[x][y]/Dmax)
            else:
                row.append(0)
        normlinkage.append(row)

    
    if val == '1':
        df = pd.DataFrame(detection)
        df.to_csv("./files/detection.csv")

        df = pd.DataFrame(Fab)
        df.to_csv("./files/Fab.csv")

        df = pd.DataFrame(linkage)
        df.to_csv("./files/linkage.csv")

        df = pd.DataFrame(normlinkage)
        df.to_csv("./files/normlinkage.csv")

        plt.figure(figsize=(50,50))
        plt.title("Normilized Linkage", size = 50)
        heat_map = sb.heatmap(normlinkage, square = True, cbar_kws=({"shrink": .80}))
        plt.savefig("./files/normlinkage.png")
        plt.clf()
    elif val == '2' or val == '3':
        NetworkCentrality(normlinkage, val, feature_data)

    else:
        return


def NetworkCentrality(normlinkage, val, feature_data):

    Lavg = np.average(normlinkage)

    totalcentrality = []
    for x in range(len(normlinkage[0])):
        cantralityx = []
        for y in range(len(normlinkage[0])):
            if(x != y):
                if normlinkage[x][y] > Lavg:
                    cantralityx.append(y)
        totalcentrality.append(cantralityx)

    average = 0
    minx = 500
    minnum = -1
    maxx = 0
    maxnum = -1

    for x in range(len(totalcentrality)):
        lenx = len(totalcentrality[x])
        if lenx < minx:
            minx = lenx
            minnum = x
        elif lenx > maxx:
            maxx = lenx
            maxnum = x
        average += lenx

    if val == '2':
        print("\n\nIndex - Edges\n")
        print("-----------------------------------------------------------------------------\n")

        for num in range(73):
            temp = num
            for x in range(len(totalcentrality)):    
                if(len(totalcentrality[x]) == temp):
                    print(" ", '%2s' % x,"  - ", totalcentrality[x], end =" ")
                    print("\n ", len(totalcentrality[x]), "/", (len(normlinkage[0]) - 1))
                    print("\n")



        print("\nStatistics\n")
        print("Average")
        print(average/len(totalcentrality), "/", (len(normlinkage[0]) - 1))
        print("\nMinimum Index - Value")
        print(minnum, "-", minx, "/", (len(normlinkage[0]) - 1))
        print("\nMaximum Index - Value")
        print(maxnum, "-", maxx, "/", (len(normlinkage[0]) - 1))

    elif val == '3':
        order = []
        for num in range(73):
            temp = num
            for x in range(len(totalcentrality)):    
                if(len(totalcentrality[x]) == temp):
                    order.append([x, totalcentrality[x]])
        Community_Detection(order, feature_data)
  
    else:
         return


def Community_Detection(order, feature_data):
    hist1 = feature_data.iloc[:,13]
    lad = feature_data.iloc[:,9]
    
    histindex = np.argwhere(hist1 == 1)
    ladindex = np.argwhere(lad == 1)

    top5 = []
    for x in range(len(order)):
        if x > len(order) - 6:
            top5.append(order[x])

    fullscale = []
    for x in top5:
        withzero = []
        for y in range(81):
            hit = False
            for z in x[1]:
                if y == z :
                    hit = True
                    withzero.append(1)
            if hit == False:
                withzero.append(0)
        fullscale.append([str(x[0]),withzero])
    fullscale.append(["HIST1",hist1])
    fullscale.append(["LAD",lad])


    print("\nReport\n")
    print("---------------------------------------------------------")
    for x in range(len(top5)):
        print("Hub ", top5[x][0])
        print("Hub size:")
        print(len(top5[x][1]))

        print("\nHub Nodes:")
        print(top5[x][1])

        print("\nHub HIST1%:")
        intersect = np.intersect1d(histindex,top5[x][1])
        print(len(intersect)/16)

        print("\nHub LAD%:")
        intersect = np.intersect1d(ladindex,top5[x][1])
        print(len(intersect)/38)
        print("\n-------------------------------------------------------------\n")


    total = []
    index = []
    for x in fullscale:
        total.append(x[1])
        index.append(x[0])

    plt.figure(figsize=(25,3))
    plt.title("Community Detection", size = 30)
    heat_map = sb.heatmap(total, yticklabels = index, cbar_kws=({"shrink": .80}))
    plt.savefig("./files/Community_Detection.png")
    plt.clf()
    return