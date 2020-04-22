import random
import math
import pandas as pd
import numpy as np
import seaborn as sb
import matplotlib.pyplot as plt
from datetime import datetime
from functools import reduce

def stats (data): #Gets the Data Number of NP by Genomic windows
    print("\nNumber of genomic windows: ")
    print(len(data.index) - 1)

    print("\nNumber of NPs: ")
    print(len(data.columns) - 3)

def colstats (colcount): #Takes the values of geomic windows like min, max, and average
    print("\nCol Average NP detected: ")
    print (np.average(colcount))

    print("\nCol Max NP detected: ")
    print (max(colcount))

    print("\nCol Min NP detected: ")
    print (min(colcount))

def rowstats (rowcount): #Takes the values of NPs like min, max, and average
    print("\nRow Average NP detected: ")
    print (np.average(rowcount))

    print("\nRow Max NP detected: ")
    print (max(rowcount))

    print("\nRow Min NP detected: ")
    print (min(rowcount))

def colposition (colcount, num, names): #sets all col to a value 1-10

    if num == '1':
        precent10 = np.argwhere(colcount < np.percentile(colcount, 10))
        print(names[precent10])

    elif num == '2':   
        precent20 = np.intersect1d(np.argwhere(colcount >= np.percentile(colcount, 10)),np.argwhere(colcount < np.percentile(colcount, 20)))
        print(names[precent20])

    elif num == '3':
        precent30 = np.intersect1d(np.argwhere(colcount >= np.percentile(colcount, 20)),np.argwhere(colcount < np.percentile(colcount, 30)))
        print(names[precent30])

    elif num == '4':
        precent40 = np.intersect1d(np.argwhere(colcount >= np.percentile(colcount, 30)),np.argwhere(colcount < np.percentile(colcount, 40)))    
        print(names[precent40])

    elif num == '5': 
        precent50 = np.intersect1d(np.argwhere(colcount >= np.percentile(colcount, 40)),np.argwhere(colcount < np.percentile(colcount, 50)))
        print(names[precent50])

    elif num == '6':
         precent60 = np.intersect1d(np.argwhere(colcount >= np.percentile(colcount, 50)),np.argwhere(colcount < np.percentile(colcount, 60)))
         print(names[precent60])

    elif num == '7': 
        precent70 = np.intersect1d(np.argwhere(colcount >= np.percentile(colcount, 60)),np.argwhere(colcount < np.percentile(colcount, 70)))
        print(names[precent70])

    elif num == '8':
        precent80 = np.intersect1d(np.argwhere(colcount >= np.percentile(colcount, 70)),np.argwhere(colcount < np.percentile(colcount, 80)))
        print(names[precent80])

    elif num == '9': 
        precent90 = np.intersect1d(np.argwhere(colcount >= np.percentile(colcount, 80)),np.argwhere(colcount < np.percentile(colcount, 90)))
        print(names[precent90])

    elif num == '10': 
        precent100 = np.argwhere(colcount > np.percentile(colcount, 90))
        print(names[precent100])

    else:
        print("Invalid Number")

def rowposition (rowcount, num, nums): #sets all rows to a value 1-5
    file = open("files/Indexrows1.txt","w") 
    
    numnames = []
    for x in nums:
        numnames.append(str(x[0]) + " " + str(x[1]) + " - " + str(x[2]))

    if num == '1':
        file = open("files/Indexrows1.txt","w") 
        precent20 = np.argwhere(rowcount < np.percentile(rowcount, 20))
        for x in precent20:
            file.write(numnames[int(x)])
            file.write("\n")
        file.close() 

    elif num == '2':
        file = open("files/Indexrows2.txt","w") 
        precent40 = np.intersect1d(np.argwhere(rowcount >= np.percentile(rowcount, 20)),np.argwhere(rowcount < np.percentile(rowcount, 40)))
        for x in precent40:
            file.write(numnames[int(x)])
            file.write("\n")
        file.close() 

    elif num == '3': 
        file = open("files/Indexrows3.txt","w") 
        precent60 = np.intersect1d(np.argwhere(rowcount >= np.percentile(rowcount, 40)),np.argwhere(rowcount < np.percentile(rowcount, 60)))
        for x in precent60:
            file.write(numnames[int(x)])
            file.write("\n")
        file.close() 

    elif num == '4': 
        file = open("files/Indexrows4.txt","w") 
        precent80 = np.intersect1d(np.argwhere(rowcount >= np.percentile(rowcount, 60)),np.argwhere(rowcount < np.percentile(rowcount, 80)))
        for x in precent80:
            file.write(numnames[int(x)])
            file.write("\n")
        file.close() 

    elif num == '5':
        file = open("files/Indexrows5.txt","w")    
        precent100 = np.argwhere(rowcount > np.percentile(rowcount, 80))
        for x in precent100:
            file.write(numnames[int(x)])
            file.write("\n")
        file.close() 

    else:
        print("Invalid Number")

def histone(histrange, rowcount, colcount):  #Gives values about the Hist genome
    histcolcount = (histrange == 1).sum(axis = 0)
    histrowcount = (histrange == 1).sum(axis = 1)
    non0histcolcount = []

    for x in histcolcount:
        if np.sum(x) > 0:
            non0histcolcount.append(x)
    histcolcount = non0histcolcount

    print("\nNumber of genomic windows")
    print(len(histrowcount)) 

    print("\nNumber of NPs")
    print(len(histcolcount))

    print("\nOn average, how many windows are present in an NP")
    print (np.average(histrowcount))

    print("\nMax number of windows present in any NP")
    print (max(histrowcount))

    print("\nMin number of windows present in any NP")
    print (min(histrowcount))

    print("\nAverage number of NPs in which a window is detected")
    print (np.average(histcolcount))

    print("\nMax number of NPs in which a window is detected")
    print (max(histcolcount))

    print("\nMin number of NPs in which a window is detected")
    print (min(histcolcount))

    precent20 = np.argwhere(rowcount < np.percentile(rowcount, 20))
    precent40 = np.intersect1d(np.argwhere(rowcount >= np.percentile(rowcount, 20)),np.argwhere(rowcount < np.percentile(rowcount, 40)))
    precent60 = np.intersect1d(np.argwhere(rowcount >= np.percentile(rowcount, 40)),np.argwhere(rowcount < np.percentile(rowcount, 60)))
    precent80 = np.intersect1d(np.argwhere(rowcount >= np.percentile(rowcount, 60)),np.argwhere(rowcount < np.percentile(rowcount, 80)))   
    precent100 = np.argwhere(rowcount > np.percentile(rowcount, 80))

    total = 0
    for x in histrowcount:
        if x < max(rowcount[precent20]):
            total += 1
        if x >= max(rowcount[precent20]) and x < max(rowcount[precent40]):
            total += 2
        if x >= max(rowcount[precent40]) and x < max(rowcount[precent60]):
            total += 3
        if x >= max(rowcount[precent60]) and x < max(rowcount[precent80]):
            total += 4
        if x >= max(rowcount[precent80]):
            total += 5
    print("\nOn Average radial position:")
    print(total/len(histrowcount))

    precent10 = np.argwhere(colcount < np.percentile(colcount, 10))
    precent20 = np.intersect1d(np.argwhere(colcount >= np.percentile(colcount, 10)),np.argwhere(colcount < np.percentile(colcount, 20)))
    precent30 = np.intersect1d(np.argwhere(colcount >= np.percentile(colcount, 20)),np.argwhere(colcount < np.percentile(colcount, 30)))
    precent40 = np.intersect1d(np.argwhere(colcount >= np.percentile(colcount, 30)),np.argwhere(colcount < np.percentile(colcount, 40)))    
    precent50 = np.intersect1d(np.argwhere(colcount >= np.percentile(colcount, 40)),np.argwhere(colcount < np.percentile(colcount, 50)))
    precent60 = np.intersect1d(np.argwhere(colcount >= np.percentile(colcount, 50)),np.argwhere(colcount < np.percentile(colcount, 60)))
    precent70 = np.intersect1d(np.argwhere(colcount >= np.percentile(colcount, 60)),np.argwhere(colcount < np.percentile(colcount, 70)))
    precent80 = np.intersect1d(np.argwhere(colcount >= np.percentile(colcount, 70)),np.argwhere(colcount < np.percentile(colcount, 80)))
    precent90 = np.intersect1d(np.argwhere(colcount >= np.percentile(colcount, 80)),np.argwhere(colcount < np.percentile(colcount, 90)))
    precent100 = np.argwhere(colcount > np.percentile(colcount, 90))

    total = 0
    for x in histcolcount:
        if x < max(colcount[precent10]):
            total += 1
        if x >= max(colcount[precent10]) and x < max(colcount[precent20]):
            total += 2
        if x >= max(colcount[precent20]) and x < max(colcount[precent30]):
            total += 3
        if x >= max(colcount[precent30]) and x < max(colcount[precent40]):
            total += 4
        if x >= max(colcount[precent40]) and x < max(colcount[precent50]):
            total += 5
        if x >= max(colcount[precent50]) and x < max(colcount[precent60]):
            total += 6
        if x >= max(colcount[precent60]) and x < max(colcount[precent70]):
            total += 7
        if x >= max(colcount[precent70]) and x < max(colcount[precent80]):
            total += 8
        if x >= max(colcount[precent80]) and x < max(colcount[precent90]):
            total += 9
        if x >= max(colcount[precent90]):
            total += 10
    print("\nOn Average compation of windows:")
    print(total/len(histcolcount))

def HistoneJaccard (histrange,names):
    histcolcount = (histrange == 1).sum(axis = 0)
    
    non0histcolcount = []
    non0names = []
    for x in range(len(histcolcount) - 1):
        if(np.sum(histrange[:,x]) > 0):
            non0histcolcount.append(histrange[:,x])
            non0names.append(names[x])
    
    Jaccardindex = []
    Jaccarddistance = []

    for x in range(len(non0histcolcount)):
        Jaccardlineindex = []
        Jaccardlinedistance = []
        for y in range(len(non0histcolcount)):
            test1 = np.argwhere(non0histcolcount[x] == 1)
            test2 = np.argwhere(non0histcolcount[y] == 1)
            
            intersect = np.intersect1d(test1,test2)
            Union = np.union1d(test1,test2)

            Jaccardlineindex.append(round((len(intersect) / len(Union))* 100, 2))

            Jaccardlinedistance.append(1 - (len(intersect) / len(Union)))

        Jaccardindex.append(Jaccardlineindex)
        Jaccarddistance.append(Jaccardlinedistance)
    
    ji = pd.DataFrame(Jaccardindex, index = non0names, columns = non0names)
    ji.to_csv("./files/Jaccardindex.csv")

    plt.figure(figsize=(30,30))
    plt.title("Jaccard Index", size = 50)
    heat_map = sb.heatmap(Jaccardindex, square = True, xticklabels = non0names, yticklabels = non0names, cbar_kws=({"shrink": .80}), cmap = "BuPu")
    plt.savefig("./files/Jaccardindexheat.png")
    plt.clf()

    ji = pd.DataFrame(Jaccarddistance, index = non0names,columns = non0names)
    ji.to_csv("./files/Jaccarddistance.csv")

    plt.figure(figsize=(30,30))
    plt.title("Jaccard Distance", size = 50)
    dheat_map = sb.heatmap(Jaccarddistance, square = True, xticklabels = non0names, yticklabels = non0names, cbar_kws=({"shrink": .80}), cmap = "BuPu")
    plt.savefig("./files/Jaccardistanceheat.png")


def cluster (histrange, names, histrangenums, feature_data):

    histcolcount = (histrange == 1).sum(axis = 0)
    
    clustergroups = []
    non0histcolcount = []
    non0names = []
    Histnums = []

    for x in histrangenums:
        Histnums.append(str(x[0]) + " - " + str(x[1]))

    for x in range(len(histcolcount) - 1):
        if(np.sum(histrange[:,x]) > 20):
            non0histcolcount.append(histrange[:,x])
            non0names.append(names[x])

    trials = 0

    random.seed(datetime.now())
    while trials < 5:
        print("\n\nTest ", trials + 1)

        np1 = random.randint(0, len(non0histcolcount)-1)
        np2 = random.randint(0, len(non0histcolcount)-1)
        np3 = random.randint(0, len(non0histcolcount)-1)
        
        while np1 == np2 or np2 == np3:
            np2 = random.randint(0, len(non0histcolcount)-1)
            np3 = random.randint(0, len(non0histcolcount)-1)

        oldnp1 = 0
        oldnp2 = 0 
        oldnp3 = 0

        i = 0
        while i < 15:
            indext = 0
            cluster1 = []
            cluster2 = []
            cluster3 = []

            j = 0
            for x in non0histcolcount:
            
                test1 = np.argwhere(x == 1)

                test2 = np.argwhere(non0histcolcount[np1] == 1)
                intersect = np.intersect1d(test1,test2)
                Union = np.union1d(test1,test2)
                dist1 = (1 - (len(intersect) / len(Union)))

                test2 = np.argwhere(non0histcolcount[np2] == 1)
                intersect = np.intersect1d(test1,test2)
                Union = np.union1d(test1,test2)
                dist2 = (1 - (len(intersect) / len(Union)))

                test2 = np.argwhere(non0histcolcount[np3] == 1)
                intersect = np.intersect1d(test1,test2)
                Union = np.union1d(test1,test2)
                dist3 = (1 - (len(intersect) / len(Union)))

                if dist1 <= dist2 and dist1 <= dist3:
                    cluster1.append([x,j])
                elif dist2 < dist1 and dist2 < dist3:
                    cluster2.append([x,j])
                elif dist3 < dist1 and dist3 < dist2:
                    cluster3.append([x,j])
                j += 1
            
            clust1dist = []
            for x in cluster1:
                dist = 0
                for y in cluster1:
                    test1 = np.argwhere(x[0] == 1)
                    test2 = np.argwhere(y[0] == 1)

                    intersect = np.intersect1d(test1,test2)
                    Union = np.union1d(test1,test2)
                    dist = (1 - (len(intersect) / len(Union)))
                clust1dist.append(dist)

            indext += np.sum(clust1dist)/ len(clust1dist)
            index = np.argmin(clust1dist)
            np1 = cluster1[index][1]

            clust2dist = []
            for x in cluster2:
                dist = 0
                for y in cluster2:
                    test1 = np.argwhere(x[0] == 1)
                    test2 = np.argwhere(y[0] == 1)

                    intersect = np.intersect1d(test1,test2)
                    Union = np.union1d(test1,test2)
                    dist = (1 - (len(intersect) / len(Union)))
                clust2dist.append(dist)

            indext += np.sum(clust2dist)/ len(clust2dist)
            index = np.argmin(clust2dist)
            np2 = cluster2[index][1]

            clust3dist = []
            for x in cluster3:
                dist = 0
                for y in cluster3:
                    test1 = np.argwhere(x[0] == 1)
                    test2 = np.argwhere(y[0] == 1)

                    intersect = np.intersect1d(test1,test2)
                    Union = np.union1d(test1,test2)
                    dist = (1 - (len(intersect) / len(Union)))
                clust3dist.append(dist)

            indext += np.sum(clust3dist)/ len(clust3dist)
            index = np.argmin(clust3dist)
            np3 = cluster3[index][1]

            i += 1

            if oldnp1 == np1  and oldnp2 == np2 and oldnp3 == np3:
                clustergroups.append([indext, cluster1, cluster2, cluster3])
                break

            oldnp1 = np1
            oldnp2 = np2
            oldnp3 = np3

        trials += 1
        print(indext)

    min = 100
    for x in clustergroups:
        if x[0] < min:
            min = x[0]
            bestnp = x



    cluster1 = bestnp[1]
    cluster2 = bestnp[2]
    cluster3 = bestnp[3]

    num1 = []
    num2 = []
    num3 = []
    names1 = []
    names2 = []
    names3 = []

    for x in cluster1:
        names1.append(non0names[x[1]])
        num1.append(x[0])
    
    for x in cluster2:
        names2.append(non0names[x[1]])
        num2.append(x[0])

    for x in cluster3:
        names3.append(non0names[x[1]])
        num3.append(x[0])


    print("length of cluster 1: ",len(cluster1))
    print("length of cluster 2: ",len(cluster2))
    print("length of cluster 3: ",len(cluster3))

    Histnums
    df1 = pd.DataFrame(num1)
    plt.figure(figsize=(30,30))
    plt.title("Cluster 1", size = 50)
    heat_map = sb.heatmap(df1, square = True, yticklabels = names1, xticklabels = Histnums, cbar_kws=({"shrink": .80}), cmap = "BuPu")
    plt.savefig("./files/clust1.png")
    plt.clf()

    df2 = pd.DataFrame(num2)
    plt.figure(figsize=(30,30))
    plt.title("Cluster 2", size = 50)
    heat_map = sb.heatmap(df2, square = True, yticklabels = names2, xticklabels = Histnums, cbar_kws=({"shrink": .80}), cmap = "BuPu")
    plt.savefig("./files/clust2.png")
    plt.clf()

    df3 = pd.DataFrame(num3)
    plt.figure(figsize=(30,30))
    plt.title("Cluster 3", size = 50)
    heat_map = sb.heatmap(df3, square = True, yticklabels = names3, xticklabels = Histnums, cbar_kws=({"shrink": .80}), cmap = "BuPu")
    plt.savefig("./files/clust3.png")
    plt.clf()

    print("\n")
    feature(cluster1,cluster2,cluster3, feature_data, non0names, Histnums)
    return

def feature(cluster1, cluster2, cluster3, feature_data, non0names, Histnums):

    hist1 = feature_data.iloc[:,13]
    lad = feature_data.iloc[:,9]
    
    histindex = np.argwhere(hist1 == 1)
    ladindex = np.argwhere(lad == 1)

    hist1clust1 = []
    ladclust1 = []
    for x in range(len(cluster1)):
        clustarg1 = np.argwhere(cluster1[x][0] == 1)

        intersect = np.intersect1d(histindex,clustarg1)
        hist1clust1.append(len(intersect)/81)

        intersect = np.intersect1d(ladindex,clustarg1)
        ladclust1.append(len(intersect)/81)

    boxplot = (hist1clust1, ladclust1)
    plt.title("Clust1 Box", size = 50)
    plt.boxplot(boxplot)
    plt.savefig("./files/clust1box.png")

    hist1clust2 = []
    ladclust2 = []
    for x in range(len(cluster2)):
        clustarg2 = np.argwhere(cluster2[x][0] == 1)

        intersect = np.intersect1d(histindex,clustarg2)
        hist1clust2.append(len(intersect)/81)

        intersect = np.intersect1d(ladindex,clustarg2)
        ladclust2.append(len(intersect)/81)

    boxplot = (hist1clust2, ladclust2)
    plt.title("Clust2 Box", size = 50)
    plt.boxplot(boxplot)
    plt.savefig("./files/clust2box.png")

    
    hist1clust3 = []
    ladclust3 = []
    for x in range(len(cluster3)):
        clustarg3 = np.argwhere(cluster3[x][0] == 1)

        intersect = np.intersect1d(histindex,clustarg3)
        hist1clust3.append(len(intersect)/81)

        intersect = np.intersect1d(ladindex,clustarg3)
        ladclust3.append(len(intersect)/81)

    boxplot = (hist1clust3, ladclust3)
    plt.title("Clust3 Box", size = 50)
    plt.boxplot(boxplot)
    plt.savefig("./files/clust3box.png")

    fake1 = [0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,69,70,71,72,73,74,75,76,77,78,79,80,81]
    fake2 = [15,16,17,18,19,20,21,22,23,24,25,26,27,28,56,57,58,59,60,61,62,63,64,65,66,67,68]
    fake3 = [29,30,31,32,33,34,35,36,37,38,39,40,41,42,43,44,45,46,47,48,49,50,51,52,53,54,55]

    print("\n\n---------------------------------------------------------------\nCluster1 top 90%:")
    printclustwindows(cluster1, Histnums, 90, 1)

    print("\nCluster2 top 90%:")
    printclustwindows(cluster2, Histnums, 90, 1)

    print("\nCluster3 top 90%:")
    printclustwindows(cluster3, Histnums, 90, 1)

    print("\n\n---------------------------------------------------------------\nCluster1")
    printclustnp(cluster1, fake1, fake2, fake3, non0names)

    print("\nCluster2")
    printclustnp(cluster2, fake1, fake2, fake3, non0names)

    print("\nCluster3")
    printclustnp(cluster3, fake1, fake2, fake3, non0names)

    c1 = printclustwindows(cluster1, Histnums, 0, 0)
    c2 = printclustwindows(cluster2, Histnums, 0, 0)
    c3 = printclustwindows(cluster3, Histnums, 0, 0)
    print("\n\n---------------------------------------------------------------\nCluster 1 and 2 / 0 Precent filter:")
    printUandI(c1, c2)

    print("\nCluster 1 and 3 / 0 Precent filter:")
    printUandI(c1, c3)

    print("\nCluster 2 and 3 / 0 Precent filter:")
    printUandI(c3, c2)

    print("\nClusters 1, 2 and 3 / 0 Precent filter:")
    test = reduce(np.intersect1d, (c1, c2, c3))
    print("Intersect Length: ", len(test))
    print("Intersect Precent: ", (len(test)/81) * 100)
    
    test = reduce(np.union1d, (c1, c2, c3))
    print("Union Length: ", len(test))
    print("Union Precent: ", (len(test)/81) * 100)

    c1 = printclustwindows(cluster1, Histnums, 50, 0)
    c2 = printclustwindows(cluster2, Histnums, 50, 0)
    c3 = printclustwindows(cluster3, Histnums, 50, 0)
    print("\n\n---------------------------------------------------------------\nCluster 1 and Cluster 2 / 50 Precent filter:")
    printUandI(c1, c2)

    print("\nCluster 1 and Cluster 3 / 50 Precent filter:")
    printUandI(c1, c3)

    print("\nCluster 2 and Cluster 3 / 50 Precent filter:")
    printUandI(c3, c2)

    print("\nClusters 1, 2 and 3 / 50 Precent filter:")
    test = reduce(np.intersect1d, (c1, c2, c3))
    print("Intersect Length: ", len(test))
    print("Intersect Precent: ", (len(test)/81) * 100)

    test = reduce(np.union1d, (c1, c2, c3))
    print("Union Length: ", len(test))
    print("Union Precent: ", (len(test)/81) * 100)

def printclustnp (cluster, fake1, fake2, fake3, non0names):

    data1 = []
    data2 = []
    data3 = []
    for x in range(len(cluster)):
        clustarg3 = np.argwhere(cluster[x][0] == 1)

        intersect1 = np.intersect1d(fake1,clustarg3)
        intersect2 = np.intersect1d(fake2,clustarg3)
        intersect3 = np.intersect1d(fake3,clustarg3)

        data1.append(len(intersect1)/len(clustarg3))
        data2.append(len(intersect2)/len(clustarg3))
        data3.append(len(intersect3)/len(clustarg3))

    print("Ideal 1 match:",np.average(data1))
    print("Ideal 2 match:",np.average(data2))
    print("Ideal 3 match:",np.average(data3),"\n")

    if np.average(data1) > .43:
        print ("Edge defined Cluster\n")
        for x in range(len(cluster)):
            clustarg1 = np.argwhere(cluster[x][0] == 1)
            intersect1 = np.intersect1d(fake1,clustarg1)
            if len(intersect1)/len(clustarg1) > .5:
                print(non0names[cluster[x][1]])

    elif np.average(data3) > .43:
        print ("Middle defined Cluster\n")
        for x in range(len(cluster)):
            clustarg1 = np.argwhere(cluster[x][0] == 1)
            intersect1 = np.intersect1d(fake3,clustarg1)
            if len(intersect1)/len(clustarg1) > .5:
                print(non0names[cluster[x][1]])
    
    else:
        print ("Both defined Cluster")
        for x in range(len(cluster)):
            clustarg1 = np.argwhere(cluster[x][0] == 1)
            intersect1 = np.intersect1d(fake2,clustarg1)
            if len(intersect1)/len(clustarg1) > .4:
                print(non0names[cluster[x][1]])

def printclustwindows(cluster, Histnums, prcent, prnt):

    total = np.zeros(81)
    clustarg = []
    for x in range(len(cluster)):
        clustarg = np.argwhere(cluster[x][0] == 1)
        for y in clustarg:
            total[int(y)] += 1
        
        precent = np.argwhere(total > np.percentile(total, prcent)) 

    if prnt == 1:
        for x in precent:
            print(Histnums[int(x)])
    else:
        return precent

def printUandI(cluster1, cluster2):
    test = np.intersect1d(cluster1, cluster2)
    print("Intersect Length: ", len(test))
    print("Intersect Precent: ", (len(test)/81) * 100)
    test = np.union1d(cluster1, cluster2)
    print("Union Length: ", len(test))
    print("Union Precent: ", (len(test)/81) * 100)