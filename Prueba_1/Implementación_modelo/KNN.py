#!/usr/local/bin/anaconda3/bin/python

import sys
import csv
import random
import math
import operator
import heapq

from datetime import datetime

if len(sys.argv)>3:
    baseTrain=int(sys.argv[1])
    ka=int(sys.argv[2])
    baseTest=int(sys.argv[3])
else:
    ka=50
    baseTrain=1
    baseTest=1
    
def loadDataset(filename, tSet=[]):
    with open(filename, 'r') as csvfile:
        lines = csv.reader(csvfile)
        dataset = list(lines)
        lon=len(dataset)
        dataset=dataset[1:lon]
        for x in range(lon-1):
            for y in range(7):
                dataset[x][y] = float(dataset[x][y])
            tSet.append(dataset[x])
 
 
def euclideanDistance(instance1, instance2):
    distance = 0
    for x in range(7):
        distance += pow((instance1[x] - instance2[x]), 2)
    return math.sqrt(distance)
 
def getNeighbors(trainingSet, testInstance, k):
    distances = []
    for x in range(len(trainingSet)):
        dist = euclideanDistance(testInstance, trainingSet[x])
        heapq.heappush(distances,(dist,trainingSet[x]))
    neighbors = []
    for j in range(k):
        elem=heapq.heappop(distances)
        neighbors.append(elem[1][-1])
    return neighbors
 
def getResponse(neighbors):
    classVotes = {}
    for x in range(len(neighbors)):
        response = neighbors[x]
        if response in classVotes:
            classVotes[response] += 1
        else:
            classVotes[response] = 1
    sortedVotes = sorted(classVotes.items(), key=operator.itemgetter(1), reverse=True)
    return sortedVotes

def main():
    # prepare data
    trainingSet=[]
    testSet=[]
    loadDataset('baseTrain'+str(baseTrain)+'_prep.csv', trainingSet)
    loadDataset('baseTest'+str(baseTest)+'_prep.csv', testSet)
    print('Train set: ' + repr(len(trainingSet)))
    print('Test set: ' + repr(len(testSet)))
    # generate predictions
    predictions=[]
    #ka = 50
    print(str(ka))
    for x in range(len(testSet)):
        neighbors = getNeighbors(trainingSet, testSet[x], ka)
        result = getResponse(neighbors)
        predictions.append(result)
        if x%1000==0:
            print('> predicted=' + repr(result[0][0]) + ', actual=' + repr(testSet[x][-1]))
            print(datetime.now())
    with open("resultados"+str(baseTrain)+"_"+str(ka)+"_"+str(baseTest)+".csv","w",newline="") as archi:
        archw=csv.writer(archi)
        lon=len(testSet)
        for i in range(lon):
            list_aux=[testSet[i][-1]]
            lon_aux=len(predictions[i])
            for j in range(lon_aux):
                list_aux.append(predictions[i][j][0])
                list_aux.append(predictions[i][j][1])
            archw.writerow(list_aux)

#main
main()
