#!/usr/local/bin/anaconda3/bin/python

import sys
import csv
import random
import math
import operator
import heapq

from datetime import datetime
from sklearn import metrics

def classification_report_csv(report):
    report_data = []
    lines = report.split('\n')
    for line in lines[2:-3]:
        row = []
        row_data = line.split('      ')
        report_data.append(row_data)
    return report_data

def main():
    dataset_completo=[]
    for baseTrain in range(1,3):
        for ka in [10,20,50,80,100,200,300,400,500]:
            res=[]
            lista=[]
            for baseTest in range(1,5):
                with open("resultados"+str(baseTrain)+"_"+str(ka)+"_"+str(baseTest)+".csv","r") as csvfile:
                    lines = csv.reader(csvfile)
                    dataset = list(lines)
                    lon=len(dataset)
                    for x in range(lon):
                        res.append(dataset[x][0])
                        lon_aux=len(dataset[x])
                        lista_aux=[]
                        for y in range(2,lon_aux,2):
                            lista_aux.append([dataset[x][y-1],int(dataset[x][y])])
                        lista.append(lista_aux)
                        ds_compl_aux=[baseTrain,ka]
                        for i in range(lon_aux):
                            ds_compl_aux.append(dataset[x][i])
                        dataset_completo.append(ds_compl_aux)
            longi=len(lista)
            lista2=[]
            for ind in range(longi):
                lista2.append(lista[ind][0][0])
            confusion_mat=metrics.confusion_matrix(res,lista2,labels=['0','1','2','3','4','5','6','7'])
            report=metrics.classification_report(res,lista2,labels=['0','1','2','3','4','5','6','7'])
            lista_report=classification_report_csv(report)
            with open("matConfusion_"+str(baseTrain)+"_"+str(ka)+".csv","w",newline="") as archiv:
                archivw=csv.writer(archiv)
                lon1=len(confusion_mat)
                lon2=len(lista_report)
                for ind1 in range(lon1):
                    archivw.writerow(confusion_mat[ind1])
                for ind2 in range(lon2):
                    archivw.writerow(lista_report[ind2])
    with open("resultados_todos.csv", "w", newline="") as archivo:
        archivow=csv.writer(archivo)
        longitud=len(dataset_completo)
        for indice in range(longitud):
            archivow.writerow(dataset_completo[indice])

main()
