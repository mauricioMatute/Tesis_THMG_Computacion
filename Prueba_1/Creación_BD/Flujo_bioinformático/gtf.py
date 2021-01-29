#!/usr/local/bin/anaconda3/bin/python 

#Este Archivo tiene funciones que se usan para obtener, para cada transcrito, las coordenadas que tienen sus CDSs. Estas coordenadas las obtiene del archivo GTF para cada especie.
#Recibe las listas por especie (listaEspecie) usando el archivo proteinortho_CDS.py. 
#En estas listas vienen los transcritos que tienen un ortólogo en otra especie estudiada para poder mapear las coordenadas de sus CDSs en la secuencia de ADN.
#SI SE QUIEREN AÑADIR NUEVAS ESPECIES AL ESTUDIO, EN PRINCIPIO NO SE NECESITARÍA MODIFICAR ESTE ARCHIVO, SIN EMBARGO, ES POSIBLE QUE ALGO EN LA PROGRAMACIÓN SE HAYA HECHO TOMANDO EN CUENTA EL FORMATO ESPECÍFICO EN LOS NOMBRES O INFORMACIÓN DE LOS ARCHIVOS DE CADA ESPECIE QUE PERTENECE AL ESTUDIO AHORITA, O QUE EL FORMATO DE LA ESPECIE A AÑADIR REQUIERA ALGUNA CONSIDERACIÓN PARTICULAR DE ESTE TIPO, COMO EJEMPLO DEJO EL CASO DEL GUSANO, CUYOS ID DE LOS TRANSCRITOS LLEVAN PUNTOS INTERMEDIOS Y NO LLEVAN LA VERSIÓN, MIENTRAS QUE PARA LAS DEMÁS ESPECIES SÍ VIENE LA VERSIÓN Y SE EXPRESA CON UN PUNTO AL FINAL DEL ID Y LUEGO EL NÚMERO DE VERSIÓN.
import random
import csv
from Bio import SeqIO
import numpy as np
from collections import deque

def cad_entre_espacios(n, n0, x):
    num1=x.find(chr(32),n,n0)
    num2=x.find(chr(9),n,n0)
    if num1<num2 and num1>=0:
        ini=num1
        while x[ini]==chr(32):
            ini+=1
        fini=x.find(chr(32),ini, n0)
        cad=x[ini:fini]
    elif num2<num1 and num2>=0:
        ini=num2
        while x[ini]==chr(9):
            ini+=1
        fini=x.find(chr(9),ini, n0)
        cad=x[ini:fini]
    else:
        print("OH, OH, ALGO NO ESTA BIEN")
        print("num1: "+str(num1))
        print("num2: "+str(num2))
        print("n: "+str(n))
        print("readline: "+x)
    if fini<0: #habrá error si no inicializó fini en alguna de las condiciones if anteriores
        print("OH, OH, ALGO NO ESTA BIEN 2")
        print("num1: "+str(num1))
        print("num2: "+str(num2))
        print("n: "+str(n))
        print("readline: "+x)
        print("fini: "+str(fini))
    return cad, fini  #cad es la cadena que nos interesa y fini es el primer espacio después de la cadena

def buscaGTF_diferentes(strGTF, listaEspecie):
    f=open("/u/scratch/mauricio/gtf/"+strGTF, "r")
    frl=f.readlines()
    band=False
    sentido=True #True='+', False='-'
    longi=len(listaEspecie)
    cola=deque()
    for x in frl:
        n=x.find("transcript")
        if n>=0 and x[n+10]!="_" and n<110: #debería jalar con n>0 porque "transcript" no va a empezar en el primer caractér
            if band: #en este caso no hubo "stop_codon" antes del siguiente transcript
                #print("NO STOP_CODON")
                band=False
                lcola=len(cola)
                suma=0
                for o in range(lcola):
                    numini=int(cola[o][0])
                    numfin=int(cola[o][1])
                    suma+=(numfin-numini+1)
                if listaEspecie[cont][7]!=suma:
                    listaEspecie[cont].append(suma)
                    listaEspecie[cont].append(cola)
                    listaEspecie[cont].append(sentido)
                else: #listaEspecie[cont][7]==suma
                    cosa=listaEspecie.pop(cont)
                    longi-=1
                cola=deque() #reinicializamos la cola para las coordenadas de los cds's del siguiente transcrito
            cont=0
            band2=True
            n0=len(x)
            n1=x.find("transcript_id ",n+10,n0)
            n2=x.find(";",n1+15,n0)
            tran=x[n1+15:n2-1]
            while cont<longi and band2:
                if listaEspecie[cont][0]==tran:
                    band2=False
                    cont-=1
                cont+=1
            if not band2:
                band=True
                nsign=x.find(".",n,n0)
                char, nalgo=cad_entre_espacios(nsign,n0,x)
                if char=='-':
                    sentido=False
                elif char=='+':
                    sentido=True
                else:
                    print("error en el sentido")
                    print("nsign: "+str(nsign))
                    print("readline: "+x)
                    print("char: "+char)
                    print("nalgo: "+str(nalgo))
        elif band:
            #analiza cds y stop_codon
            n=x.find("stop_codon")
            if n>=0 and n<110:
                band=False #siempre debería haber un "stop_codon" antes del siguiente "transcript"
                n0=len(x)
                inicio, n_separa=cad_entre_espacios(n, n0, x) #inicio es la cadena que nos interesa y n_separa es el lugar inmediato después de la cadena
                fin, nalgo=cad_entre_espacios(n_separa, n0, x)
                infi=[inicio,fin]
                if sentido:
                    cola.append(infi)
                else:
                    cola.appendleft(infi)
                lcola=len(cola)
                suma=0
                for o in range(lcola):
                    numini=int(cola[o][0])
                    numfin=int(cola[o][1])
                    suma+=(numfin-numini+1)
                if listaEspecie[cont][7]!=suma:
                    listaEspecie[cont].append(suma)
                    listaEspecie[cont].append(cola)
                    listaEspecie[cont].append(sentido)
                else: #listaEspecie[cont][7]==suma
                    cosa=listaEspecie.pop(cont)
                    longi-=1
                cola=deque() #reinicializamos la cola para las coordenadas de los cds's del siguiente transcrito
            else:
                #analiza cds
                n=x.find("CDS")
                if n>=0 and n<110: #n<110 es porque se llegó a encontrar una subcadena que contenía "CDS" en un lugar 369 de x, obviamente no era el "CDS" que estábamos buscando
                    n0=len(x)
                    inicio, n_separa=cad_entre_espacios(n, n0, x) #inicio es la cadena que nos interesa y n_separa es el lugar inmediato después de la cadena
                    fin, nalgo=cad_entre_espacios(n_separa, n0, x)
                    infi=[inicio,fin]
                    if sentido:
                        cola.append(infi)
                    else:
                        cola.appendleft(infi)
    f.close()

def buscaGTF(strGTF, listaEspecie):
    f=open("/u/scratch/mauricio/gtf/"+strGTF, "r")
    frl=f.readlines()
    band=False
    sentido=True #True='+', False='-'
    longi=len(listaEspecie)
    cola=deque()
    for x in frl:
        n=x.find("transcript")
        if n>=0 and x[n+10]!="_" and n<110: #debería jalar con n>0 porque "transcript" no va a empezar en el primer caractér
            if band: #en este caso no hubo "stop_codon" antes del siguiente transcript
                band=False
                nume1=x.find(chr(32))
                nume2=x.find(chr(9))
                if nume1>0 and nume2>0:
                    nume=min(nume1,nume2)
                else:
                    nume=max(nume1,nume2)
                crom=x[:nume]
                listaEspecie[cont].append(crom)
                listaEspecie[cont].append(cola)
                cola=deque() #reinicializamos la cola para las coordenadas de los cds's del siguiente transcrito
            cont=0
            band2=True
            n0=len(x)
            n1=x.find("transcript_id ",n+10,n0)
            n2=x.find(";",n1+15,n0)
            tran=x[n1+15:n2-1]
            while cont<longi and band2:
                if listaEspecie[cont][0]==tran:
                    band2=False
                    cont-=1
                cont+=1
            if not band2:
                band=True
                nsign=x.find(".",n,n0)
                char, nalgo=cad_entre_espacios(nsign,n0,x)
                if char=='-':
                    sentido=False
                elif char=='+':
                    sentido=True
                else:
                    print("error en el sentido")
                    print("nsign: "+str(nsign))
                    print("readline: "+x)
                    print("char: "+char)
                    print("nalgo: "+str(nalgo))
        elif band:
            #analiza stop_codon
            n=x.find("stop_codon")
            if n>=0 and n<110:
                band=False #siempre debería haber un "stop_codon" antes del siguiente "transcript"
                n0=len(x)
                inicio, n_separa=cad_entre_espacios(n, n0, x) #inicio es la cadena que nos interesa y n_separa es el lugar inmediato después de la cadena
                fin, nalgo=cad_entre_espacios(n_separa, n0, x)
                infi=[inicio,fin]
                if sentido:
                    cola.append(infi)
                else:
                    cola.appendleft(infi)
                nume1=x.find(chr(32))
                nume2=x.find(chr(9))
                if nume1>0 and nume2>0:
                    nume=min(nume1,nume2)
                else:
                    nume=max(nume1,nume2)
                crom=x[:nume]
                listaEspecie[cont].append(crom)
                listaEspecie[cont].append(cola)
                cola=deque() #reinicializamos la cola para las coordenadas de los cds's del siguiente transcrito
            else:
                #analiza cds
                n=x.find("CDS")
                if n>=0 and n<110: #n<110 es porque se llegó a encontrar una subcadena que contenía "CDS" en un lugar 369 de x, obviamente no era el "CDS" que estábamos buscando
                    n0=len(x)
                    inicio, n_separa=cad_entre_espacios(n, n0, x) #inicio es la cadena que nos interesa y n_separa es el lugar inmediato después de la cadena
                    fin, nalgo=cad_entre_espacios(n_separa, n0, x)
                    infi=[inicio,fin]
                    if sentido:
                        cola.append(infi)
                    else:
                        cola.appendleft(infi)
    f.close()
