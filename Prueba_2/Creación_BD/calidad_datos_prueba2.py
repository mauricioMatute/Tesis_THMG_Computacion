#!/usr/local/bin/anaconda3/bin/python

#Este archivo recibe las secuencias de ADN de cada cromosoma de cada especie ya "limpia", es decir sin genes (transcritos) ortólogos.
#Con esto crea la base de datos con mediciones de las 9 variables en ventanas de tamaño 300 pares de bases (por ahora).
#SI SE AÑADEN NUEVAS ESPECIES AL ESTUDIO ENTONCES ESTE ARCHIVO SE TIENE QUE MODIFICAR PARA INCLUIRLAS.
import random
import csv
from Bio import SeqIO
import numpy as np

#DESDE AQUÍ SE DEBE(N) AÑADIR LA(S) NUEVA(S) ESPECIE(S)
cromosomas_acinetobacter=1

#defino arreglos con los nombres de los cromosomas por cada especie en el estudio.
a_acinetobacter=[""]

arr_as=[a_acinetobacter]

Especies=["Acinetobacter_baumannii_ncgm_237_gca_000828795.ASM82879v1.dna.chromosome.Chromosome"]

#FUNCIONES
def selec_crom_esp(esp):
    cromosomas=0
    res="/scr/k61san/mmauricio/prueba_2/Acinetobacter_baumannii/"
    if esp=="Acinetobacter_baumannii_ncgm_237_gca_000828795.ASM82879v1.dna.chromosome.Chromosome":
        cromosomas=cromosomas_acinetobacter
        #res2="Acinetobacter_baumannii"
        res2=esp
    else:
        print("ERROR, NO SE DIO LA ESPECIE CORRECTA")
    return cromosomas, res, res2

def obten_info_secuencia(esp, a):
    cromosomas, cad, espe=selec_crom_esp(esp)
    res=[i for i in range(cromosomas)]
    i=0
    for ind in a:
        arch=cad+esp+str(ind)+".fa"
        res[i]=lee_sec(arch)
        i+=1
    return res, cromosomas, espe

def lee_sec(arch):
    fasta_sequences=SeqIO.parse(arch, 'fasta')
    for seq_record in fasta_sequences:
        cad=str(seq_record.seq)
        lon=len(cad)
        N=cad.count("N")
        A=cad.count("A")
        G=cad.count("G")
        C=cad.count("C")
        T=cad.count("T")
        Ast=cad.count("*")
        N=N/lon
        A=A/lon
        G=G/lon
        C=C/lon
        T=T/lon
        Ast=Ast/lon
    return [N,A,G,C,T,Ast,lon]

#MAIN----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#0->levadura 
#1->ciona_intestinalis
#2->mosca
#3->pez_zebra
#4->pollo
#5->raton
#6->humano
#7->gusano
lon=len(arr_as)
for indi in range(lon):
    lista, crom, espe=obten_info_secuencia(Especies[indi], arr_as[indi])
    with open('info_calidad_Abaumannii.csv','w') as newFile:
        newFileWriter=csv.writer(newFile)
        newFileWriter.writerow(["N","A","G","C","T","*","Cromosoma"])
        for ind in range(crom):
            lista[ind].append(arr_as[indi][ind])
            newFileWriter.writerow(lista[ind])


