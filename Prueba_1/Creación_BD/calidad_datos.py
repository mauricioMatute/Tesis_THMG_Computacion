#!/usr/local/bin/anaconda3/bin/python

#Este archivo recibe las secuencias de ADN de cada cromosoma de cada especie ya "limpia", es decir sin genes (transcritos) ortólogos.
#Con esto crea la base de datos con mediciones de las 9 variables en ventanas de tamaño 300 pares de bases (por ahora).
#SI SE AÑADEN NUEVAS ESPECIES AL ESTUDIO ENTONCES ESTE ARCHIVO SE TIENE QUE MODIFICAR PARA INCLUIRLAS.
import random
import csv
from Bio import SeqIO
import numpy as np

#DESDE AQUÍ SE DEBE(N) AÑADIR LA(S) NUEVA(S) ESPECIE(S)
cromosomas_ciona=14
cromosomas_mosca=7
cromosomas_pez_zebra=25
cromosomas_gusano=6
cromosomas_homo=24
cromosomas_levadura=16
cromosomas_pollo=34 #33 cromosomas excepto el 29, más el W y el Z
cromosomas_raton=21

#defino arreglos con los nombres de los cromosomas por cada especie en el estudio.
a_ciona=[i for i in range(1,cromosomas_ciona+1)]
a_gusano=['I','II','III','IV','V','X']
a_homo=[i for i in range(1,cromosomas_homo-1)]
a_homo.append('X')
a_homo.append('Y')
a_levadura=['I','II','III','IV','V','VI','VII','VIII','IX','X','XI','XII','XIII','XIV','XV','XVI']
a_mosca=['2L','2R','3L','3R',4,'X','Y']
a_pollo=[1,2,3,4,5,6,7,8,9,10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,30,31,32,33,'W','Z']
a_raton=[i for i in range(1,cromosomas_raton-1)]
a_raton.append('X')
a_raton.append('Y')
a_pez=[i for i in range(1,cromosomas_pez_zebra+1)]

arr_as=[a_levadura, a_ciona, a_mosca, a_pez, a_pollo, a_raton, a_homo, a_gusano]

Especies=["sec_Saccharomyces_cerevisiae_crom_","sec_Ciona_intestinalis_crom_","sec_Drosophila_melanogaster_crom_","sec_Danio_rerio_crom_","sec_Gallus_gallus_crom_","sec_Mus_musculus_crom_","sec_Homo_sapiens_crom_","sec_Caenorhabditis_elegans_crom_"]

#FUNCIONES
def selec_crom_esp(esp):
    cromosomas=0
    res="/u/scratch/mauricio/code/secuencias/"
    if esp=="sec_Ciona_intestinalis_crom_":
        cromosomas=cromosomas_ciona
        res2="ciona_intestinalis"
    elif esp=="sec_Drosophila_melanogaster_crom_":
        cromosomas=cromosomas_mosca
        res2="mosca"
    elif esp=="sec_Danio_rerio_crom_":
        cromosomas=cromosomas_pez_zebra
        res2="zebra_fish"
    elif esp=="sec_Caenorhabditis_elegans_crom_":
        cromosomas=cromosomas_gusano
        res2="gusano"
    elif esp=="sec_Homo_sapiens_crom_":
        cromosomas=cromosomas_homo
        res2="homo_sapiens"
    elif esp=="sec_Saccharomyces_cerevisiae_crom_":
        cromosomas=cromosomas_levadura
        res2="levadura"
    elif esp=="sec_Gallus_gallus_crom_":
        cromosomas=cromosomas_pollo
        res2="pollo"
    elif esp=="sec_Mus_musculus_crom_":
        cromosomas=cromosomas_raton
        res2="raton"
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
    with open('info_calidad_'+espe+'.csv','w') as newFile:
        newFileWriter=csv.writer(newFile)
        newFileWriter.writerow(["N","A","G","C","T","*","Cromosoma"])
        for ind in range(crom):
            lista[ind].append(arr_as[indi][ind])
            newFileWriter.writerow(lista[ind])


