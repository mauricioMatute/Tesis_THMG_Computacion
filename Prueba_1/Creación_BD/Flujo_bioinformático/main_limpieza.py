#!/usr/local/bin/anaconda3/bin/python

#Este archivo recibe la base donde vienen los transcritos con las coordenadas de sus CDSs dentro de la secuencia del cromosoma correspondiente a cada transcrito.
#Esta base se forma usando los scripts proteinortho_CDS.py & gtf.py
#Con esta base se mapean las coordenadas de todos los CDSs de todos los transcritos para eliminarlas de la secuencia de los cromosomas correspondientes. 
#Esto con ayuda de las funciones en limpia_secuencias.py
#Esto se hace con el objetivo de tener secuencias de ADN que no tengan genes ortólogos entre las especies del estudio, para poder clasificarlos de mejor manera.
#Al final lo que se busca es clasificar material genético en general, y prestando particular atención al material genético no codificante, por tanto retirar genes no es tan relevante en el estudio.
#SI SE AÑADEN NUEVAS ESPECIES AL ESTUDIO ENTONCES ESTE ARCHIVO SE TIENE QUE MODIFICAR PARA INCLUIRLAS.
import random
import csv
from Bio import SeqIO
#from Bio.Alphabet import generic_alphabet   #DEPRECATED
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
import numpy as np
from collections import deque
import proteinortho_CDS as prot
import gtf as gtf
import limpia_secuencias as limp

cromosomas_ciona=14
cromosomas_mosca=7
cromosomas_pez_zebra=25
cromosomas_gusano=6
cromosomas_homo=24
cromosomas_levadura=16
cromosomas_pollo=34 #33 cromosomas excepto el 29, más el W y el Z
cromosomas_raton=21
def lee_sec_desc(arch):
    fasta_sequences = SeqIO.parse(arch,'fasta')
    for seq_record in fasta_sequences:
        secuencia=str(seq_record.seq).upper()
        longi=len(seq_record.seq)
        desc=seq_record.description
    return longi, secuencia, desc

def selec_crom_esp(esp):
    res="/u/scratch/mauricio/"
    if esp=="Ciona_intestinalis.KH.dna.chromosome.":
        cromosomas=cromosomas_ciona
        res+="ciona_intestinalis/"
    elif esp=="Drosophila_melanogaster.BDGP6.28.dna.chromosome.":
        cromosomas=cromosomas_mosca
        res+="mosca/"
    elif esp=="Danio_rerio.GRCz11.dna.chromosome.":
        cromosomas=cromosomas_pez_zebra
        res+="zebra_fish/"
    elif esp=="Caenorhabditis_elegans.WBcel235.dna.chromosome.":
        cromosomas=cromosomas_gusano
        res+="gusano/"
    elif esp=="Homo_sapiens.GRCh38.dna.chromosome.":
        cromosomas=cromosomas_homo
        res+="homo_sapiens/"
    elif esp=="Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.":
        cromosomas=cromosomas_levadura
        res+="levadura/"
    elif esp=="Gallus_gallus.GRCg6a.dna.chromosome.":
        cromosomas=cromosomas_pollo
        res+="pollo/"
    elif esp=="Mus_musculus.GRCm38.dna.chromosome.":
        cromosomas=cromosomas_raton
        res+="raton/"
    else:
        print("ERROR, NO SE DIO LA ESPECIE CORRECTA")
    return cromosomas, res

def obten_sequencia_desc(esp, a):
    cromosomas, cad=selec_crom_esp(esp)
    secuencias=[i for i in range(cromosomas)]
    longitudes=np.arange(cromosomas)
    descripciones=[i for i in range(cromosomas)]
    i=0
    for ind in a:
        arch=cad+esp+str(ind)+".fa"
        longitudes[i], secuencias[i], descripciones[i]=lee_sec_desc(arch)
        i+=1
    return longitudes, secuencias, descripciones

#MAIN----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#defino los archivos .gtf de cada especie en el estudio.
gtfs=["Saccharomyces_cerevisiae.R64-1-1.99.gtf","Ciona_intestinalis.KH.99.gtf","Drosophila_melanogaster.BDGP6.28.99.gtf","Danio_rerio.GRCz11.99.chr_patch_hapl_scaff.gtf","Gallus_gallus.GRCg6a.99.gtf","Mus_musculus.GRCm38.99.chr_patch_hapl_scaff.gtf","Homo_sapiens.GRCh38.99.chr_patch_hapl_scaff.gtf","Caenorhabditis_elegans.WBcel235.99.gtf"]

#Extraigo los transcritos de cada especie contenidos en el archivo mauricio.Proteinortho-graph y les pego la información adicional que viene en los archivos .cds de cada especie,
# generando una lista de transcritos con su información adicional por cada especie.
listaLevadura, listaCiona, listaMosca, listaPez, listaPollo, listaRaton, listaHumano, listaGusano=prot.proteinortho_graph()
contador=0
listaEspecies=[listaLevadura, listaCiona, listaMosca, listaPez, listaPollo, listaRaton, listaHumano, listaGusano] 

#A las listas generadas de cada especie se le adjuntan las coordenadas de las regiones codificantes (CDSs) de cada transcrito, sacadas de los archivos .gtf de cada especie.
for ind in gtfs:
    gtf.buscaGTF(ind, listaEspecies[contador])
    contador+=1

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

#defino, para cada especie en el estudio, el prefijo para los archivos que contienen los cromosomas de cada especie.
Especies=["Saccharomyces_cerevisiae.R64-1-1.dna.chromosome.","Ciona_intestinalis.KH.dna.chromosome.","Drosophila_melanogaster.BDGP6.28.dna.chromosome.","Danio_rerio.GRCz11.dna.chromosome.","Gallus_gallus.GRCg6a.dna.chromosome.","Mus_musculus.GRCm38.dna.chromosome.","Homo_sapiens.GRCh38.dna.chromosome.","Caenorhabditis_elegans.WBcel235.dna.chromosome."]

#con las listas, por cada especie del estudio, que contienen los transcritos, el cromosoma en el que se encuentran y las coordenadas de sus CDSs, se procede a extraerlos de la secuencia correspondiente.
c=0
for esp in Especies:
    longitudes, secuencias, desc=obten_sequencia_desc(esp, arr_as[c])
    lon=len(listaEspecies[c])
    n=esp.find(".")
    l_u=[i for i in range(len(arr_as[c]))] #lista de los intervalos unidos de cada cromosoma
    for j in range(len(l_u)):
        l_u[j]=[]
    for indice in range(lon):
        #listaEspecie[c][indice]:= t, version, crom, coord_cds_gtf
        crom=listaEspecies[c][indice][2]
        numCrom=limp.obten_numCrom(crom, arr_as[c])
        if numCrom>=0:
            lista=listaEspecies[c][indice][3]
            limp.union_intervalos(lista, l_u[numCrom])
    for ind in range(len(l_u)):
        descript=desc[ind]+"_REDUCED"
        sec_limpia=limp.remueve_transcrito(longitudes[ind],secuencias[ind],l_u[ind])
        simple_seq=Seq(sec_limpia)
        #simple_seq.alphabet=generic_alphabet      #DEPRECATED
        simple_seq_r=SeqRecord(simple_seq)
        simple_seq_r.description=descript
        list_seq=[simple_seq_r]
        SeqIO.write(list_seq,"sec_"+esp[:n]+"_crom_"+str(arr_as[c][ind])+".fa","fasta")
    c+=1
