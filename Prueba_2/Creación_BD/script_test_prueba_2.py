#!/usr/local/bin/anaconda3/bin/python

#Este archivo realiza la prueba 2, inserta aleatoriamente un segmento de ADN de bacteria (ralstonia picketti) en un segmento más largo de ADN de homo sapiens y recorre el segmento formado
#con ventanas secuenciales de 300pb
#SI SE ANIADEN NUEVAS ESPECIES AL ESTUDIO ENTONCES ESTE ARCHIVO SE TIENE QUE MODIFICAR PARA INCLUIRLAS.
import random
import csv
from Bio import SeqIO
import numpy as np

#DESDE AQUI SE DEBE(N) ANIADIR LA(S) NUEVA(S) ESPECIE(S)
cromosomas_homo=24
cromosomas_acinetobacter=1

#defino arreglos con los nombres de los cromosomas por cada especie en el estudio.
a_homo=[i for i in range(1,cromosomas_homo-1)]
a_homo.append('X')
a_homo.append('Y')
a_acinetobacter=[""]

arr_as=[a_acinetobacter, a_homo]

Especies=["Acinetobacter_baumannii_ncgm_237_gca_000828795.ASM82879v1.dna.chromosome.Chromosome", "Homo_sapiens.GRCh38.dna.chromosome."]

#FUNCIONES
def lee_sec(arch):
    fasta_sequences = SeqIO.parse(arch,'fasta')
    for seq_record in fasta_sequences:
        secuencia=str(seq_record.seq).upper()
        longi=len(seq_record.seq)
    return longi, secuencia

def selec_crom_esp(esp):
    cromosomas=0
    #res="/scr/k61san/mmauricio/secuencias/"
    res=""
    if esp=="Acinetobacter_baumannii_ncgm_237_gca_000828795.ASM82879v1.dna.chromosome.Chromosome":
        cromosomas=cromosomas_acinetobacter
        res+="Acinetobacter_baumannii/"
    elif esp=="Homo_sapiens.GRCh38.dna.chromosome.":
        cromosomas=cromosomas_homo
        res+="../homo_sapiens/"
    else:
        print("ERROR, NO SE DIO LA ESPECIE CORRECTA")
    return cromosomas, res

def obten_sequencia(esp, a):
    cromosomas, cad=selec_crom_esp(esp)
    secuencias=[i for i in range(cromosomas)]
    longitudes=np.arange(cromosomas)
    i=0
    for ind in a:
        arch=cad+esp+str(ind)+".fa"
        longitudes[i], secuencias[i]=lee_sec(arch)
        i+=1
    return longitudes, secuencias, cromosomas

def randNumCrom(longitudes,chrom):
    longTot=0
    acumulado=np.arange(chrom)
    for i in range(chrom):
        longTot+=longitudes[i]
        acumulado[i]=longTot
    nume=random.randint(1,longTot)
    j=True
    ind=0
    while j:
        j=(nume>acumulado[ind])
        ind+=1
    return ind

def calcula_sensores(secuencia, inicio):
    respuesta=[i for i in range(9)]
    window=300
    CGcont=0
    CpGcont=0
    Nww=0
    Nss=0
    Nws=0
    Nsw=0
    Nw=0
    Ns=0
    Nrr=0
    Nyy=0
    Nry=0
    Nyr=0
    Nr=0
    Ny=0
    Nmm=0
    Nkk=0
    Nmk=0
    Nkm=0
    Nm=0
    Nk=0
    LTPcont=0
    HTPcont=0
    VTPcont=0
    ITPcont=0
    cAux='E'
    band=False
    contN=0
    div_window=window-1
    for j in range(inicio,inicio+window):
        c=secuencia[j]
        if c=='*':
            band=True
            break
        elif c=='N':
            contN+=1
            div_window-=1
            if (contN/window)>=(1/6):
                band=True
                break
            if cAux in ['A','G','C','T']:
                div_window-=1
            if j==inicio+window-1:
                div_window+=1
        else:
            if c=='C':
                CGcont+=1
                Ns+=1
                Ny+=1
                Nm+=1
                if cAux=='C':
                    Nss+=1
                    Nyy+=1
                    Nmm+=1
                    LTPcont+=1
                elif cAux=='G':
                    Nss+=1
                    Nry+=1
                    Nkm+=1
                    HTPcont+=1
                elif cAux=='A':
                    Nws+=1
                    Nry+=1
                    Nmm+=1
                    ITPcont+=1
                elif cAux=='T':
                    Nws+=1
                    Nyy+=1
                    Nkm+=1
                    HTPcont+=1
            elif c=='G':
                CGcont+=1
                Ns+=1
                Nr+=1
                Nk+=1
                if cAux=='C':
                    CpGcont+=1
                    Nss+=1
                    Nyr+=1
                    Nmk+=1
                    VTPcont+=1
                elif cAux=='G':
                    Nss+=1
                    Nrr+=1
                    Nkk+=1
                    LTPcont+=1
                elif cAux=='A':
                    Nws+=1
                    Nrr+=1
                    Nmk+=1
                    LTPcont+=1
                elif cAux=='T':
                    Nws+=1
                    Nyr+=1
                    Nkk+=1
                    VTPcont+=1
            elif c=='A':
                Nw+=1
                Nr+=1
                Nm+=1
                if cAux=='C':
                    Nsw+=1
                    Nyr+=1
                    Nmm+=1
                    VTPcont+=1
                elif cAux=='G':
                    Nsw+=1
                    Nrr+=1
                    Nkm+=1
                    HTPcont+=1
                elif cAux=='A':
                    Nww+=1
                    Nrr+=1
                    Nmm+=1
                    LTPcont+=1
                elif cAux=='T':
                    Nww+=1
                    Nyr+=1
                    Nkm+=1
                    VTPcont+=1
            elif c=='T':
                Nw+=1
                Ny+=1
                Nk+=1
                if cAux=='C':
                    Nsw+=1
                    Nyy+=1
                    Nmk+=1
                    LTPcont+=1
                elif cAux=='G':
                    Nsw+=1
                    Nry+=1
                    Nkk+=1
                    ITPcont+=1
                elif cAux=='A':
                    Nww+=1
                    Nry+=1
                    Nmk+=1
                    ITPcont+=1
                elif cAux=='T':
                    Nww+=1
                    Nyy+=1
                    Nkk+=1
                    LTPcont+=1
        cAux=c
    if band:
        respuesta[0]=0
        respuesta[1]=0
        respuesta[2]=0
        respuesta[3]=0
        respuesta[4]=0
        respuesta[5]=0
        respuesta[6]=0
        respuesta[7]=0
        respuesta[8]=0
    else:
        CGcontent=CGcont/(window-contN)
        if CGcont!=0:
            CpGratio=CpGcont/(((CGcont/2)**2)/(window-contN))
        else:
            CpGratio=0
        if Ns!=0 and Nw!=0:
            dws=(Nww*Nss-Nsw*Nws)/(Ns*Nw)
        else:
            dws=0
        if Nr!=0 and Ny!=0:
            dry=(Nrr*Nyy-Nyr*Nry)/(Nr*Ny)
        else:
            dry=0
        if Nm!=0 and Nk!=0:
            dmk=(Nmm*Nkk-Nkm*Nmk)/(Nm*Nk)
        else:
            dmk=0
        LTP=LTPcont/(div_window) 
        HTP=HTPcont/(div_window) 
        VTP=VTPcont/(div_window) 
        ITP=ITPcont/(div_window) 

        respuesta[0]=CGcontent
        respuesta[1]=CpGratio
        respuesta[2]=dws
        respuesta[3]=dry
        respuesta[4]=dmk
        respuesta[5]=LTP
        respuesta[6]=HTP
        respuesta[7]=VTP
        respuesta[8]=ITP

    return respuesta
            

#MAIN----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#0->bacteria
#1->humano

with open('base_test2_prueba_2.csv','w') as newFile:
    newFileWriter=csv.writer(newFile)
    newFileWriter.writerow(['GC_content','CpG_ratio','dws','dry','dmk','LTP','HTP','VTP','ITP'])
    v1, v2, v3=obten_sequencia(Especies[0], arr_as[0]) #longitudes, secuencias, cromosomas
    b1, b2, b3=obten_sequencia(Especies[1], arr_as[1])
    #numCrom0=randNumCrom(v1,v3) #la bacteria no necesita obtenerse un cromosoma al azar porque tiene un solo cromosoma.
    numCrom1=randNumCrom(b1,b3)
    inicio0=random.randint(0,v1[0]-400)
    inicio1=random.randint(0,b1[numCrom1-1]-700)
    sec0=v2[0][inicio0:inicio0+400]
    sec1=b2[numCrom1-1][inicio1:inicio1+700]
    sec=sec1[:350]+sec0+sec1[350:]
    for j in range(801):
        res=calcula_sensores(sec,j)
        newFileWriter.writerow(res)

