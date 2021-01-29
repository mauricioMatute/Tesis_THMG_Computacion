#!/usr/local/bin/anaconda3/bin/python

#Este archivo recibe las secuencias de ADN de cada cromosoma de cada especie ya "limpia", es decir sin genes (transcritos) ortologos.
#Con esto crea la base de datos con mediciones de las 9 variables en ventanas de tamanio 300 pares de bases (por ahora).
#SI SE ANIADEN NUEVAS ESPECIES AL ESTUDIO ENTONCES ESTE ARCHIVO SE TIENE QUE MODIFICAR PARA INCLUIRLAS.
import random
import csv
from Bio import SeqIO
import numpy as np

#DESDE AQUI SE DEBE(N) ANIADIR LA(S) NUEVA(S) ESPECIE(S)
cromosomas_ciona=14
cromosomas_mosca=7
cromosomas_pez_zebra=25
cromosomas_gusano=6
cromosomas_homo=24
cromosomas_levadura=16
cromosomas_pollo=34 #33 cromosomas excepto el 29, mas el W y el Z
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
#arr_ite=[41000,250000,250000,450000,350000,850000,1000000,150000]
#arr_ite=[1000,2500,2500,4500,3500,8500,10000,2500]
arr_ite=[250,625,625,1125,875,2125,2500,625]

Especies=["sec_Saccharomyces_cerevisiae_crom_","sec_Ciona_intestinalis_crom_","sec_Drosophila_melanogaster_crom_","sec_Danio_rerio_crom_","sec_Gallus_gallus_crom_","sec_Mus_musculus_crom_","sec_Homo_sapiens_crom_","sec_Caenorhabditis_elegans_crom_"]

#FUNCIONES
def lee_sec(arch):
    fasta_sequences = SeqIO.parse(arch,'fasta')
    for seq_record in fasta_sequences:
        secuencia=str(seq_record.seq).upper()
        longi=len(seq_record.seq)
    return longi, secuencia

def selec_crom_esp(esp):
    cromosomas=0
    res="/scr/k61san/mmauricio/secuencias/"
    if esp=="sec_Ciona_intestinalis_crom_":
        cromosomas=cromosomas_ciona
    elif esp=="sec_Drosophila_melanogaster_crom_":
        cromosomas=cromosomas_mosca
    elif esp=="sec_Danio_rerio_crom_":
        cromosomas=cromosomas_pez_zebra
    elif esp=="sec_Caenorhabditis_elegans_crom_":
        cromosomas=cromosomas_gusano
    elif esp=="sec_Homo_sapiens_crom_":
        cromosomas=cromosomas_homo
    elif esp=="sec_Saccharomyces_cerevisiae_crom_":
        cromosomas=cromosomas_levadura
    elif esp=="sec_Gallus_gallus_crom_":
        cromosomas=cromosomas_pollo
    elif esp=="sec_Mus_musculus_crom_":
        cromosomas=cromosomas_raton
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

def calcula_sensores(longitudes, secuencia, a, cromosomas, ind_esp):
    respuesta=[i for i in range(12)]
    window=300  #random.randint(100,1000)
    numCrom=randNumCrom(longitudes,cromosomas)
    inicio=random.randint(0,longitudes[numCrom-1]-window)
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
        c=secuencia[numCrom-1][j]
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
        respuesta[9]=a[numCrom-1]
        respuesta[10]=inicio
        respuesta[11]=ind_esp

    return respuesta
            

#MAIN----------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------

#0->levadura 
#1->ciona_intestinalis
#2->mosca
#3->pez_cebra
#4->pollo
#5->raton
#6->humano
#7->gusano
lon=len(arr_as)
with open('baseTest4.csv','w') as newFile:
    newFileWriter=csv.writer(newFile)
    newFileWriter.writerow(['GC_content','CpG_ratio','dws','dry','dmk','LTP','HTP','VTP','ITP','cromosoma','inicio_ventana','especie'])
    for indi in range(lon):
        v1, v2, v3=obten_sequencia(Especies[indi], arr_as[indi])
        iteraciones=arr_ite[indi]
        for k in range(iteraciones):
            bande=True
            while bande:
                res=calcula_sensores(v1, v2, arr_as[indi], v3,indi)
                bande=(res[0]==0 and res[1]==0 and res[2]==0 and res[3]==0 and res[4]==0 and res[5]==0 and res[6]==0 and res[7]==0 and res[8]==0)
            newFileWriter.writerow(res)


