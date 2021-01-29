#!/usr/local/bin/anaconda3/bin/python 


#Aquí hay una función que se usa para extraer los transcritos que tienen un ortólogo en alguna especie de las del estudio, esto lo hace del archivo mauricio.Proteinortho-graph
#Adicionalmente hay funciones para que una vez que tiene el identificador del transcrito, los busca en los archivos .cds de cada especie para extraer info adicional como el cromosoma donde se encuentra el transcrito.
#Esto con el fin de hacer una comparación entre las secuencias de los archivos .cds y las obtenidas de las coordenadas codificantes de los archivos .gtf, pero solo se hicieron para comprobar el correcto entendimiento de la información.
#SI SE QUIEREN AÑADIR NUEVAS ESPECIES AL ESTUDIO ENTONCES SE TIENE QUE CREAR UN NUEVO ARCHIVO mauricio.Proteinortho-graph QUE ENCUENTRE TRANSCRITOS ORTÓLOGOS CONSIDERANDO LA NUEVA ESPECIE.
import random
import csv
from Bio import SeqIO
import numpy as np


def lee_sec_gen_cds(arch, trans):
    fasta_sequences = SeqIO.parse(arch,'fasta')
    lon=len(trans)
    conta=0
    for seq_record in fasta_sequences:
        desc=seq_record.description
        n0=desc.find(">")
        n1=desc.find(" ")
        tr=desc[n0+1:n1]
        ind=0
        band=True
        while ind<lon and band:
            if len(trans[ind][1])>0:
                trans_aux=trans[ind][0]+"."+trans[ind][1]
            else:
                trans_aux=trans[ind][0]
            if tr==trans_aux:
                band=False
            ind+=1
        if not band : #tr==trans[ind-1][0] or tr==trans[ind-1][0]+"."+trans[ind-1][1]                                                                                                            
            conta+=1
            n2=desc.find("cds chromosome:")
            cte_n2=len("cds chromosome:") #=15 
            if n2<=0:
                n2=desc.find("cds scaffold:")
                cte_n2=len("cds scaffold:") #=13                 
            n3=desc.find("gene:")
            n4=desc.find(":",n2+cte_n2,n3)
            n5=desc.find(":",n4+1,n3)
            n6=desc.find(":",n5+1,n3)
            n7=desc.find(":",n6+1,n3)
            n8=desc.find(":",n7+1,n3)
            n9=desc.find(" ",n3+5)
            crom=desc[n4+1:n5]
            inicio=desc[n5+1:n6]
            fin=desc[n6+1:n7]
            gen=desc[n3+5:n9]
            sequ=seq_record.seq
            l_sec=len(sequ)
            trans[ind-1].append(crom)
            trans[ind-1].append(inicio)
            trans[ind-1].append(fin)
            trans[ind-1].append(gen)
            trans[ind-1].append(sequ)
            trans[ind-1].append(l_sec)
        if conta==lon:
            break

def llena_arr_especie(arrL, arrC, arrM, arrPe, arrPo, arrR, arrH, arrG, esp, t, ver):
    li=[t,ver]                       
    if esp=="Homo_sapiens.GRCh38.cds.all.fa":
        ele=len(arrH)
        res=encuentra_t(arrH, ele, t)
        if res==ele:
            arrH.append(li)
    elif esp=="Mus_musculus.GRCm38.cds.all.fa":
        ele=len(arrR)
        res=encuentra_t(arrR, ele, t)
        if res==ele:
            arrR.append(li)
    elif esp=="Danio_rerio.GRCz11.cds.all.fa":
        ele=len(arrPe)
        res=encuentra_t(arrPe, ele, t)
        if res==ele:
            arrPe.append(li)
    elif esp=="Caenorhabditis_elegans.WBcel235.cds.all.fa":
        ele=len(arrG)
        res=encuentra_t(arrG, ele, t)
        if res==ele:
            arrG.append(li)
    elif esp=="Drosophila_melanogaster.BDGP6.28.cds.all.fa":
        ele=len(arrM)
        res=encuentra_t(arrM, ele, t)
        if res==ele:
            arrM.append(li)
    elif esp=="Gallus_gallus.GRCg6a.cds.all.fa":
        ele=len(arrPo)
        res=encuentra_t(arrPo, ele, t)
        if res==ele:
            arrPo.append(li)
    elif esp=="Ciona_intestinalis.KH.cds.all.fa":
        ele=len(arrC)
        res=encuentra_t(arrC, ele, t)
        if res==ele:
            arrC.append(li)
    elif esp=="Saccharomyces_cerevisiae.R64-1-1.cds.all.fa":
        ele=len(arrL)
        res=encuentra_t(arrL, ele, t)
        if res==ele:
            arrL.append(li)

def encuentra_t(arre,ele,t):
    #ele=len(arre)   
    m=0
    bande=True
    while m<ele and bande:
        if t==arre[m][0]:
            bande=False
            m-=1
        m+=1
    return m

def proteinortho_graph_cds():
    f=open("/u/maribel/genomes/mauricio.proteinortho-graph","r")
    frl=f.readlines()
    listaLevadura=[]
    listaCiona=[]
    listaMosca=[]
    listaPez=[]
    listaPollo=[]
    listaRaton=[]
    listaHumano=[]
    listaGusano=[]
    l=0
    for x in frl:
        c=x[0]
        if c=="#":
            n0=len(x)
            n1=x.find(" ")
            n2=x.find(chr(9))
            e1_aux=x[n1+1:n2]
            n3=x.find(chr(9),n2+1,n0)
            e2_aux=x[n2+1:n3]
            if e1_aux.find(".cds.all.fa")>=0 and e2_aux.find(".cds.all.fa")>=0:
                e1=e1_aux
                e2=e2_aux
        else:
            n4=x.find(chr(9))
            t1=x[:n4]
            n5=x.find(chr(9),n4+1,len(x))
            t2=x[n4+1:n5]
            transcrito1=t1
            versionT1=""
            npunt1=t1.find(".")
            if e1!="Caenorhabditis_elegans.WBcel235.cds.all.fa" and npunt1>=0:
                transcrito1=t1[:npunt1]
                versionT1=t1[npunt1+1:]
            llena_arr_especie(listaLevadura, listaCiona, listaMosca, listaPez, listaPollo, listaRaton, listaHumano, listaGusano, e1, transcrito1, versionT1)
            transcrito2=t2
            versionT2=""
            npunt2=t2.find(".")
            if e2!="Caenorhabditis_elegans.WBcel235.cds.all.fa" and npunt2>=0:
                transcrito2=t2[:npunt2]
                versionT2=t2[npunt2+1:]
            llena_arr_especie(listaLevadura, listaCiona, listaMosca, listaPez, listaPollo, listaRaton, listaHumano, listaGusano, e2, transcrito2, versionT2)
    f.close()
    contador=0
    listaCDSs=["Saccharomyces_cerevisiae.R64-1-1.cds.all.fa", "Ciona_intestinalis.KH.cds.all.fa", "Drosophila_melanogaster.BDGP6.28.cds.all.fa", "Danio_rerio.GRCz11.cds.all.fa", "Gallus_gallus.GRCg6a.cds.all.fa", "Mus_musculus.GRCm38.cds.all.fa", "Homo_sapiens.GRCh38.cds.all.fa", "Caenorhabditis_elegans.WBcel235.cds.all.fa"]
    listaEspecies=[listaLevadura, listaCiona, listaMosca, listaPez, listaPollo, listaRaton, listaHumano, listaGusano]
    for ind in listaEspecies:
        espe=listaCDSs[contador]
        lee_sec_gen_cds("/u/maribel/genomes/"+espe, ind)
        contador+=1
    return listaLevadura, listaCiona, listaMosca, listaPez, listaPollo, listaRaton, listaHumano, listaGusano

def proteinortho_graph():
    f=open("/u/maribel/genomes/mauricio.proteinortho-graph","r")
    frl=f.readlines()
    listaLevadura=[]
    listaCiona=[]
    listaMosca=[]
    listaPez=[]
    listaPollo=[]
    listaRaton=[]
    listaHumano=[]
    listaGusano=[]
    l=0
    for x in frl:
        c=x[0]
        if c=="#":
            n0=len(x)
            n1=x.find(" ")
            n2=x.find(chr(9))
            e1_aux=x[n1+1:n2]
            n3=x.find(chr(9),n2+1,n0)
            e2_aux=x[n2+1:n3]
            if e1_aux.find(".cds.all.fa")>=0 and e2_aux.find(".cds.all.fa")>=0:
                e1=e1_aux
                e2=e2_aux
        else:
            n4=x.find(chr(9))
            t1=x[:n4]
            n5=x.find(chr(9),n4+1,len(x))
            t2=x[n4+1:n5]
            transcrito1=t1
            versionT1=""
            npunt1=t1.find(".")
            if e1!="Caenorhabditis_elegans.WBcel235.cds.all.fa" and npunt1>=0:
                transcrito1=t1[:npunt1]
                versionT1=t1[npunt1+1:]
            llena_arr_especie(listaLevadura, listaCiona, listaMosca, listaPez, listaPollo, listaRaton, listaHumano, listaGusano, e1, transcrito1, versionT1)
            transcrito2=t2
            versionT2=""
            npunt2=t2.find(".")
            if e2!="Caenorhabditis_elegans.WBcel235.cds.all.fa" and npunt2>=0:
                transcrito2=t2[:npunt2]
                versionT2=t2[npunt2+1:]
            llena_arr_especie(listaLevadura, listaCiona, listaMosca, listaPez, listaPollo, listaRaton, listaHumano, listaGusano, e2, transcrito2, versionT2)
    f.close()
    return listaLevadura, listaCiona, listaMosca, listaPez, listaPollo, listaRaton, listaHumano, listaGusano
