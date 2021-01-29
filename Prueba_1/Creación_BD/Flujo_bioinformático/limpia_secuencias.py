#!/usr/local/bin/anaconda3/bin/python 

#Este archivo contiene las funciones necesarias para que, una vez que se tenga para cada especie del estudio una lista de transcritos con sus correspondientes coordenadas codificantes (CDSs),
#se cree la unión de dichas coordenadas para generar intervalos los cuales serán removidos de la secuencia completa (correspondiente a cada cromosoma y a cada especie).
#Adicionalmente contiene funciones que se usaron para verificar que se extrageran correctamente los segmentos de CDSs de cada transcrito y coincidiera con la secuencia en los archivos .cds, esto considerando si el transcrito está en la dirección 5'->3' o viceversa (3'->5')
#SI SE QUIERE AÑADIR ALGUNA ESPECIE NUEVA AL ESTUDIO, PARECE QUE ESTE ARCHIVO NO NECESITARÍA TENER MODIFICACIONES.
#SI SE QUIERE UN TAMAÑO DE VENTANA FIJO EN 300, OTRO TAMAÑO FIJO O TAMAÑO VARIABLE, SE REQUIERE CAMBIAR EN LA FUNCIÓN union_intervalos().
import numpy as np

def obten_numCrom(crom, a):
    try:
        numCrom=a.index(crom)
    except:
        try:
            crom_aux=int(crom)
            numCrom=a.index(crom_aux)
        except:
            numCrom=-1
    return numCrom

def obten_cds(longitudes, secuencia, ini, fi, numCrom):
    res="ERROR"
    inicio=int(ini) #estas dos líneas las usamos porque las coordenadas vienen en formato String desde el método anterior, buscaGTF
    fin=int(fi)
    if numCrom>=0:
        if inicio<=fin:
            if fin<=longitudes[numCrom]:
                res=secuencia[numCrom][inicio-1]
                for j in range(inicio,fin):
                    c=secuencia[numCrom][j]
                    res+=c
            else:
                print("fin>longitudDelCromosoma")
                print("fin: "+str(fin))
                print("longitudDelCromosoma: "+str(longitudes[numCrom]))
        else:
            print("inicio>fin")
            print("inicio: "+str(inicio))
            print("fin: "+str(fin))
    return res

def compl_sec(sec):
    n=len(sec)
    res=compl_nuc(sec[0])
    for i in range(1,n):
        res+=compl_nuc(sec[i])
    return res

def inv_sec(sec):
    n=len(sec)
    res=sec[n-1]
    for i in range(2,n+1):
        res+=sec[n-i]
    return res

def compl_nuc(char):
    res='E'
    if char=='A':
        res='T'
    elif char=='T':
        res='A'
    elif char=='G':
        res='C'
    elif char=='C':
        res='G'
    return res

def imprime_cds(longitudes, secuencia, sentido, arr, numCrom):
    longi=len(arr)
    res=obten_cds(longitudes, secuencia, arr[0][0],arr[0][1],numCrom)
    for i in range(1,longi):
        res+=obten_cds(longitudes, secuencia, arr[i][0], arr[i][1], numCrom)
    if not sentido:
        stri=compl_sec(res)
        str2=inv_sec(stri)
        res=str2
    return res

#lista_nueva: es una cola que tiene como elementos listas de la forma [inicio, fin]
#lista_unida: es una lista que contendrá listas de la forma [inicio, fin] que serán unión de las contenidas en todas las colas "lista_nueva" que le lleguen
def union_intervalos(lista_nueva, lista_unida):
    lon_u=len(lista_unida)
    lon_n=len(lista_nueva)
    l_mix=[]
    for i in range(lon_n):
        ini_i=int(lista_nueva[i][0]) #en la cola los inicios y fines vienen como Strings
        fin_i=int(lista_nueva[i][1])
        j=0
        banw1=True
        banw2=True
        while j<lon_u and banw1 and banw2:
            #Es +-100 en estas dos primeras condiciones para dejar de considerar espacios entre asteríscos que sean más pequeños que el tamaño mínimo de ventana.
            #Es +-300 en estas dos primeras condiciones para dejar de considerar espacios entre asteríscos que sean más pequeños que el tamaño FIJO de ventana.
            if ini_i>lista_unida[j][1]+300:
                j+=1
            elif fin_i<lista_unida[j][0]-300:
                banw1=False
            else:
                inicio=min(ini_i,lista_unida[j][0])
                fin=max(fin_i,lista_unida[j][1])
                l_mix.append(inicio)
                l_mix.append(fin)
                banw2=False
        if j==lon_u:
            lista_unida.append([ini_i,fin_i])
        elif not banw1:
            lista_unida.insert(j,[ini_i,fin_i])
        else: #not bandw2:
            lista_unida.pop(j)
            lista_unida.insert(j,l_mix)

def remueve_transcrito(longitud,secuencia,l_u):
    res=""
    l=len(l_u)
    if l>0:
        for i in range(l):
            if i==l-1:
                res+='*'
                res+=secuencia[l_u[i-1][1]+1:l_u[i][0]]
                if l_u[i][1]<longitud-1:
                    res+=secuencia[l_u[i][1]+1:longitud]
            elif i>0:
                res+='*'
                res+=secuencia[l_u[i-1][1]+1:l_u[i][0]]
            else: #i==0
                res+=secuencia[:l_u[i][0]]
    else:
        res=secuencia
    return res
        
