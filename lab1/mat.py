import numpy as np
import matplotlib.pyplot as plt

def leer_alineacion(nombre_archivo):
    with open(nombre_archivo, 'r') as f:
        lineas = f.readlines()
    return lineas[0].strip(), lineas[1].strip()

def generar_matriz_puntos(seq1, seq2):
    matriz = np.zeros((len(seq1), len(seq2)))
    for i in range(len(seq1)):
        for j in range(len(seq2)):
            if seq1[i] == seq2[j] and seq1[i] != "-":
                matriz[i, j] = 1
    return matriz

def graficar_matriz_puntos(matriz):
    plt.figure(figsize=(8, 8))
    plt.imshow(matriz, cmap='gray_r', interpolation='nearest')
    plt.xlabel("Secuencia 2")
    plt.ylabel("Secuencia 1")
    plt.title("Matriz de puntos de la alineación")
    plt.show()

# Leer la alineación
doc = "alineacion.txt"
seq1, seq2 = leer_alineacion(doc)

# Generar y graficar la matriz de puntos
matriz = generar_matriz_puntos(seq1, seq2)
graficar_matriz_puntos(matriz)
