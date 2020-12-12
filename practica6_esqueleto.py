#! /usr/bin/python

# Complementos Matemáticos I
# TP Final

import argparse
import matplotlib.pyplot as plt
from math import sqrt
from random import randint


class Par:
    '''
    Representa un par de reales
    '''
    x = 0
    y = 0


class LayoutGraph:

    def __init__(self, grafo, iters, refresh, c1, c2, temp, verbose = False):
        '''
        Parámetros:
        grafo: grafo en formato lista
        iters: cantidad de iteraciones a realizar
        refresh: cada cuántas iteraciones graficar. Si su valor es cero, entonces debe graficarse solo al final.
        c1: constante de repulsión
        c2: constante de atracción
        verbose: si está encendido, activa los comentarios
        '''

        # Guardo el grafo
        self.grafo = grafo

        # Inicializo estado
        self.posiciones = {}
        self.fuerzas    = {}

        # Guardo opciones
        self.iters   = iters
        self.refresh = refresh
        self.c1      = c1
        self.c2      = c2
        self.temp    = temp
        self.verbose = verbose

        # Guardo constantes
        self.H     = 100    # Altura (Height) del frame
        self.W     = 100    # Ancho (Width) del frame
        self.area  = (self.H * 2) * (self.W * 2)
        self.cantV = len(self.grafo[0])
        self.k     = sqrt(self.area / self.cantV)

    def layout(self):
        '''
        Aplica el algoritmo de Fruchtermann-Reingold para obtener (y mostrar)
        un layout
        '''
        self.grafo_pos_generate()

        if self.refresh == 0:
            for i in range(self.iters):
                self.step()

        else:
            plt.ion()

            for i in range(self.iters):
                if i % self.refresh == 0:
                    self.show()
                self.step()

            plt.ioff()

        self.show()
        return

    def step(self):
        '''
        Realiza cada paso de la iteración, para inicializar las fuerzas (atracción,
        repulsión y de gravedad) y actualizarlas, y actualizar las posiciones de
        los vértices
        '''
        self.initialize_forces()
        self.compute_attraction_forces()
        self.compute_repulsion_forces()
        self.compute_gravity_forces()
        self.update_positions()
        return

    def grafo_pos_generate(self):
        '''
        Genera posiciones aleatorias en el plano, para cada vértice del grafo,
        y se almacenan en un diccionario
        '''
        V = self.grafo[0]
        posiciones = {}

        for v in V:
            '''
            Para cada vértice, creamos un objeto Par, que representa sus coordenadas
            en el plano, y generamos su posición aleatoriamente
            '''
            coor = Par()
            coor.x = randint((-self.W)+1, self.W-1)
            coor.y = randint((-self.H)+1, self.H-1)
            posiciones[v] = coor

        self.posiciones = posiciones
        return

    def initialize_forces(self):
        '''
        Inicializa las fuerzas de cada vértice en 0
        '''
        V = self.grafo[0]

        for v in V:
            '''
            Para cada vértice, creamos un objeto Par, que representa las coordenadas
            x e y de la fuerza, inicializadas en 0
            '''
            fuerza = Par()
            fuerza.x = 0
            fuerza.y = 0
            self.fuerzas[v] = fuerza

        return

    def f_a(self, distancia):
        '''
        Defino la función f_a como la función que calcula la fuerza de atracción
        '''
        return distancia ** 2 / (self.k * self.c2)

    def compute_attraction_forces(self):
        '''
        Para cada arista, actualiza las fuerzas de atracción de sus vértices
        '''
        E = self.grafo[1]

        for v1,v2 in E:
            
            coor_v1 = self.posiciones[v1]
            coor_v2 = self.posiciones[v2]

            distancia = sqrt((coor_v1.x - coor_v2.x) ** 2 + (coor_v1.y - coor_v2.y) ** 2)
            if distancia >= 0.05:
                # Se actualizan las fuerzas sólo si la distancia es mayor o igual a 0.05
                mod_fa = self.f_a(distancia)

                fx = mod_fa * (coor_v2.x - coor_v1.x) / distancia
                fy = mod_fa * (coor_v2.y - coor_v1.y) / distancia

                self.fuerzas[v1].x += fx
                self.fuerzas[v1].y += fy
                self.fuerzas[v2].x -= fx
                self.fuerzas[v2].y -= fy

        return

    def f_r(self, distancia):
        '''
        Defino la función f_r como la función que calcula la fuerza de repulsión
        '''
        return (self.k * self.c1) ** 2 / distancia

    def compute_repulsion_forces(self):
        '''
        Para cada vértice, actualiza las fuerzas de repulsión respecto a todos los demás
        vértices del grafo
        '''
        V = self.grafo[0]

        for v1 in V:
            for v2 in V:
                if v1 != v2:
                    # Para cada par de vértices distintos, actualizo las fuerzas de repulsión
                    coor_v1 = self.posiciones[v1]
                    coor_v2 = self.posiciones[v2]

                    distancia = sqrt((coor_v1.x - coor_v2.x) ** 2 + (coor_v1.y - coor_v2.y) ** 2)
                    if distancia >= 0.05:
                        # Se actualizan las fuerzas sólo si la distancia es mayor o igual a 0.05
                        mod_fr = self.f_r(distancia)

                        fx = mod_fr * (coor_v1.x - coor_v2.x) / distancia
                        fy = mod_fr * (coor_v1.y - coor_v2.y) / distancia

                        self.fuerzas[v1].x += fx
                        self.fuerzas[v1].y += fy
                        self.fuerzas[v2].x -= fx
                        self.fuerzas[v2].y -= fy

        return

    def compute_gravity_forces(self):
        '''
        Para cada vértice del grafo, calcula la fuerza de gravedad
        '''
        V = self.grafo[0]
        P = self.posiciones

        for v in V:
            distancia = sqrt(P[v].x ** 2 + P[v].y **2)
            if distancia != 0:
                fx = 0.5 * (0 - P[v].x) / distancia
                fy = 0.5 * (0 - P[v].y) / distancia

                self.fuerzas[v].x += fx
                self.fuerzas[v].y += fy

        return

    def update_positions(self):
        '''
        Actualiza las posiciones de los vértices en función de la temperatura
        '''
        V = self.grafo[0]

        for v in V:
            fx = self.fuerzas[v].x
            fy = self.fuerzas[v].y
            modulo = sqrt(fx ** 2 + fy ** 2)

            '''
            Para cada coordenada, actualizo su posición. Si la temperatura es menor al módulo
            del vector de la fuerza, se suma la fuerza sobre su módulo multiplicada por la
            temperatura. Si no, se suma la fuerza.
            '''
            self.posiciones[v].x += fx / modulo * min(modulo, self.temp)
            self.posiciones[v].y += fy / modulo * min(modulo, self.temp)

        return

    def show(self):
        '''
        Grafica el grafo con las posiciones actuales de sus vértices
        '''
        plt.clf()
        V = self.grafo[0]
        E = self.grafo[1]
        P = self.posiciones

        for v in V:
            plt.scatter(P[v].x, P[v].y, s = 200, c = 1, alpha = 0.5)
            plt.text(P[v].x, P[v].y, v)

        for v1,v2 in E:
            plt.plot([P[v1].x, P[v2].x], [P[v1].y, P[v2].y], linewidth = 2)

        plt.show()
        plt.pause(0.0001)
        return


def lee_grafo_archivo(file_path):
    '''
    Crea un grafo a partir de un archivo
    '''
    V = []
    E = []

    with open(file_path, 'r') as f:
        n = int(f.readline())
        for i in range(n):
            V.append(f.readline().strip())

        for line in f:
            edge = line.strip().split()
            edge = tuple(edge)
            E.append(edge)

    return (V, E)


def main():
    # Definimos los argumentos de línea de comando que aceptamos
    parser = argparse.ArgumentParser()

    # Verbosidad, opcional, False por defecto
    parser.add_argument(
        '-v', '--verbose',
        action = 'store_true',
        help = 'Muestra más información al correr el programa'
    )
    # Cantidad de iteraciones, opcional, 50 por defecto
    parser.add_argument(
        '--iters',
        type = int,
        help = 'Cantidad de iteraciones a efectuar',
        default = 50
    )
    # Temperatura inicial
    parser.add_argument(
        '--temp',
        type = float,
        help = 'Temperatura inicial',
        default = 100.0
    )
    # Refresh, opcional, 1 por defecto
    parser.add_argument(
        '-r', '--refresh',
        type = int,
        help = 'Cada cuántas iteraciones graficar. Si es 0, se grafica sólo al final',
        default = 1
    )
    # C1 (constante de repulsión), opcional, 0.1 por defecto
    parser.add_argument(
        '--c1',
        type = float,
        help = 'Constante de repulsión para modificar el esparcimiento',
        default = 0.1
    )
    # C2 (constante de atracción), opcional, 10.0 por defecto
    parser.add_argument(
        '--c2',
        type = float,
        help = 'Constante de atracción para modificar el esparcimiento',
        default = 10.0
    )
    # Archivo del cual leer el grafo
    parser.add_argument(
        'file_name',
        help = 'Archivo del cual leer el grafo a dibujar'
    )

    args = parser.parse_args()

    # Creamos nuestro objeto LayoutGraph
    layout_gr = LayoutGraph(
        grafo   = lee_grafo_archivo(args.file_name),
        iters   = args.iters,
        refresh = args.refresh,
        c1      = args.c1,
        c2      = args.c2,
        temp    = args.temp,
        verbose = args.verbose
    )

    # Ejecutamos el layout
    layout_gr.layout()
    return


if __name__ == '__main__':
    main()
