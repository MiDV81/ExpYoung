'''
-------------------------------
Ondas EM y Experimento de Young
-------------------------------
Miguel Durán Vera
'''
# Este código representa varias configuraciones en una misma plot pero distintas subplots

import numpy as np
import matplotlib.pyplot as plt


'''
El objetivo de esta práctica es simular el experimento de Young.
'''


'''
Datos del sistema
'''

c = 3e8  # velocidad de la luz

size = 401  # número de puntos

dx = 2  # Espaciado en el mallado de posición
dy = dx
dt = dx*1e-9/(2*c)  # Espaciado temporal
xpos = np.arange(0, size*dx, dx)*1e-3  # en micras
ypos = np.arange(0, size*dy, dy)*1e-3  # en micras
yy, xx = np.meshgrid(ypos, xpos)  # mallado de posiciones

Ez = np.zeros((size, size))  # Iniciamos array vacíos para las ondas
Hx = np.zeros((size, size))
Hy = np.zeros((size, size))

sumI = np.zeros(size)  # array vacío donde iremos sumando las intensidades

'''
Datos de la simulación
'''

nsim = 1000  # número de simulaciones/pasos temporales
nrep = 10  # número de simulaciones/pasos entre representaciones
trep = 0.01  # tiempo entre representaciones


'''
Datos rendija
'''
#valores de los distintos parametros que estudiaremos
n_s = [1,2,3,4,5,6]  # numero de rendijas
a_s = [20]  # distancia entre rendijas (en nº de veces dy)
d_s = [100]  # distancia a la pantalla de las rendijas (en nº de veces dx)

L = 2  # anchura de las rendijas (en nº de veces dy)
namejpg = f'Variasrendijasv3.jpg'  # nombre que pondremos a la imagen resultante

posx = 200  # posicion donde se coloca la pared con la rendijas

'''
Pulso Inicial
'''

E0 = 1
dxp = 40  # Anchura del pulso
dtp = dxp*1e-9/c  # desviación típica del pulso
t0p = 5*dtp  # Tiempo inicial del pulso


'''
Ecs. propagación
'''

'''
Ey(k) = Ey(k)-1/2*(Hz(k)-Hz(k-1))
Hz(k) = Hz(k)-1/2*(Ey(k+1)-Ey(k))
'''

'''
Datos Medios (Por ahora solo puesto en el vacío)
'''

frontx = 100  # frontera entre ambos medios
fronty = 200

eps0 = 8.854e-12  # F/m
eps1 = 1
eps2 = 1  # Hecho para tener 4 zonas pero ahora solo dos: izqda y derecha
eps3 = 1
eps4 = 1

coef = np.ones((size, size))  # Array con la permitividad del espacio
coef[:frontx, :fronty] = eps2
coef[:frontx, fronty:] = eps1
coef[frontx:, :fronty] = eps3
coef[frontx:, fronty:] = eps4

sigma1 = 0
sigma2 = 0  # Hecho para tener 4 zonas pero ahora solo dos: izqda y derecha
sigma3 = 0
sigma4 = 0

# Array con la conductividad# Array con la permitividad del espacio
sigma = np.ones((size, size))
sigma[:frontx, :fronty] = sigma2
sigma[:frontx, fronty:] = sigma1
sigma[frontx:, :fronty] = sigma3
sigma[frontx:, fronty:] = sigma4

ctc = sigma*dt/(2*eps0*coef)  # Coeficiente tiempo conductancia

Z0 = 1/(c*eps0)
# Función para el calculo de intensidad
def I(E, Z=Z0, n=1): return n/(2*Z)*E**2

'''
Simulación
'''
Ezant = np.zeros(
    (4, 2, size))  # Sin fronteras ajustadas a los medios dieléctricos


def Simulacion(n_list, d_list, a_list):
    fig2 = plt.figure(figsize=(10, 8))
    fig2.suptitle(f'Comparacion entre distintas configuraciones')
    plt.subplots_adjust(wspace=0.2,
                        hspace=0.3) # modificamos el espacio entre subplots
    splot = 231 # indice que determinará la posición de la subplot


    for d in d_list:
        for a in a_list:
            for numrendijas in n_list:
                # Iniciamos array vacíos para las ondas
                Ez = np.zeros((size, size))

                Hx = np.zeros((size, size))
                Hy = np.zeros((size, size))

                sumI = np.zeros(size)
                start = int((size-(numrendijas-1)*a-numrendijas*L)/2)
                for sim in range(0, nsim):

                    t = sim*dt  # Actualizamos el tiempo en el que

                    Ez[1:, 1:] = Ez[1:, 1:]*(1-ctc[1:, 1:])/(1+ctc[1:, 1:]) + 1/(2*coef[1:, 1:]*(1+ctc[1:, 1:]))*(
                        (Hy[1:, 1:]-Hy[:-1, 1:])-(Hx[1:, 1:]-Hx[1:, :-1]))  # Calculamos el nuevo campo eléctrico

                    Ez[:, -1] = Ezant[0, 0, :]  # Arriba
                    Ez[:, 0] = Ezant[1, 0, :]  # Abajo
                    Ez[0, :] = Ezant[2, 0, :]  # Izqda
                    Ez[-1, :] = Ezant[3, 0, :]  # Drch

                    # Movemos las copias anteriores una posición
                    Ezant[0, 0, :] = Ezant[0, 1, :]
                    Ezant[1, 0, :] = Ezant[1, 1, :]
                    Ezant[2, 0, :] = Ezant[2, 1, :]
                    Ezant[3, 0, :] = Ezant[3, 1, :]

                    # Cogemos los datos actuales en la ultima copia
                    Ezant[0, 1, :] = Ez[:, -2]
                    Ezant[1, 1, :] = Ez[:, 1]
                    Ezant[2, 1, :] = Ez[1, :]
                    Ezant[3, 1, :] = Ez[-2, :]

                    # Introducimos el pulso inicial
                    Ez[10, 50:350] = E0*np.sin(-1/2*(t-t0p)/dtp)
                    # Ez[102, 212] = E0*np.sin(-1/2*(t-t0p)/dtp-np.pi/2)

                    '''Introducimos la posición de las rendijas en función de los parámetros introducidos'''

                    Ez[posx:posx+1, 0:start] = 0
                    ind = start+L
                    for i in range(numrendijas-1):
                        Ez[posx:posx+1, ind:ind+a] = 0
                        ind += a+L
                    Ez[posx:posx+1, ind:] = 0

                    # Calculamos los nuevos campos magnéticos
                    Hx[:, :-1] = Hx[:, :-1] - 1/2*(Ez[:, 1:]-Ez[:, :-1])
                    Hy[:-1, :] = Hy[:-1, :] + 1/2*(Ez[1:, :]-Ez[:-1, :])

                    # Vamos modificando la intensidad que representaremos
                    sumI += I(Ez[posx+d, :])

                ax3=fig2.add_subplot(splot)
                splot += 1 # cambiamos el número de subplot para la siguiente
                ax3.set_title(f'{numrendijas} rendija(s), d={d}, a={a}', y=1.05, weight='bold', size=10)
                ax3.plot(ypos, sumI)

    plt.savefig(namejpg)
    plt.show()


Simulacion(n_s, d_s, a_s)
