'''
-------------------------------
Ondas EM y Experimento de Young
-------------------------------
Miguel Durán Vera
'''
# Este código representa solo una configuración

import numpy as np
import matplotlib.pyplot as plt


'''
El objetivo de esta práctica es simular el experimento de Young.
'''


'''
Datos del sistema
'''

c = 3e8 # velocidad de la luz

size = 401 # número de puntos

dx = 2 # Espaciado en el mallado de posición
dy = dx
dt = dx*1e-9/(2*c) # Espaciado temporal
xpos = np.arange(0,size*dx,dx)*1e-3 # en micras
ypos = np.arange(0,size*dy,dy)*1e-3 # en micras
yy,xx = np.meshgrid(ypos,xpos) # mallado de posiciones

Ez = np.zeros((size,size)) # Iniciamos array vacíos para las ondas
Hx = np.zeros((size,size))
Hy = np.zeros((size,size))

sumI = np.zeros(size) # array vacío donde iremos sumando las intensidades

'''
Datos de la simulación
'''

nsim = 1000 # número de simulaciones/pasos temporales
nrep = 10 # número de simulaciones/pasos entre representaciones
trep = 0.01 # tiempo entre representaciones


'''
Datos rendija
'''

numrendijas = 1 # número de rendijas que estudiaremos

a = 20 #distancia entre rendijas (en nº de veces dy)
L = 2  # anchura de las rendijas (en nº de veces dy)
d = 100 # distancia a la pantalla de las rendijas (en nº de veces dx)

namejpg = f'{numrendijas}rend_a{a}_d{d}.jpg' #nombre que pondremos a la imagen resultante

start = int((size-(numrendijas-1)*a-numrendijas*L)/2) # indice inicial para colocar las rendijas
posx = 200 # posicion donde se coloca la pared con la rendijas

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

frontx = 100 # frontera entre ambos medios
fronty = 200

eps0 = 8.854e-12 # F/m
eps1 = 1 
eps2 = 1 # Hecho para tener 4 zonas pero ahora solo dos: izqda y derecha
eps3 = 1
eps4 = 1

coef = np.ones((size,size)) # Array con la permitividad del espacio
coef[:frontx,:fronty] = eps2
coef[:frontx,fronty:] = eps1
coef[frontx:,:fronty] = eps3
coef[frontx:,fronty:] = eps4

sigma1 = 0 
sigma2 = 0 # Hecho para tener 4 zonas pero ahora solo dos: izqda y derecha
sigma3 = 0
sigma4 = 0

sigma = np.ones((size,size)) # Array con la conductividad# Array con la permitividad del espacio
sigma[:frontx,:fronty] = sigma2
sigma[:frontx,fronty:] = sigma1
sigma[frontx:,:fronty] = sigma3
sigma[frontx:,fronty:] = sigma4

ctc = sigma*dt/(2*eps0*coef) # Coeficiente tiempo conductancia

Z0 = 1/(c*eps0) 
I = lambda E, Z = Z0, n = 1: n/(2*Z)*E**2 #Función para el calculo de intensidad


'''
Simulación
'''

levels = np.linspace(-0.1,0.1,21)
cmap = 'cool'

fig1 = plt.figure(figsize=(10, 8)) # Creamos figura para representar las ondas EM
fig1.suptitle('Simulación Ondas EM') 
ax1 = fig1.add_subplot(111)
cs = ax1.contourf(xx,yy,np.clip(Ez,-0.1,0.1),levels, cmap = cmap) # Creamos plot de contorno para representar el mapa en 2D y la intensidad de campo eléctrico con un mapa de colores
bar = plt.colorbar(cs)
ax1.set_xlim(0,size*dx*1e-3)  
ax1.set_ylim(0,size*dx*1e-3)
ax1.set_xlabel('x (micras)')
ax1.set_ylabel('y (micras)')

Ezant = np.zeros((4,2,size)) # Sin fronteras ajustadas a los medios dieléctricos

for sim in range(0,nsim):
    
    t = sim*dt # Actualizamos el tiempo en el que 
    
    Ez[1:,1:] = Ez[1:,1:]*(1-ctc[1:,1:])/(1+ctc[1:,1:]) + 1/(2*coef[1:,1:]*(1+ctc[1:,1:]))*((Hy[1:,1:]-Hy[:-1,1:])-(Hx[1:,1:]-Hx[1:,:-1])) # Calculamos el nuevo campo eléctrico
    

    Ez[:,-1] = Ezant[0,0,:]# Arriba
    Ez[:,0] = Ezant[1,0,:]#Abajo    
    Ez[0,:] = Ezant[2,0,:]# Izqda
    Ez[-1,:] = Ezant[3,0,:]# Drch

    Ezant[0,0,:] = Ezant[0,1,:] # Movemos las copias anteriores una posición
    Ezant[1,0,:] = Ezant[1,1,:]
    Ezant[2,0,:] = Ezant[2,1,:]
    Ezant[3,0,:] = Ezant[3,1,:]
    
    Ezant[0,1,:] = Ez[:,-2] # Cogemos los datos actuales en la ultima copia
    Ezant[1,1,:] = Ez[:,1]
    Ezant[2,1,:] = Ez[1,:]
    Ezant[3,1,:] = Ez[-2,:] 
         
    Ez[10, 50:350] = E0*np.sin(-1/2*(t-t0p)/dtp) # Introducimos el pulso inicial
    # Ez[102, 212] = E0*np.sin(-1/2*(t-t0p)/dtp-np.pi/2)
    
    '''Introducimos la posición de las rendijas en función de los parámetros introducidos'''
    Ez[posx:posx+1, 0:start] = 0 
    ind = start+L
    for i in range(numrendijas-1):
        Ez[posx:posx+1, ind:ind+a] = 0  
        ind += a+L    
    Ez[posx:posx+1, ind:] = 0

    Hx[:,:-1] = Hx[:,:-1] - 1/2*(Ez[:,1:]-Ez[:,:-1]) # Calculamos los nuevos campos magnéticos
    Hy[:-1,:] = Hy[:-1,:] + 1/2*(Ez[1:,:]-Ez[:-1,:])
       
    sumI += I(Ez[posx+d,:]) # Vamos modificando la intensidad que representaremos

    if sim % nrep == 0: # Representamos solo algunas de las simulaciones
        
        ax1.cla()
        cs = ax1.contourf(xx,yy,np.clip(Ez,-0.1,0.1),levels, cmap = cmap)
        title = round(t*1e15)
        ax1.title.set_text(f'Time: {title} fs') # Título que nos indica el tiempo en fs

        '''Dibujamos las rendijas'''
        ax1.axvspan(xpos[posx], xpos[posx+1], ymin=ypos[0]/ypos[-1], ymax=ypos[start]/ypos[-1], color='k')  
        ind = start+L
        for i in range(numrendijas-1):
            ax1.axvspan(xpos[posx], xpos[posx+1], ymin=ypos[ind]/ypos[-1], ymax=ypos[ind+a]/ypos[-1], color='k')
            ind += a+L
        ax1.axvspan(xpos[posx], xpos[posx+1], ymin=ypos[ind]/ypos[-1], ymax=ypos[-1]/ypos[-1], color='k')
        ax1.axvspan(xpos[posx+d], xpos[posx+d+1], ymin=0, ymax=1, color='silver', linestyle = 'dashed')
        
        plt.pause(trep)   

fig2 = plt.figure(figsize=(10,8)) # Figura para la representación de la intensidad
fig2.suptitle(f'{numrendijas} rendija(s),d={d},a={a}', weight = 'bold', size = 16)

ax3 = fig2.add_subplot()
ax3.plot(ypos, sumI)
plt.savefig(namejpg)
plt.show()

    