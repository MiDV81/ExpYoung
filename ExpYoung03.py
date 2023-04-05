'''
--------------------
Ondas EM
--------------------
Miguel Durán Vera
'''


import numpy as np
import matplotlib.pyplot as plt


'''
El objetivo de esta práctica es simular un pulso electromagnetico propagandose por un plano 2D en el vacio.
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
yy,xx = np.meshgrid(ypos,xpos)

print(ypos[400])
Ez = np.zeros((size,size)) # Iniciamos array vacíos para las ondas
Hx = np.zeros((size,size))
Hy = np.zeros((size,size))

sumI = np.zeros(size)
'''
Datos de la simulación
'''

nsim = 1000 # número de simulaciones/pasos temporales
nrep = 10 # número de simulaciones/pasos entre representaciones
trep = 0.01 # tiempo entre representaciones


'''
Gaussiana Inicial
'''

E0 = 1 # 
dxp = 40 # Anchura del pulso
dtp = dxp*1e-9/c # desviación típica de la gaussiana
t0p = 5*dtp # Tiempo inicial de la gaussiana
pos0px = 100 # Posición donde metemos el pulso x
pos0py = 200 # Posición donde metemos el pulso y

'''
Ecs. propagación
'''

'''
Ey(k) = Ey(k)-1/2*(Hz(k)-Hz(k-1))
Hz(k) = Hz(k)-1/2*(Ey(k+1)-Ey(k))
'''

'''
Datos Medios
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
I = lambda E, Z = Z0, n = 1: n/(2*Z)*E**2
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
# ax1.axvline(frontx*dx*1e-3, color = 'black', linestyle = 'dashed')

# ax2 = fig1.add_subplot(212)
# # ax2.set_ylim(0, 0.2)
# pantalla0, = ax2.plot(ypos,Ez[150,:])
# pantalla1, = ax2.plot(ypos,Ez[200,:])
# pantalla2, = ax2.plot(ypos, Ez[250, :])


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
         
    Ez[10, 50:350] = E0*np.sin(-1/2*(t-t0p)/dtp)
    # Ez[102, 212] = E0*np.sin(-1/2*(t-t0p)/dtp-np.pi/2)
    
    Ez[200:201, 0:188] = 0
    Ez[200:201, 190:210] = 0
    Ez[200:201, 212:] = 0

    Hx[:,:-1] = Hx[:,:-1] - 1/2*(Ez[:,1:]-Ez[:,:-1]) # Calculamos los nuevos campos magnéticos
    Hy[:-1,:] = Hy[:-1,:] + 1/2*(Ez[1:,:]-Ez[:-1,:])
       
    sumI += I(Ez[250,:])

    if sim % nrep == 0: # Representamos solo algunas de las simulaciones
        
        ax1.cla()
        cs = ax1.contourf(xx,yy,np.clip(Ez,-0.1,0.1),levels, cmap = cmap)
        title = round(t*1e15)
        ax1.title.set_text(f'Time: {title} fs') # Título que nos indica el tiempo en fs
        # ax1.axvline(frontx*dx*1e-3, color = 'black', linestyle = 'dashed') # Frontera entre las zonas
        # ax1.axhline(fronty*dy*1e-3, color = 'black', linestyle = 'dashed')
        # ax2.set_ylim(0, np.max(I(Ez[150, :])*1e10+1e-10))
        # pantalla0.set_data(ypos, I(Ez[150,:])*1e10)
        # pantalla1.set_data(ypos, I(Ez[200,:])*1e10)
        # pantalla2.set_data(ypos, I(Ez[250, :])*1e10)
        ax1.axvspan(xpos[200], xpos[201], ymin=ypos[0]/ypos[400],ymax=ypos[188]/ypos[400], color='k')       
        ax1.axvspan(xpos[200], xpos[201], ymin=ypos[190]/ypos[400], ymax=ypos[210]/ypos[400], color='k')
        ax1.axvspan(xpos[200], xpos[201], ymin=ypos[212]/ypos[400], ymax = ypos[400]/ypos[400], color='k')

        plt.pause(trep)   


fig2 = plt.figure(figsize=(10,8))
fig2.suptitle('Pantallas')

ax3 = fig2.add_subplot()
ax3.plot(ypos, sumI)

plt.show()

    