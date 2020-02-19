# coding: utf-8

#Jose Benavides (jose.astroph@gmail.com)
"""
Recibe como argumentos las posiciones X,Y,Z de lo puntos a los que se les desea ajustar un plano.
Devuelve varias características del plano:

solution:		solución de los valores (a,b,c) de la ecuación del plano z = ax + by +c	
dev2:			calcula de desviación cuadrática definida como sumatoria(d_i**2)
angle_to_plane:		genera un vector normal al plano calculado y al plano de referencia XY, luego calcula y devuelve una tupla con el ángulo entre los planos en grados y radianes
distance_to_plane:	calcula la distancia perpendicular de todos los puntos al plano ajustado, devuelve un array con los valores de las distancias
plot_plane:		muestra una figura 3D de los datos y el plano ajustado.	

=================
Deteminamos la ecuación de un plano que mejor ajuste una distribución de puntos.
Buscamos una ecuacion de la forma:
z  = ax + by + c

Debemos ajustar una seria de puntos (xi,yi,zi) tal que su distancia de al plano sea mínima.

f(a,b,c)=sumatoria[( axi + byi + c - zi )**2]

Minimizando:

df/da = 0   ;   df/db = 0   ; df/dc = 0

llegamos a:

a*Sxx + b*Sxy + c*Sx = Sxz
a*Sxy + b*Syy + c*Sy = Syz
a*Sx  + b*Sy  + c*N  = Sz

donde Suv = sumatoria[u_i*v_i] y Su = sumatoria[u_i]

=====================
Example:

data=np.genfromtxt("datos_test3.txt",skip_header=1)
#0:i
#1:Xi
#2:Yi
#3:Zi

X_data=data[:,1]
Y_data=data[:,2]
Z_data=data[:,3]
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
#from matplotlib.patches import ConnectionPatch


def solution(X_data,Y_data,Z_data):
	"""
	Recibe como argumentos las posiciones X,Y,Z de lo puntos a los que se les desea ajustar un plano.
	Devuelve solución de los valores (a,b,c) de la ecuación del plano z = ax + by +c	
	"""

	#def solution(X_data,Y_data,Z_data):

	#sumatorias de los datos (Xi,Yi,Zi):
	Sx=X_data.sum()
	Sy=Y_data.sum()
	Sz=Z_data.sum()
	Sxx=(X_data**2).sum()
	Syy=(Y_data**2).sum()
	Sxy=(X_data*Y_data).sum()
	Sxz=(X_data*Z_data).sum()
	Syz=(Y_data*Z_data).sum()
	N=len(X_data)

	"""
	El sistema de ecuaciones se puede reescribir en términos de una matriz y un vector solución de la forma:

	A[a,b,c]=v

	    |Sxx Sxy Sx|
	A = |Sxy Syy Sy|
	    |Sx  Sy  N |

	    |Sxz|
	v = |Syz|
	    |Sz |
	"""
	#------- 1: Vamos a generar un arreglo para "acomodar" todos los coeficientes de las ecuaciones: -------

	Matriz_coef=np.array([[Sxx,Sxy,Sx],[Sxy,Syy,Sy],[Sx,Sy,N]]) #matriz que almacena los coeficientes: A
	vect_sol=np.array([Sxz,Syz,Sz]) #vector que representa la solucion, es decir el lado derecho del sistema: v

	#------------------------- 2: Ahora vamos a resolver el sistema: -------------------------
	#para ello usaremos "linalg" que es una libreria de algebra lineal, tiene un monto de herramientas que vale la pena revisar:
	#norma, inversa, solucion, determinates, ...
	#tambien tiene algunas cosas de Tensores

	sol_abc=np.linalg.solve(Matriz_coef,vect_sol)  #usaremos la funcion solve para resolver el sistema

	#print "La solución para los valores de (a,b,c) son:\na=",sol_abc[0],"\nb=",sol_abc[1],"\nc=",sol_abc[2]

	a=float("{0:.2f}".format(sol_abc[0]))
	b=float("{0:.2f}".format(sol_abc[1]))
	c=float("{0:.2f}".format(sol_abc[2]))

	#print "\nPor tanto la ecuación del plano es: Z = "+str(a)+" X + "+str(b)+" Y + "+str(c)+"\n"

	return sol_abc


def dev2(X_data,Y_data,Z_data):
	"""
	Recibe como argumentos las posiciones X,Y,Z de lo puntos a los que se les desea ajustar un plano.
	Devuelve la desviación cuadrática de los puntos.
	"""
	sol_abc = solution(X_data,Y_data,Z_data)

	#------------------- 3: Determinamos la desviacion cuadratica d**2: -------------------------
	#calculamos los nuevos valores de zeta: Zi(Xi,Yi) usando la ecuacion del plano:
	Zi_in_plano=sol_abc[0]*X_data + sol_abc[1]*Y_data + sol_abc[2]

	#d**2:
	desv2 = ((Zi_in_plano - Z_data)**2).sum()

	#print "El valor de la desviación cuadrática es:", desv2
	
	return desv2


def angle_to_plane(X_data,Y_data,Z_data):
	"""
	Recibe como argumentos las posiciones X,Y,Z de lo puntos a los que se les desea ajustar un plano.
	Devuelve una tupla con el ángulo que forma el plano que ajusta los puntos con el plano XY.
	El primer valor corresponde al ángulo en radianes, el segundo en grados.
	"""

	sol_abc = solution(X_data,Y_data,Z_data)

	#---------------------------- 4: vector normal y ángulo del plano: ---------------------------
	#para determinar el ángulo entre dos planos debemos definir un vector normal a cada uno y luego calcular el angulo entre estos vectores por su producto escalar
	#Defininos tres puntos A,B,C sobre el plano y con ellos los vectores AB, AC, con estos hallamos su producto vectorial que da, por definición, un vector normal al plano
	#podemos los valores (x,y) => (0,1), (1,1) , (1,0) y calculamos sus respectivos valores Z:
	#pp=0.05
	#XA=np.mean(Y_data)
	#YA=np.mean(Z_data)

	#se pueden colocar cualquier valores, no depende de ello:
	XA=5
	YA=5

	XB=2*XA + 5
	YB=2*YA - 4

	XC=1.5*XA -2
	YC=3*YA+2
	
	A=np.array([XA,YA,sol_abc[0]*XA + sol_abc[1]*YA + sol_abc[2]])
	B=np.array([XB,YB,sol_abc[0]*XB + sol_abc[1]*YB + sol_abc[2]])
	C=np.array([XC,YC,sol_abc[0]*XC + sol_abc[1]*YC + sol_abc[2]])
	ABC=np.vstack((A,B,C))

	V_AB = np.array([B[0] - A[0],B[1] - A[1],B[2] - A[2]])
	V_AC = np.array([C[0] - A[0],C[1] - A[1],C[2] - A[2]])

	#el producto vectorial con la funcion cross de python:
	vector_normal=np.cross(V_AB,V_AC)
	vector_normal=np.cross(V_AB,V_AC) / np.linalg.norm(vector_normal) #vector normal unitario

	#consideramos el vector sobre el plano de referencia dirigido en la direccion z:
	vector_k=np.array([0,0,1])

	#calculamos ahora el angulo entre los vectores a partir del producto escalar, usando las funciones np.dot y np.linalg.norm:
	Theta = np.arccos(np.dot(vector_k,vector_normal) / (np.linalg.norm(vector_k) * np.linalg.norm(vector_normal)))

	Theta_grad=Theta*180./np.pi
	angle_grad=float("{0:.2f}".format(Theta_grad))

	#print "El ángulo estimado con el eje z es:", angle_grad

	return (Theta,Theta_grad)


def distances_to_plane(X_data,Y_data,Z_data):
	"""
	Recibe como argumentos las posiciones X,Y,Z de lo puntos a los que se les desea ajustar un plano.
	Devuelve un array con las distancias de los puntos al plano z = ax + by +c.

	========
	
	La idea es determinar la distancia mínima de cada punto al plano (distancia perpendicular al plano).
	Para ello usamos la descripción mostrada en: https://www.ditutor.com/distancias/punto_plano.html 

	siendo el punto P(Xo,Yo,Zo) y el plano AX + BY + CZ +D =0

	d(P,plano) = |A*Xo + B*Yo + C*Zo + D| / srqt(A**2 + B**2 + C**2)

	siendo:
	A = -a
	B = -b
	C =  1
	D = -c

	"""

	sol_abc = solution(X_data,Y_data,Z_data)

	#------------------------- 5: distancia de los puntos al plano -------------------------------
	#punto P de prueba:
	#Xo = 2.
	#Yo = 4.
	#Zo = 3.

	Xo = X_data
	Yo = Y_data
	Zo = Z_data

	AA = -sol_abc[0]
	BB = -sol_abc[1]
	CC = 1.
	DD = -sol_abc[2]

	dist_P_to_plane = np.abs(AA*Xo + BB*Yo + CC*Zo + DD) / np.sqrt(AA**2 + BB**2 + CC**2)
	#print "Las distancias de cada punto al plano son: \n",dist_P_to_plane

	return dist_P_to_plane


def plot_plane(X_data,Y_data,Z_data):
	"""
	Recibe como argumentos las posiciones X,Y,Z de lo puntos a los que se les desea ajustar un plano.
	Muestra el una figura 3D del plano z = ax + by +c.
	"""

	sol_abc = solution(X_data,Y_data,Z_data)

	#-------------------- 6: plot distribucion de puntos y plano de ajuste: ----------------------
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.set_aspect('equal')
	#ax.set_xlim( -5, 5)
	#ax.set_ylim( -5, 5)
	#ax.set_zlim( -5, 5)
	
	xx, yy = np.meshgrid(np.linspace(min(X_data),max(X_data),10), np.linspace(min(Y_data),max(Y_data),10))
	Z_plano = sol_abc[0]*xx + sol_abc[1]*yy + sol_abc[2]

	ax.plot_surface(xx, yy, Z_plano,color="b", alpha=0.6)
	ax.plot(X_data,Y_data,Z_data,"r.")
	
	#se incluye la localizacion de la galaxia central:
	center=np.vstack(([0,0,0],[0,0,0]))
	ax.plot(center[:,0],center[:,1],center[:,2],"ko")
	#------------------------------------------------	

	ax.set_xlabel("x",fontsize=7)
	ax.set_ylabel("y",fontsize=7)
	ax.set_zlabel("z",fontsize=7)
	plt.savefig("distribucion_de_puntos_y_plano_de_ajuste.png",dpi=200)

	plt.show()
	plt.close()
