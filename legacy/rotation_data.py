# coding: utf-8

#Jose Benavides (jose.astroph@gmail.com)
"""
Recibe como argumentos las posiciones X,Y,Z de un conjunto de datos y dos ángulos para realizar la rotación.
Los ángulos se pasan por referencia como angle_XY_plane y angle_Z_axis, por definición se incluyen como 0 y su valor debe introducirse en grados.

=================

Utilizamos las matrices de rotación presentadas por Beauge C. Notas de Mecánica Celeste 2018 ec.(2.16) pág.71

|X'| | cos(w)  sin(w)  0| |  1     0      0   | | cos(Omega)  sin(Omega)  0| |X|
|Y'|=|-sin(w)  cos(w)  0| |  0  cos(I)  sin(I)| |-sin(Omega)  cos(Omega)  0| |Y|
|Z'| |   0       0     1| |  0 -sin(I)  cos(I)| |      0          0       1| |Z|

donde:	w es el ángulo del árgumento de pericentro
	I es el ángulo de inclinación respecto plano XY, que es equivalente al ángulo formado con el eje Z
	Omega ángulo de rotación en el plano XY

Solo usamos los ángulos I y Omega
"""

import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D


def rotate(X_data,Y_data,Z_data,angle_XY_plane=0, angle_Z_axis=0):
	"""
	Recibe como argumentos las posiciones X,Y,Z de un conjunto de datos y dos ángulos para realizar la rotación.
	Los ángulos se pasan por referencia como angle_XY_plane y angle_Z_axis, por definición se incluyen como 0 y su valor debe introducirse en grados.
	
	Aplicamos las matrices de rotación Todo_sobre_Bauge pag 71 pdf:
	con w=0 y cuya expresión queda de la forma:
	
	x'= x*cos(Omega) + y*sin(Omega)
	y'= -x*cos(I)*sin(Omega) + y*cos(I)*cos(Omega) + z*sin(I)
	z'= x*sin(I)*sin(Omega) - y*sin(I)*cos(Omega) + z*cos(I)
	"""
	Omega_grad = angle_XY_plane
	I_grad = angle_Z_axis

	Omega=Omega_grad*np.pi/180.
	I=I_grad*np.pi/180.
	
	x_prima = X_data*np.cos(Omega) + Y_data*np.sin(Omega)
	y_prima = -X_data*np.cos(I)*np.sin(Omega) + Y_data*np.cos(I)*np.cos(Omega) + Z_data*np.sin(I)
	z_prima = X_data*np.sin(I)*np.sin(Omega) - Y_data*np.sin(I)*np.cos(Omega) + Z_data*np.cos(I)

	POS_prima = np.column_stack((x_prima,y_prima,z_prima))

	return POS_prima


def plot_rotated_system(X_data,Y_data,Z_data,angle_XY_plane=0, angle_Z_axis=0):
	"""
	Recibe como argumentos las posiciones X,Y,Z de un conjunto de datos y dos ángulos para realizar la rotación.
	Los ángulos se pasan por referencia como angle_XY_plane y angle_Z_axis, por definición se incluyen como 0 y su valor debe introducirse en grados.
	Muestra una figura 3D de los puntos rotados.
	"""

	POS_prima = rotate(X_data,Y_data,Z_data,angle_XY_plane=0, angle_Z_axis=0)

	#-------------------- 6: plot distribucion de puntos y plano de ajuste: ----------------------
	fig = plt.figure()
	ax = fig.add_subplot(111, projection='3d')
	ax.set_aspect('equal')
	#ax.set_xlim( -5, 5)
	#ax.set_ylim( -5, 5)
	#ax.set_zlim( -5, 5)
			
	ax.plot(POS_prima[:,0],POS_prima[:,1],POS_prima[:,2],"r.")
	
	ax.set_xlabel("x",fontsize=7)
	ax.set_ylabel("y",fontsize=7)
	ax.set_zlabel("z",fontsize=7)
	plt.savefig("rotated_system.png",dpi=200)

	plt.show()
	plt.close()

