# coding: utf-8

import numpy as np
import matplotlib.pyplot as plt
import find_plane
import rotation_data 

data_gruposFoF=np.genfromtxt('info_data_gruposFoF_y_rotated.txt')
#0:ID
#1:GrNr
#2:M200 en Msun
#3:R200
#4:Nsubh
#5-8:fotometria en los filtros U,B,V,K
#9:Rgal en kpc
#10:RhalfMass calculado en kpc
#11:Mstar calculada con Rgal en Msun
#12:Dispersion de velocidades in RhalfMass en km/s
#13:angulo Omega en grados
#14:angulo I en grados
#15:Energia cinetica total in Rgal in Msun*km**2/s**2
#16:Energia rotacional total in Rgal Msun*km**2/s**2
#17:Parametro de rotacion normalizado Kappa_rot

#ii=3

median_distances_glx=[]
#for ii in range(10):
for ii in range(len(data_gruposFoF)):
	#print "ii:",ii
	
	#if data_glxs.shape
	sh_ID=data_gruposFoF[ii,0].astype("int32")
	GrNr=data_gruposFoF[ii,1].astype("int32")
	M200=data_gruposFoF[ii,2]
	R200=data_gruposFoF[ii,3]
	N_subh=data_gruposFoF[ii,4].astype("int32")
	Ms_glx_host=data_gruposFoF[ii,11]
	dispVEL=data_gruposFoF[ii,12] #dispersion de velocidades
	Omega=data_gruposFoF[ii,13]
	I=data_gruposFoF[ii,14]
	kappa=data_gruposFoF[ii,17]

	if N_subh <= 3:
		continue

	data_glxs=np.genfromtxt("info_subhalos_gruposFoF/info_sh_grupo_FoF_"+str(data_gruposFoF[ii,1].astype("int32"))+".txt")
	#0:sh_ID
	#1:Mstar en M_sun
	#2:Mgas en M_sun
	#3:DM en M_sun
	#4:SFR
	#5,6,7:POS x,y,z en kpc
	#8,9,10:VEL Vx,Vy,Vz en km/s
	#11-18: fotometria U,B,V,K,g,r,i,z

	#---------- datos centrados en la primera galaxia de la lista: -----------------------
	X=data_glxs[:,5] - data_glxs[0,5]
	Y=data_glxs[:,6] - data_glxs[0,6]
	Z=data_glxs[:,7] - data_glxs[0,7]

	r_POS=np.sqrt(X**2 + Y**2 + Z**2)
	rows_near1=np.where(r_POS <= 300)[0]
	glx_near_first=data_glxs[rows_near1,:]

	X1=glx_near_first[:,5] - glx_near_first[0,5]
	Y1=glx_near_first[:,6] - glx_near_first[0,6]
	Z1=glx_near_first[:,7] - glx_near_first[0,7]

	#hace falta rotar los datos!!!
	#newPOS=rotation_data.rotate(X,Y,Z,Omega,I)
	#rotation_data.plot_rotated_system(X,Y,Z,Omega,I) #realiza un plot de los daos rotados

	newPOS1=rotation_data.rotate(X1,Y1,Z1,Omega,I)
	if len(newPOS1) <= 3:
		#print "hay menos de tres puntos, no se puede obtener un plano"
		median_distances1 = -1
		std_distances1 = -1
	elif len(newPOS1) > 3:
		#find_plane.angle_to_plane(newPOS1[:,0],newPOS1[:,1],newPOS1[:,2])
		#find_plane.plot_plane(newPOS1[:,0],newPOS1[:,1],newPOS1[:,2])
	
		median_distances1=np.median(find_plane.distances_to_plane(newPOS1[:,0],newPOS1[:,1],newPOS1[:,2]))
		std_distances1=np.std(find_plane.distances_to_plane(newPOS1[:,0],newPOS1[:,1],newPOS1[:,2]))
		angle1=find_plane.angle_to_plane(newPOS1[:,0],newPOS1[:,1],newPOS1[:,2])[1]

	#---------- datos centrados en la segunda galaxia de la lista: -----------------------
	X=data_glxs[:,5] - data_glxs[1,5]
	Y=data_glxs[:,6] - data_glxs[1,6]
	Z=data_glxs[:,7] - data_glxs[1,7]

	r_POS=np.sqrt(X**2 + Y**2 + Z**2)
	rows_near2=np.where(r_POS <= 300)[0]
	glx_near_second=data_glxs[rows_near2,:]

	X2=glx_near_second[:,5] - glx_near_second[0,5]
	Y2=glx_near_second[:,6] - glx_near_second[0,6]
	Z2=glx_near_second[:,7] - glx_near_second[0,7]

	newPOS2=rotation_data.rotate(X2,Y2,Z2,Omega,I)
	if len(newPOS2) <= 3:
		#print "hay menos de tres puntos, no se puede obtener un plano"
		median_distances2 = -1
		std_distances2 = -1
	elif len(newPOS2) > 3:
		#find_plane.angle_to_plane(newPOS2[:,0],newPOS2[:,1],newPOS2[:,2])
		#find_plane.plot_plane(newPOS2[:,0],newPOS2[:,1],newPOS2[:,2])

		median_distances2=np.median(find_plane.distances_to_plane(newPOS2[:,0],newPOS2[:,1],newPOS2[:,2]))
		std_distances2=np.std(find_plane.distances_to_plane(newPOS2[:,0],newPOS2[:,1],newPOS2[:,2]))

	median_distances_glx.append((GrNr,median_distances1,median_distances2,std_distances1,std_distances2,N_subh,len(newPOS1),Ms_glx_host,kappa,angle1))


median_distances_glx=np.array(median_distances_glx)
#0:GrNr
#1:mediana de las distancias para la primera galaxia
#2:mediana de las distancias para la segunda galaxia
#3:dispersion de las distancias para la primera galaxia
#4:dispersion de las distancias para la segunda galaxia
#5:numero de galaxias en el grupo
#6:numero de galaxias dentro del radio limitado (300 kpc)
#7:Mstar de la primera galaxia central
#8:parametro de rotacion normalizado
#9:angulo entre el plano de la galaxia y el plano ajustado, en grados

#hacemos un corte en Krot >= 0.5 para quedarnos con galaxias preferentemente rotantes:
median_distances_glx = median_distances_glx[median_distances_glx[:,8]>=0.5]


#-----------------histogramas de los angulos entre el plano de galaxias y el plano ajustado plano: ------------------

#First galaxy
fig = plt.figure(figsize=(6,6))

llss=10
tMl=5
tml=3
plt.locator_params(axis="y",tight=True, nbins=6)
plt.tick_params(which="major",length=tMl,top=True, left=True, right=True,direction="in",labelsize=llss)

plt.minorticks_on()
plt.tick_params(which="minor",length=tml,top=False,bottom=False, left=True, right=True,direction="in",labelsize=llss)
#ax[0,0].locator_params(axis="both",tight=True, nbins=7)

plt.locator_params(axis="x",tight=True, nbins=11)
plt.locator_params(axis="y",tight=True, nbins=11)

first_plane=median_distances_glx[median_distances_glx[:,1]>0]

#plt.hist(color_gi_all,histtype="step",log=False,bins=17,range=[-0.3,1.4],color="grey",linestyle="-",orientation="horizontal")
plt.hist(first_plane[:,9],histtype="step",bins=10,range=[0,100],color="red",linestyle="-",label=r"Respect First plane")

plt.legend(loc="upper right",ncol=1,frameon=False,fontsize=10)
plt.xlabel(r"angle between planes [grad]",fontsize=14)
plt.ylabel(r"$\rm{N}$",fontsize=14)
#plt.xlim(-2,102)
#plt.ylim(0,10)

plt.savefig("plots/histograma_angle_between_planes_for_first_galaxy_rotantes.png",bbox_inches='tight',dpi=200)
plt.show()
plt.close()

#-----------------Relacion entre el numero de galaxias y los angulos entre el plano de galaxias y el plano ajustado plano: ------------------

#First galaxy
fig = plt.figure(figsize=(6,6))

llss=10
tMl=5
tml=3
plt.locator_params(axis="y",tight=True, nbins=6)
plt.tick_params(which="major",length=tMl,top=True, left=True, right=True,direction="in",labelsize=llss)

plt.minorticks_on()
plt.tick_params(which="minor",length=tml,top=True,bottom=True, left=True, right=True,direction="in",labelsize=llss)
#ax[0,0].locator_params(axis="both",tight=True, nbins=7)

plt.locator_params(axis="x",tight=True, nbins=11)
plt.locator_params(axis="y",tight=True, nbins=11)

first_plane=median_distances_glx[median_distances_glx[:,1]>0]

plt.plot(first_plane[:,6],first_plane[:,9],"r.")

plt.legend(loc="upper right",ncol=1,frameon=False,fontsize=10)
plt.xlabel(r"$\rm{N_{gal \, in \, group}}$",fontsize=14)
plt.ylabel(r"angle between planes [grad]",fontsize=14)
#plt.xlim(-2,102)
#plt.ylim(0,10)

plt.savefig("plots/relacion_Ngal_in_group_y_angle_between_planes_for_first_galaxy_rotantes.png",bbox_inches='tight',dpi=200)
plt.show()
plt.close()


#-----------------Relacion entre el las distancias y los angulos entre el plano de galaxias y el plano ajustado plano: ------------------

#First galaxy
fig = plt.figure(figsize=(6,6))

llss=10
tMl=5
tml=3
plt.locator_params(axis="y",tight=True, nbins=6)
plt.tick_params(which="major",length=tMl,top=True, left=True, right=True,direction="in",labelsize=llss)

plt.minorticks_on()
plt.tick_params(which="minor",length=tml,top=True,bottom=True, left=True, right=True,direction="in",labelsize=llss)
#ax[0,0].locator_params(axis="both",tight=True, nbins=7)

plt.locator_params(axis="x",tight=True, nbins=11)
plt.locator_params(axis="y",tight=True, nbins=11)

first_plane=median_distances_glx[median_distances_glx[:,1]>0]

plt.plot(first_plane[:,1],first_plane[:,9],"r.")

plt.legend(loc="upper right",ncol=1,frameon=False,fontsize=10)
plt.xlabel(r"$\rm{median(d_{gal-plane}) [kpc]}$",fontsize=14)
plt.ylabel(r"angle between planes [grad]",fontsize=14)
#plt.xlim(-2,102)
#plt.ylim(0,10)

plt.savefig("plots/relacion_distancias_y_angle_between_planes_for_first_galaxy_rotantes.png",bbox_inches='tight',dpi=200)
plt.show()
plt.close()

#-----------------Relacion entre el las distancias al plano ajustado plano y el numero de miembros: ------------------

#First galaxy
fig = plt.figure(figsize=(6,6))

llss=10
tMl=5
tml=3
plt.locator_params(axis="y",tight=True, nbins=6)
plt.tick_params(which="major",length=tMl,top=True, left=True, right=True,direction="in",labelsize=llss)

plt.minorticks_on()
plt.tick_params(which="minor",length=tml,top=True,bottom=True, left=True, right=True,direction="in",labelsize=llss)
#ax[0,0].locator_params(axis="both",tight=True, nbins=7)

plt.locator_params(axis="x",tight=True, nbins=11)
plt.locator_params(axis="y",tight=True, nbins=11)

first_plane=median_distances_glx[median_distances_glx[:,1]>0]

plt.plot(first_plane[:,6],first_plane[:,1],"r.")

plt.legend(loc="upper right",ncol=1,frameon=False,fontsize=10)
plt.xlabel(r"$\rm{N_{gal \, in \, group}}$",fontsize=14)
plt.ylabel(r"$\rm{median(d_{gal-plane}) [kpc]}$",fontsize=14)
#plt.xlim(-2,102)
#plt.ylim(0,10)

plt.savefig("plots/relacion_distancias_y_N_miembros_for_first_galaxy_rotantes.png",bbox_inches='tight',dpi=200)
plt.show()
plt.close()

"""
#-----------------histogramas de las medianas de las distancias de las galaxias al plano: ------------------

#First galaxy
fig = plt.figure(figsize=(6,6))

llss=10
tMl=5
tml=3
plt.locator_params(axis="y",tight=True, nbins=6)
plt.tick_params(which="major",length=tMl,top=True, left=True, right=True,direction="in",labelsize=llss)

plt.minorticks_on()
plt.tick_params(which="minor",length=tml,top=False,bottom=False, left=True, right=True,direction="in",labelsize=llss)
#ax[0,0].locator_params(axis="both",tight=True, nbins=7)

plt.locator_params(axis="x",tight=True, nbins=11)
plt.locator_params(axis="y",tight=True, nbins=11)

first_plane=median_distances_glx[median_distances_glx[:,1]>0]

#plt.hist(color_gi_all,histtype="step",log=False,bins=17,range=[-0.3,1.4],color="grey",linestyle="-",orientation="horizontal")
plt.hist(first_plane[:,1],histtype="step",bins=10,range=[0,100],color="red",linestyle="-",label=r"Distances respect First plane")

plt.legend(loc="upper right",ncol=1,frameon=False,fontsize=10)
plt.xlabel(r"$\rm{median(d_{gal-plane}) [kpc]}$",fontsize=14)
plt.ylabel(r"$\rm{N}$",fontsize=14)
plt.xlim(-2,102)
#plt.ylim(0,10)

#plt.savefig("plots/histograma_mediana_distancias_galaxies_at_first_galaxy.png",bbox_inches='tight',dpi=200)
plt.savefig("plots/histograma_mediana_distancias_galaxies_at_first_galaxy_rotantes.png",bbox_inches='tight',dpi=200)
plt.show()
plt.close()
#----------------------------------------------------

#First galaxy
fig = plt.figure(figsize=(6,6))

llss=10
tMl=5
tml=3
plt.locator_params(axis="y",tight=True, nbins=6)
plt.tick_params(which="major",length=tMl,top=True, left=True, right=True,direction="in",labelsize=llss)

plt.minorticks_on()
plt.tick_params(which="minor",length=tml,top=False,bottom=False, left=True, right=True,direction="in",labelsize=llss)
#ax[0,0].locator_params(axis="both",tight=True, nbins=7)

plt.locator_params(axis="x",tight=True, nbins=11)
plt.locator_params(axis="y",tight=True, nbins=11)

second_plane=median_distances_glx[median_distances_glx[:,2]>0]

#plt.hist(color_gi_all,histtype="step",log=False,bins=17,range=[-0.3,1.4],color="grey",linestyle="-",orientation="horizontal")
plt.hist(second_plane[:,2],histtype="step",bins=10,range=[0,100],color="red",linestyle="-",label=r"Distances respect Second plane")

plt.legend(loc="upper right",ncol=1,frameon=False,fontsize=10)
plt.xlabel(r"$\rm{median(d_{gal-plane}) [kpc]}$",fontsize=14)
plt.ylabel(r"$\rm{N}$",fontsize=14)
plt.xlim(-2,102)
#plt.ylim(0,10)

#plt.savefig("plots/histograma_mediana_distancias_galaxies_at_second_galaxy.png",bbox_inches='tight',dpi=200)
plt.savefig("plots/histograma_mediana_distancias_galaxies_at_second_galaxy_rotantes.png",bbox_inches='tight',dpi=200)
plt.show()
plt.close()


#-----------------histogramas de las dispersiones de las distancias de las galaxias al plano: ------------------

#First galaxy
fig = plt.figure(figsize=(6,6))

llss=10
tMl=5
tml=3
plt.locator_params(axis="y",tight=True, nbins=6)
plt.tick_params(which="major",length=tMl,top=True, left=True, right=True,direction="in",labelsize=llss)

plt.minorticks_on()
plt.tick_params(which="minor",length=tml,top=False,bottom=False, left=True, right=True,direction="in",labelsize=llss)
#ax[0,0].locator_params(axis="both",tight=True, nbins=7)

plt.locator_params(axis="x",tight=True, nbins=11)
plt.locator_params(axis="y",tight=True, nbins=11)

first_plane=median_distances_glx[median_distances_glx[:,3]>0]

#plt.hist(color_gi_all,histtype="step",log=False,bins=17,range=[-0.3,1.4],color="grey",linestyle="-",orientation="horizontal")
plt.hist(first_plane[:,3],histtype="step",bins=10,range=[0,100],color="red",linestyle="-",label=r"Distances respect First plane")

plt.legend(loc="upper right",ncol=1,frameon=False,fontsize=10)
plt.xlabel(r"$\rm{\sigma(d_{gal-plane}) [kpc]}$",fontsize=14)
plt.ylabel(r"$\rm{N}$",fontsize=14)
plt.xlim(-2,102)
#plt.ylim(0,10)

#plt.savefig("plots/histograma_dispersion_distancias_galaxies_at_first_galaxy.png",bbox_inches='tight',dpi=200)
plt.savefig("plots/histograma_dispersion_distancias_galaxies_at_first_galaxy_rotantes.png",bbox_inches='tight',dpi=200)
plt.show()
plt.close()
#----------------------------------------------------

#First galaxy
fig = plt.figure(figsize=(6,6))

llss=10
tMl=5
tml=3
plt.locator_params(axis="y",tight=True, nbins=6)
plt.tick_params(which="major",length=tMl,top=True, left=True, right=True,direction="in",labelsize=llss)

plt.minorticks_on()
plt.tick_params(which="minor",length=tml,top=False,bottom=False, left=True, right=True,direction="in",labelsize=llss)
#ax[0,0].locator_params(axis="both",tight=True, nbins=7)

plt.locator_params(axis="x",tight=True, nbins=11)
plt.locator_params(axis="y",tight=True, nbins=11)

second_plane=median_distances_glx[median_distances_glx[:,4]>0]

#plt.hist(color_gi_all,histtype="step",log=False,bins=17,range=[-0.3,1.4],color="grey",linestyle="-",orientation="horizontal")
plt.hist(second_plane[:,4],histtype="step",bins=10,range=[0,100],color="red",linestyle="-",label=r"Distances respect Second plane")

plt.legend(loc="upper right",ncol=1,frameon=False,fontsize=10)
plt.xlabel(r"$\rm{\sigma(d_{gal-plane}) [kpc]}$",fontsize=14)
plt.ylabel(r"$\rm{N}$",fontsize=14)
plt.xlim(-2,102)
#plt.ylim(0,10)

#plt.savefig("plots/histograma_dispersion_distancias_galaxies_at_second_galaxy.png",bbox_inches='tight',dpi=200)
plt.savefig("plots/histograma_dispersion_distancias_galaxies_at_second_galaxy_rotantes.png",bbox_inches='tight',dpi=200)
plt.show()
plt.close()
"""

"""
#calculamos el R200 de la galaxia central a partir de la dispersion de velocidades y la DarkMatter+Mgas+Mstar:
Mass=data_glxs[0,1]+data_glxs[0,2]+data_glxs[0,3]
sigmaVEL = dispVEL*(1e3)/(3.086e19) #conversion kpc/s
G_const=4.514e-39

Rvir_calc=Mass*G_const/(sigmaVEL**2)

#find_plane.plot_plane(X,Y,Z)
#find_plane.plot_plane(newPOS[:,0],newPOS[:,1],newPOS[:,2])
#find_plane.angle_to_plane(newPOS[:,0],newPOS[:,1],newPOS[:,2])
"""

#---------------------------------
print "--------------------\nEl programa ha terminado con exito, muy bien Jose"
