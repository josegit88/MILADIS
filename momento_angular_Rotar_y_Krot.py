import numpy as np
import matplotlib.pyplot as plt
import Rotate
plt.close()

data_gruposFoF=np.genfromtxt('info_gruposFoF_filtrado.txt')
#0:ID del sh central del FoF
#1:GrNr del grupoFoF a z=0
#2:M200 del FoF
#3:R200 del FoF
#4:Numero de subhalos en el FoF
#5:Numero de subhalos filtrados por Ms >= 1e8 Msun
#6-9: fotometria U,B,V,K

data_gruposFoF_rotated=[]
#for ii in range(100,102):
for ii in range(len(data_gruposFoF)):

	sh_ID=data_gruposFoF[ii,0].astype("int32")
	GrNr_FoF=data_gruposFoF[ii,1].astype("int32")
	M200=data_gruposFoF[ii,2]
	R200=data_gruposFoF[ii,3]
	Nsubh=data_gruposFoF[ii,5].astype("int32")
	phot_U=data_gruposFoF[ii,6]
	phot_B=data_gruposFoF[ii,7]
	phot_V=data_gruposFoF[ii,8]
	phot_K=data_gruposFoF[ii,9]

	print "sh_ID:",sh_ID

	try:
		data_sh=np.load("/media/jbenavides/HDD4T/glx_in_sistemas_tipo_LG/TNG100-1_gal_ID_"+str(sh_ID)+"_particle_type_4.npy")
		#0:particle MASS en 1e10 Msun
		#1-3:particle POS en kpc
		#4-6:particle VEL en km/s
		#7:particle ID
	except:
		print "El sh_ID:",sh_ID,"no tiene estrellas, next..."
		continue

	Rgal=min(0.15*R200,50)
	particle_masses=data_sh[:,0]*1e10
	radial_pos=np.sqrt(data_sh[:,1]**2 + data_sh[:,2]**2 + data_sh[:,3]**2)
	masas_y_radios=np.column_stack((particle_masses,radial_pos))
	what_row_leq_rgal=np.where(masas_y_radios[:,1] <= Rgal)[0]
	data_sh_in_Rgal=data_sh[what_row_leq_rgal,:]
	Mstar_with_Rgal=(data_sh_in_Rgal[:,0]*1e10).sum()
	sort_Rparticles=masas_y_radios[masas_y_radios[:,1].argsort()]
	
	sum_MassParticles=sort_Rparticles[0,0]
	for jj in range(1,len(sort_Rparticles)+1):
		sum_MassParticles = sum_MassParticles + sort_Rparticles[jj,0]
		if sum_MassParticles < 0.5*Mstar_with_Rgal:
			continue
		elif sum_MassParticles >= 0.5*Mstar_with_Rgal:
			RhalfMass=sort_Rparticles[jj,1]
			break

	what_row_RhallMass=np.where(radial_pos <= RhalfMass)[0]
	data_sh_inRadHallMass=data_sh[what_row_RhallMass,:]

	Vx = data_sh_inRadHallMass[:,4]  
	Vy = data_sh_inRadHallMass[:,5]
	Vz = data_sh_inRadHallMass[:,6]
	
	Vx_mean=np.mean(data_sh_inRadHallMass[:,4])
	Vy_mean=np.mean(data_sh_inRadHallMass[:,5])
	Vz_mean=np.mean(data_sh_inRadHallMass[:,6])

	Delta_Vx_cuadr=((Vx - Vx_mean)**2).sum(0)
	Delta_Vy_cuadr=((Vy - Vy_mean)**2).sum(0)
	Delta_Vz_cuadr=((Vz - Vz_mean)**2).sum(0)

	Delta_V=((Delta_Vx_cuadr + Delta_Vy_cuadr + Delta_Vz_cuadr) / (len(data_sh_inRadHallMass)))**0.5
		
	#------------------------------ Calculo de Krot: ---------------------------

	#1) rotar la galaxia en direccion al mayor momento angular (redefinimos pos y vel). Usa la funcion Rotate que te la adjunto:
	in_3rhalf = np.where(radial_pos < 3.*RhalfMass)[0]
	data_sh_in_3Rfalf=data_sh[in_3rhalf,:]

	#*************************************************************************************************
	#calculo del momento angular y angulos de rotacion a lo Jose
	"""
	Momento angular L a partir del producto vectorial:
	L = r x p = r x mv = m[r x v] = m[(Y*Vz - Z*Vy)i + (Z*Vx - X*Vz)j + (X*Vy - Y*Vx)k]
	donde se ha absorvido el signo menos del segundo termino dentro del parentesis
	"""
	Lx_i = (data_sh_in_3Rfalf[:,0]*1e10)*(data_sh_in_3Rfalf[:,2]*data_sh_in_3Rfalf[:,6] - data_sh_in_3Rfalf[:,3]*data_sh_in_3Rfalf[:,5])
	Ly_i = (data_sh_in_3Rfalf[:,0]*1e10)*(data_sh_in_3Rfalf[:,3]*data_sh_in_3Rfalf[:,4] - data_sh_in_3Rfalf[:,1]*data_sh_in_3Rfalf[:,6])
	Lz_i = (data_sh_in_3Rfalf[:,0]*1e10)*(data_sh_in_3Rfalf[:,1]*data_sh_in_3Rfalf[:,5] - data_sh_in_3Rfalf[:,2]*data_sh_in_3Rfalf[:,4])

	L_i = np.column_stack((Lx_i,Ly_i,Lz_i))
	L_tot = L_i.sum(0) #momento angular total

	"""
	Los angulos Omega -> phi y I(inclinacion) -> theta se calcularon siguiendo a Todo_sobre_Bauge

	Omega = arctan(Lx/-Ly)
	Inclinacion = arctan(sqrt(Lx**2 + Ly**2)/Lz)
	"""
	Omega = np.arctan2(L_tot[0] , (-L_tot[1]))
	I = np.arctan2(np.sqrt(L_tot[0]**2 + L_tot[1]**2) , L_tot[2])

	Omega_grad=Omega*180./np.pi
	I_grad=I*180./np.pi
	#*************************************************************************************************


	jx=(data_sh_in_3Rfalf[:,0]*1e10)*(data_sh_in_3Rfalf[:,2]*data_sh_in_3Rfalf[:,6] - data_sh_in_3Rfalf[:,3]*data_sh_in_3Rfalf[:,5])
	jy=(data_sh_in_3Rfalf[:,0]*1e10)*(data_sh_in_3Rfalf[:,3]*data_sh_in_3Rfalf[:,4] - data_sh_in_3Rfalf[:,1]*data_sh_in_3Rfalf[:,6])
	jz=(data_sh_in_3Rfalf[:,0]*1e10)*(data_sh_in_3Rfalf[:,1]*data_sh_in_3Rfalf[:,5] - data_sh_in_3Rfalf[:,2]*data_sh_in_3Rfalf[:,4])

	jx_tot= jx.sum() / (data_sh_in_3Rfalf[:,0]*1e10).sum()
	jy_tot= jy.sum() / (data_sh_in_3Rfalf[:,0]*1e10).sum()
	jz_tot= jz.sum() / (data_sh_in_3Rfalf[:,0]*1e10).sum()

	jtot = np.sqrt(jx_tot**2 + jy_tot**2 + jz_tot**2)

	jx_tot = jx_tot/jtot
	jy_tot = jy_tot/jtot
	jz_tot = jz_tot/jtot	

	jx=0
	jy=0
	jz=0

	# -- compute rotation vectors:
	[e1,e2,e3] = Rotate.exyz(jx_tot,jy_tot,jz_tot)
	npar=len(data_sh[:,0])

	#defino nuevas posiciones y velocidades:
	rot_pos = np.zeros(shape=(npar,3))
	rot_vel = np.zeros(shape=(npar,3))

	[rot_pos[:,0],rot_pos[:,1],rot_pos[:,2]] = Rotate.rotate(npar,data_sh[:,1],data_sh[:,2],data_sh[:,3],e1,e2,e3)
	[rot_vel[:,0],rot_vel[:,1],rot_vel[:,2]] = Rotate.rotate(npar,data_sh[:,4],data_sh[:,5],data_sh[:,6],e1,e2,e3)

	#2) Usando las nuevas posiciones y velocidades, calculamos el Krot:
	#Nos quedamos solo con las particulas adentro de Rgal:
	son = np.where(radial_pos < Rgal)[0]
	mass = data_sh[son,0]*1e10
	rot_pos = rot_pos[son,:]
	rot_vel = rot_vel[son,:]
	npar = len(mass)

	# kinetic energy in ordered rotation
	jz = np.zeros(npar)
	E_rot = np.zeros(npar)
	jz[:] = rot_pos[:,0]*rot_vel[:,1] - rot_pos[:,1]*rot_vel[:,0]
	iii = (rot_pos[:,0] == 0)
	rot_pos[iii,:] = 1.e-10
	jz[iii] = 0.
	E_rot[:] = mass[:] * jz[:]**2 / (rot_pos[:,0]**2 + rot_pos[:,1]**2)

	# kinetic energy (total):
	E_kin = np.zeros(npar)
	E_kin[:] = mass[:]*(rot_vel[:,0]**2 + rot_vel[:,1]**2 + rot_vel[:,2]**2)

	# defino el kappa:
	E_rot_Tot=np.sum(E_rot[:],dtype="float64")
	E_kin_Tot=np.sum(E_kin[:],dtype="float64")

	kappa = E_rot_Tot / E_kin_Tot
	#---------------------------------------------------------------------------


	data_gruposFoF_rotated.append((sh_ID,GrNr_FoF,M200,R200,Nsubh,phot_U,phot_B,phot_V,phot_K,Rgal,RhalfMass,mass.sum(),Delta_V,Omega_grad,I_grad,E_kin_Tot,E_rot_Tot,kappa))
	

data_gruposFoF_rotated=np.array(data_gruposFoF_rotated)
np.savetxt("info_data_gruposFoF_y_rotated.txt",data_gruposFoF_rotated, fmt=["%7.f","%7.f","%13.6e","%13.6e","%4.f","%8.2f","%8.2f","%8.2f","%8.2f","%8.2f","%8.2f","%13.6e","%8.2f","%8.2f","%8.2f","%13.6e","%13.6e","%8.4f"])
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


#---- Energia de rotacion normalizada vs Mstar: ----------
fig = plt.figure(figsize=(6,6))
llss=10
tMl=7
tml=4
#plt.locator_params(axis="y",tight=True, nbins=2)
plt.tick_params(which="major",length=tMl,top=True, left=True, right=True,direction="in",labelsize=llss)

plt.minorticks_on()
plt.tick_params(which="minor",length=tml,top=True, left=True, right=True,direction="in",labelsize=llss)

plt.locator_params(axis="x",tight=True, nbins=8)
plt.locator_params(axis="y",tight=True, nbins=8)

plt.plot(np.log10(data_gruposFoF_rotated[:,11]),data_gruposFoF_rotated[:,17], color="r", marker="o",fillstyle="none",ms=5,ls="none",label=r"Sample Groups: $\rm{12 \, \leq log(M_{200} / M_{\odot}) \, \leq 13}$")
#plt.plot([0,4],[0,4],"k-",lw=1)

plt.legend(loc="lower left",ncol=1,frameon=False,fontsize=9,handletextpad=0.001)
plt.xlabel(r"$\rm{log(M^\bigstar_{R_{gal}}) \ [M_{\odot}]}$",fontsize=14)
plt.ylabel(r"$\it{ \mathcal{K}_{rot}}$",fontsize=14)
#plt.xlim(1.45,3.3)
#plt.ylim(1.45,3.3)

plt.savefig("plots/Mstar_in_Rgal_vs_Krot_grupos_Sample.png",bbox_inches='tight',dpi=200)
plt.show()
plt.close()

#---- Energia de rotacion normalizada vs Color B-V: ----------
fig = plt.figure(figsize=(6,6))
llss=10
tMl=7
tml=4
#plt.locator_params(axis="y",tight=True, nbins=2)
plt.tick_params(which="major",length=tMl,top=True, left=True, right=True,direction="in",labelsize=llss)

plt.minorticks_on()
plt.tick_params(which="minor",length=tml,top=True, left=True, right=True,direction="in",labelsize=llss)

plt.locator_params(axis="x",tight=True, nbins=8)
plt.locator_params(axis="y",tight=True, nbins=8)

plt.plot(data_gruposFoF_rotated[:,6] - data_gruposFoF_rotated[:,7] ,data_gruposFoF_rotated[:,17], color="r", marker="o",fillstyle="none",ms=5,ls="none",label=r"Sample Groups: $\rm{12 \, \leq log(M_{200} / M_{\odot}) \, \leq 13}$")
#plt.plot([0,4],[0,4],"k-",lw=1)

plt.legend(loc="lower left",ncol=1,frameon=False,fontsize=9,handletextpad=0.001)
plt.xlabel(r"B - V",fontsize=14)
plt.ylabel(r"$\it{ \mathcal{K}_{rot}}$",fontsize=14)
#plt.xlim(1.45,3.3)
#plt.ylim(1.45,3.3)

plt.savefig("plots/color_vs_Krot_grupos_Sample.png",bbox_inches='tight',dpi=200)
plt.show()
plt.close()

#---------------------------------
print "--------------------\nEl programa ha terminado con exito, muy bien Jose"
