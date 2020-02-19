
# coding: utf-8

# In[5]:

"""
Functions to rotate positions to a system oriented according to lx,ly,lz 
exyz: creates the rotation vectors e1, e2, e3.
     Input: lx,ly,lz : x,y,z of the vector in the direction of wanted rotation
     Output: e1(3), e2(3), e3(3): rotation matrix
rotate: rotates single or multiple coordinates to new system.
     Input: n = dimension of array to rotate
            x(n),y(n),z(n) = components in old coordinate system (should be np.arrays)
     Output: rx(n),ry(n),rz(n) = coordinates rotated
"""
#funcion Rotate

import numpy as np

def exyz(lx,ly,lz):
    e1 = np.zeros(3)
    e2 = np.zeros(3)
    e3 = np.zeros(3)

    e1x =  ly*lz
    e1y = -lx*lz
    e1z = 0.
    me1  = np.sqrt(e1x**2 + e1y**2 + e1z**2)
    if me1 > 0:
        e1x = e1x/me1
        e1y = e1y/me1
        e1z = e1z/me1
    e1[0]=e1x
    e1[1]=e1y
    e1[2]=e1z

    e2x = lx*lz**2
    e2y = ly*lz**2
    e2z = -(lx**2+ly**2)*lz
    me2 = np.sqrt(e2x**2 + e2y**2 + e2z**2)
    if me2 > 0: 
        e2x = e2x/me2 
        e2y = e2y/me2 
        e2z = e2z/me2 
    e2[0] = e2x
    e2[1] = e2y
    e2[2] = e2z

    l = np.sqrt(lx**2 + ly**2 + lz**2)
    if (l > 0):
        lx = lx/l
        ly = ly/l
        lz = lz/l

    e3x = lx
    e3y = ly
    e3z = lz      # already normalized
    e3[0] = e3x
    e3[1] = e3y
    e3[2] = e3z

    return [e1,e2,e3]    



def rotate(n,x,y,z,e1,e2,e3):

    rx = np.zeros(n)
    ry = np.zeros(n)
    rz = np.zeros(n)
    
    rx[:] = e1[0]*x[:] + e1[1]*y[:] + e1[2]*z[:]
    ry[:] = e2[0]*x[:] + e2[1]*y[:] + e2[2]*z[:]
    rz[:] = e3[0]*x[:] + e3[1]*y[:] + e3[2]*z[:]

#    rx[:] = e1[0]*x[:] + e2[0]*y[:] + e3[0]*z[:]
#    ry[:] = e1[1]*x[:] + e2[1]*y[:] + e3[1]*z[:]
#    rz[:] = e1[2]*x[:] + e2[2]*y[:] + e3[2]*z[:]
 
    return [rx, ry, rz]

# -- computes the angular momentum using pos and vel, rotates and returns rotated pos vel ----#
def rotate_pos_vel(pos, vel, mass):

    npar = len(mass)
    angmom = np.cross(pos*np.reshape(mass,(npar,1)),vel)    #LVS: provided mass is shape=(len,1)
    angmomtot = np.sum(angmom, axis=0)
    mtot = np.sum(mass,dtype="float64")
    angmomtot /= mtot

    [e1,e2,e3] = exyz(angmomtot[0],angmomtot[1],angmomtot[2])

    [rotx,roty,rotz] = rotate(npar,pos[:,0],pos[:,1],pos[:,2],e1,e2,e3)
    [rotvx,rotvy,rotvz] = rotate(npar,vel[:,0],vel[:,1],vel[:,2],e1,e2,e3)

    return[ np.array([rotx,roty,rotz]).T ,  np.array([rotvx,rotvy,rotvz]).T]
