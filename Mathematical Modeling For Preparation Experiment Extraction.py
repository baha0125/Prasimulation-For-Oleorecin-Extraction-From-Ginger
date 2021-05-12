# -*- coding: utf-8 -*-
"""
Created on Wed Nov 11 23:35:47 2020

@author: Baha Tegar Ramadhan
"""

import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import odeint

Nr = 30
Ca0 = 1.671
De = 1e-5
eps = 4e-1
kr = 2e-8
Vl = 2e3
R = 2
Nb = 1e2
tspan = np.linspace(0,36000,401)
rspan = np.linspace(0,R,Nr)
dr = R/(Nr-1)

Ca_init = np.zeros(Nr+1)
Ca_init[:Nr] = np.ones(Nr)*Ca0

def func (C,t):
    dcdt = np.zeros(len(C))
    C[0]=(4*C[1]-C[2])/3
    C[-2]= C[-1]
    for i in range (1,Nr-1):
        dcdt[i]=De/eps*((C[i+1]-2*C[i]+C[i-1])/dr**2+2/rspan[i]*(C[i+1]-C[i-1])/(2*dr)-kr/De*C[i]**1.5)
    dcdt[-1]=-De/Vl*4*np.pi*R**2*Nb*(3*C[-2]-4*C[-3]+C[-4])/(2*dr)
    return dcdt

solv = odeint(func,Ca_init,tspan)
Ca_ball = solv[:,:Nr-1]
Ca_f = solv[:,-1]

Ca_ball[:,0]=(4*Ca_ball[:,1]-Ca_ball[:,2])/3
Ca_ball[1:,-1]=Ca_f[1:]

plt.figure(0,figsize=(10,10))
plt.imshow(np.transpose(Ca_ball),cmap="jet",extent=[0,tspan[-1]/3600,R,0],aspect=tspan[-1]/3600/R,interpolation='bicubic')
plt.xlabel("time (hour)")
plt.ylabel("position radial (cm)")
plt.colorbar()

plt.figure(1)
plt.plot(tspan/3600,Ca_f)
plt.ylabel("Konsentrasi (g/cc)")
plt.xlabel("Waktu (hour)")
plt.ylim(0,0.85)
