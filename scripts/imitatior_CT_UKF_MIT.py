import time
import numpy as np
import math
import matplotlib.pyplot as plt
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from scipy.stats import chi2, poisson, uniform, norm
from IPython.display import clear_output
import random

from models import Target
from filterpy.kalman import UnscentedKalmanFilter
from filterpy.kalman import MerweScaledSigmaPoints
from filterpy.common import Q_discrete_white_noise

from copy import deepcopy
from math import log, exp, sqrt
import sys
import numpy as np
from numpy import eye, zeros, dot, isscalar, outer
from scipy.linalg import cholesky
from filterpy.kalman import unscented_transform
from filterpy.stats import logpdf
from filterpy.common import pretty_str

dt = 0.25
pd = 1.0

R = np.diag([10000.0, (1.146496815*(np.pi/180.0))**2, (1.146496815*(np.pi/180.0)**2)]) # дисперсии

plt.rcParams['figure.figsize'] = [10, 6]
fig2 = make_subplots(rows=1, cols=1, specs=[[{'type': 'scatter3d'}]])

   

# =============== Блок 1 ===================
    
# ИНИЦИАЛИЗАЦИЯ МОДЕЛИ ДВИЖЕНИЯ
tg1 = Target()

tg1.init_state({'x':0.0,'y':0.0,'z':0.0, 'vx':200.0,'vy':0.0,'vz':0.0,'w':0.098})
# tg1.init_state({'x':0.0,'y':0.0,'z':0.0, 'vx':200.0,'vy':0.0,'vz':0.0,'w':0.098})


def remove_zero_columns(arr):

    zero_columns = np.all(arr == 0, axis=0)
        
    # Удаляем столбцы, заполненные нулями
    arr = arr[:, ~zero_columns]
        
    return arr

def make_pass(X_true_data, pd): #при pd = 1 пропусков не будет
    xtmp = X_true_data.copy()
    
    for i in range (2,xtmp.shape[1]):
        if uniform.rvs() < (1 - pd):
           xtmp [:,i] = 0
    pass_columns = np.all(xtmp == 0, axis=0)
    return xtmp, pass_columns


def make_true (tg1,n): 
    # n = 125 # количество измерении 
    # n = np.pi*r
    x1=[];y1=[];z1=[]; vx1=[];vy1=[];vz1=[]; w=[]  
    for i in range(n):
        state1 = tg1.CT(dt)
        x1.append(state1['x'])
        y1.append(state1['y'])
        z1.append(state1['z'])
        vx1.append(state1['vx'])
        vy1.append(state1['vy'])
        vz1.append(state1['vz'])
        w.append(state1['w'])

    X_true_data_not_pass = np.array([x1,vx1,y1,vy1,z1,vz1,w])              
    # plt.plot(x1,y1,'b')
    # plt.grid(True)
    return (X_true_data_not_pass)
n=125
X_true_data_not_pass = make_true(tg1,n)

X_true_data_with_pass, pass_index = make_pass(X_true_data_not_pass, pd)


with_pass = remove_zero_columns(X_true_data_with_pass)


# # ================= Блок 2 =================
# # создание истинной зашумленной траектории
process_var = 100
Qp = np.diag([process_var, process_var, 1.0, 1.0 *(np.pi/180.0)]) 
G = np.array([[dt**2/2,  0.0,         0.0, 0.0],
             [dt,       0.0,          0.0, 0.0],
             [0.0,      dt**2/2,      0.0, 0.0],
             [0.0,      dt,           0.0, 0.0],
             [0.0,      0.0,      dt**2/2, 0.0],
             [0.0,      0.0,      dt     , 0.0],
             [0.0,      0.0,      0.0    , 1.0]])

Q = G@Qp@G.T
def add_process_noise(X, Var):
    X_true_plus_ProcNoise = X + np.sqrt(Var) @ np.random.normal(loc=0, scale=1.0, size=(X.shape[0], X.shape[1]))
    return X_true_plus_ProcNoise
X_true_plus_ProcNoise = add_process_noise(X_true_data_not_pass, Q) # движение цели по модели

X_true_plus_ProcNoise_with_pass = add_process_noise(X_true_data_with_pass, Q) # движение цели, но с пропусками

#=================================================================
#построение графика истинной зашумленной траектории с пропусками

X_true_plus_ProcNoise_not_pass = X_true_plus_ProcNoise_with_pass
X_true_plus_ProcNoise_not_pass[:,pass_index] = 0
X_true_plus_ProcNoise_not_pass = X_true_plus_ProcNoise_with_pass[:, ~pass_index] # вырезаю пропуски, что построить график
# print("X_true_plus_ProcNoise_not_pass" , X_true_plus_ProcNoise_not_pass)
plt.plot(X_true_plus_ProcNoise_not_pass[0],X_true_plus_ProcNoise_not_pass[2], label='X_true_plus_ProcNoise+pass', linestyle='-', marker='o')
plt.legend()
# print("X_true+noise", X_true_plus_ProcNoise)


# ================= Блок 3 =====================

def Zsph2cart(Z):
    
    x = Z[0] * np.cos(Z[1]) * np.cos(Z[2])
    y = Z[0] * np.sin(Z[1]) * np.cos(Z[2])
    z = Z[0] * np.sin(Z[2])
    Z_cart = np.vstack((x,y,z))
    return Z_cart

def do_measurement(X_plusProcNoise,R, pass_index):

    X_plusProcNoise[:, pass_index] = 0
    r_true_with_noise = np.sqrt(np.array(X_plusProcNoise[0])**2 + np.array(X_plusProcNoise[2])**2 + np.array(X_plusProcNoise[4])**2)
    az_true_with_noise = np.arctan2(np.array(X_plusProcNoise[2]),np.array(X_plusProcNoise[0])) # Азимут
    um_true_with_noise = np.arctan2(np.array(X_plusProcNoise[4]),np.sqrt(np.array(X_plusProcNoise[0])**2+np.array(X_plusProcNoise[2])**2)) # Угол места

    vr_with_noise = (np.array(X_plusProcNoise[0])*np.array(X_plusProcNoise[1]) + np.array(X_plusProcNoise[2])*np.array(X_plusProcNoise[3]) + np.array(X_plusProcNoise[4])* np.array(X_plusProcNoise[5])) / r_true_with_noise 
    vr_with_noise = np.nan_to_num(vr_with_noise, nan=0.)

    Zm = np.zeros((R.shape[0], X_plusProcNoise.shape[1]))
    for i in range(Zm.shape[1]):
        Zm[0,i] = r_true_with_noise[i]
        Zm[1,i] = az_true_with_noise[i]
        Zm[2,i] = um_true_with_noise[i]
        # Zm[3,i] = vr_with_noise[i]

    Z_plus_noise = Zm + np.sqrt(R) @ np.random.normal(loc=0, scale=math.sqrt(1.0), size=(Zm.shape[0], Zm.shape[1]))
    Z_plus_noise[:, pass_index] = 0

    return Z_plus_noise

Z = do_measurement (X_true_plus_ProcNoise_with_pass, R, pass_index)
# #==================Отрисовка==================
Z_not_pass = Z
Z_not_pass = remove_zero_columns(Z_not_pass)
Zc = Zsph2cart(Z_not_pass)
plt.plot(Zc[0],Zc[1], label='do_meas', linestyle='-', marker='x')
plt.legend()

# # ================= Блок 4 ===================


def fx(x, dt):

    w = x[6]
    if w == 0.0:
        w= 0.0000001
    F = np.array([[1.0,  1/w*np.sin(w*dt),     0.0,      -1/w*(1-np.cos(w*dt)),  0.0,    0.0,     0.0],
             [0.0,  np.cos(w*dt),         0.0,      -np.sin(w*dt),          0.0,    0.0,     0.0],
             [0.0,  1/w*(1-np.cos(w*dt)), 1.0,      1/w*np.sin(w*dt),       0.0,    0.0,     0.0],
             [0.0,  np.sin(w*dt),         0.0,       np.cos(w*dt),          0.0,    0.0,     0.0],
             [0.0,  0.0,                  0.0,      0.0,                    1.0,     dt,     0.0],
             [0.0,  0.0,                  0.0,      0.0,                    0.0,    1.0,     0.0],
             [0.0,  0.0,                  0.0,      0.0,                    0.0,    0.0,     1.0]], dtype=float)
    return np.dot(F, x)


def hx(x):

    range = np.sqrt(np.power(x[0], 2) + pow(x[2], 2) + pow(x[4], 2))
    az = np.arctan2(x[2], x[0])
    el = np.arctan2(x[4], np.sqrt(np.power(x[0], 2) + np.power(x[2], 2)))

    return np.array([range, az, el])

points = MerweScaledSigmaPoints(7, alpha=.1, beta=2., kappa=-1)

Z_cart = Zsph2cart(Z)

w1 = 0.0000000
kf = UnscentedKalmanFilter(dim_x=7, dim_z=3, dt=dt, fx=fx,hx=hx, points=points)
kf.x = np.array(([Z_cart[0,0], 0., Z_cart[1,0], 0., Z_cart[2,0], 0.,w1]))

kf.P = 100.0*(np.pi/180.0)
z_std = 0.1 *(np.pi/180.0)
kf.R = R
kf.Q = Q#_discrete_white_noise(dim=4, dt=dt, var=0.5**2, block_size=3)


def estimate (Z):
    
    X_c = np.array([])
    for i in range (Z.shape[1]-1):
        kf.predict()
        kf.update(Z[:,i])
        print("kf.x", kf.x)
        if X_c.size == 0:
            X_c = kf.x
        X_c = np.vstack((X_c,kf.x))
    return X_c 



X_c = estimate(Z)

# print(X_c)

# def err1(X_c,X_true_plus_ProcNoise):

#     er = X_c[:,:] - X_true_plus_ProcNoise [:,1:]
#     return er

# e = err1(X_c, X_true_plus_ProcNoise)
# # print("e =", e[0]) #Распечатка ошибок по Х одной трассы


# #==================Отрисовка==================
Z_cart = Zsph2cart(Z)
plt.figure()
plt.plot(X_c[:,0], X_c[:,2], label='Correct', marker='o')
plt.plot(X_true_plus_ProcNoise[0],X_true_plus_ProcNoise[2], label='truth', marker='x')
plt.plot(Zc[0], Zc[1], label='Meas',marker='o')
plt.legend()


# # ================= Блок 5 ===================
# # СБОР СТАТИСТИКИ

def calc_err(X):

    Xn = add_process_noise(X, Q)
    X_pass, pass_id = make_pass(Xn,pd)
    Zn = do_measurement(X_pass, R, pass_id)
    X_c = estimate(Zn)

    print("x_c",X_c)
    print("x_n",Xn)
    err = X_c[:,:] - Xn [:,:].T # ошибка вычисляется со второго столбца.

    # print("ошибка в статистике calc_err",err[0])

    return err

from tqdm import tqdm

def calc_std_err(X):
    num_iterations = 2
    var_err = np.zeros((X.shape[0], X.shape[1])).T

    for i in tqdm(range(num_iterations)):
        err = calc_err(X)
        var_err += err ** 2

    var_err /= num_iterations
    return np.sqrt(var_err)

# tg2G = Target()
# tg2G.init_state({'x':0.0,'y':0.0,'z':0.0, 'vx':200.0,'vy':0.0,'vz':0.0,'w':0.098})
# n=125
# X_true_data_not_pass_2G = make_true(tg2G,n)
# w = 0.00000001
# std_err_2G = calc_std_err(X_true_data_not_pass_2G)

# tg5G = Target()
# tg5G.init_state({'x':0.0,'y':0.0,'z':0.0, 'vx':200.0,'vy':0.0,'vz':0.0,'w':0.245})
# n=50
# X_true_data_not_pass_5G = make_true(tg5G,n)
# w = 0.245
# std_err_5G = calc_std_err(X_true_data_not_pass_5G, w)

# tg8G = Target()
# tg8G.init_state({'x':0.0,'y':0.0,'z':0.0, 'vx':200.0,'vy':0.0,'vz':0.0,'w':0.392})
# n=31
# X_true_data_not_pass_8G = make_true(tg8G,n)
# w = 0.392
# std_err_8G = calc_std_err(X_true_data_not_pass_8G,w)


# plt.figure(num="2G")
# plt.subplot(6, 1, 1)
# plt.plot((np.arange(len(std_err_2G[0, :]))+1)*dt, std_err_2G[0, :])
# plt.xlabel('Time,s')
# plt.ylabel('std_x, met')
# plt.grid(True)
# plt.subplot(6, 1, 2)
# plt.plot((np.arange(len(std_err_2G[1, :]))+1)*dt, std_err_2G[1, :])
# plt.grid(True)
# plt.xlabel('Time,s')
# plt.ylabel('std_vx, m/s')
# plt.subplot(6, 1, 3)
# plt.plot((np.arange(len(std_err_2G[2, :]))+1)*dt, std_err_2G[2, :])
# plt.grid(True)
# plt.xlabel('Time,s')
# plt.ylabel('std_y, met')
# plt.subplot(6, 1, 4)
# plt.plot((np.arange(len(std_err_2G[3, :]))+1)*dt, std_err_2G[3, :])
# plt.grid(True)
# plt.xlabel('Time,s')
# plt.ylabel('std_vy, m/s')
# plt.subplot(6, 1, 5)
# plt.plot((np.arange(len(std_err_2G[4, :]))+1)*dt, std_err_2G[4, :])
# plt.grid(True)
# plt.xlabel('Time,s')
# plt.ylabel('std_z, met')
# plt.subplot(6, 1, 6)
# plt.plot((np.arange(len(std_err_2G[5, :]))+1)*dt, std_err_2G[5, :])
# plt.grid(True)
# plt.xlabel('Time,s')
# plt.ylabel('std_vz, m/s')
# plt.subplot(7, 1, 7)
# plt.plot((np.arange(len(std_err_2G[6, :]))+1)*dt, std_err_2G[6, :])
# plt.xlabel('Time,s')
# plt.ylabel('std_w, rad')
# plt.grid(True)
# plt.subplots_adjust(wspace=12.0, hspace=1.0)


# plt.figure(num="5G")
# plt.subplot(6, 1, 1)
# plt.plot((np.arange(len(std_err_5G[0, 0:50])))*dt, std_err_5G[0, 0:50])
# plt.xlabel('Time,s')
# plt.ylabel('std_x, met')
# plt.grid(True)
# plt.subplot(6, 1, 2)
# plt.plot((np.arange(len(std_err_5G[1, 0:50])))*dt, std_err_5G[1, 0:50])
# plt.grid(True)
# plt.xlabel('Time,s')
# plt.ylabel('std_vx, m/s')
# plt.subplot(6, 1, 3)
# plt.plot((np.arange(len(std_err_5G[2, 0:50])))*dt, std_err_5G[2, 0:50])
# plt.grid(True)
# plt.xlabel('Time,s')
# plt.ylabel('std_y, met')
# plt.subplot(6, 1, 4)
# plt.plot((np.arange(len(std_err_5G[3, 0:50])))*dt, std_err_5G[3, 0:50])
# plt.grid(True)
# plt.xlabel('Time,s')
# plt.ylabel('std_vy, m/s')
# plt.subplot(6, 1, 5)
# plt.plot((np.arange(len(std_err_5G[4, 0:50])))*dt, std_err_5G[4, 0:50])
# plt.grid(True)
# plt.xlabel('Time,s')
# plt.ylabel('std_z, met')
# plt.subplot(6, 1, 6)
# plt.plot((np.arange(len(std_err_5G[5, 0:50])))*dt, std_err_5G[5, 0:50])
# plt.grid(True)
# plt.xlabel('Time,s')
# plt.ylabel('std_vz, m/s')
# plt.subplots_adjust(wspace=12.0, hspace=1.0)

# plt.figure(num="8G")
# plt.subplot(6, 1, 1)
# plt.plot((np.arange(len(std_err_8G[0, 0:31])))*dt, std_err_8G[0, 0:31])
# plt.xlabel('Time,s')
# plt.ylabel('std_x, met')
# plt.grid(True)
# plt.subplot(6, 1, 2)
# plt.plot((np.arange(len(std_err_8G[1, 0:31])))*dt, std_err_8G[1, 0:31])
# plt.grid(True)
# plt.xlabel('Time,s')
# plt.ylabel('std_vx, m/s')
# plt.subplot(6, 1, 3)
# plt.plot((np.arange(len(std_err_8G[2, 0:31])))*dt, std_err_8G[2, 0:31])
# plt.grid(True)
# plt.xlabel('Time,s')
# plt.ylabel('std_y, met')
# plt.subplot(6, 1, 4)
# plt.plot((np.arange(len(std_err_8G[3, 0:31])))*dt, std_err_8G[3, 0:31])
# plt.grid(True)
# plt.xlabel('Time,s')
# plt.ylabel('std_vy, m/s')
# plt.subplot(6, 1, 5)
# plt.plot((np.arange(len(std_err_8G[4, 0:31])))*dt, std_err_8G[4, 0:31])
# plt.grid(True)
# plt.xlabel('Time,s')
# plt.ylabel('std_z, met')
# plt.subplot(6, 1, 6)
# plt.plot((np.arange(len(std_err_8G[5, 0:31])))*dt, std_err_8G[5, 0:31])
# plt.grid(True)
# plt.xlabel('Time,s')
# plt.ylabel('std_vz, m/s')
# plt.subplots_adjust(wspace=12.0, hspace=1.0)

# plt.figure(num="Together")
# plt.subplot(6, 1, 1)
# plt.plot((np.arange(len(std_err_2G[0, :]))+1)*dt, std_err_2G[0, :], label='2G')
# plt.plot((np.arange(len(std_err_5G[0, 0:50])))*dt, std_err_5G[0, 0:50], label='5G')
# plt.plot((np.arange(len(std_err_8G[0, 0:31])))*dt, std_err_8G[0, 0:31], label='8G')
# plt.xlabel('Time,s')
# plt.ylabel('std_x, met')
# plt.grid(True)
# plt.subplot(6, 1, 2)
# plt.plot((np.arange(len(std_err_2G[1, :]))+1)*dt, std_err_2G[1, :], label='2G')
# plt.plot((np.arange(len(std_err_5G[1, 0:50])))*dt, std_err_5G[1, 0:50], label='5G')
# plt.plot((np.arange(len(std_err_8G[1, 0:31])))*dt, std_err_8G[1, 0:31], label='8G')
# plt.grid(True)
# plt.xlabel('Time,s')
# plt.ylabel('std_vx, m/s')
# plt.subplot(6, 1, 3)
# plt.plot((np.arange(len(std_err_2G[2, :]))+1)*dt, std_err_2G[2, :], label='2G')
# plt.plot((np.arange(len(std_err_5G[2, 0:50])))*dt, std_err_5G[2, 0:50], label='5G')
# plt.plot((np.arange(len(std_err_8G[2, 0:31])))*dt, std_err_8G[2, 0:31], label='8G')
# plt.grid(True)
# plt.xlabel('Time,s')
# plt.ylabel('std_y, met')
# plt.subplot(6, 1, 4)
# plt.plot((np.arange(len(std_err_2G[3, :]))+1)*dt, std_err_2G[3, :], label='2G')
# plt.plot((np.arange(len(std_err_5G[3, 0:50])))*dt, std_err_5G[3, 0:50], label='5G')
# plt.plot((np.arange(len(std_err_8G[3, 0:31])))*dt, std_err_8G[3, 0:31], label='8G')
# plt.grid(True)
# plt.xlabel('Time,s')
# plt.ylabel('std_vy, m/s')
# plt.subplot(6, 1, 5)
# plt.plot((np.arange(len(std_err_2G[4, :]))+1)*dt, std_err_2G[4, :], label='2G')
# plt.plot((np.arange(len(std_err_5G[4, 0:50])))*dt, std_err_5G[4, 0:50], label='5G')
# plt.plot((np.arange(len(std_err_8G[4, 0:31])))*dt, std_err_8G[4, 0:31], label='8G')
# plt.grid(True)
# plt.xlabel('Time,s')
# plt.ylabel('std_z, met')
# plt.subplot(6, 1, 6)
# plt.plot((np.arange(len(std_err_2G[5, :]))+1)*dt, std_err_2G[5, :], label='2G')
# plt.plot((np.arange(len(std_err_5G[5, 0:50])))*dt, std_err_5G[5, 0:50], label='5G')
# plt.plot((np.arange(len(std_err_8G[5, 0:31])))*dt, std_err_8G[5, 0:31], label='8G')
# plt.grid(True)
# plt.xlabel('Time,s')
# plt.ylabel('std_vz, m/s')
# plt.subplots_adjust(wspace=12.0, hspace=1.0)
# plt.legend(bbox_to_anchor=(1, 1), loc='upper left')
plt.show()

