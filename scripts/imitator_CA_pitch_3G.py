import time
import numpy as np
import math
import matplotlib.pyplot as plt
import plotly.graph_objs as go
from plotly.subplots import make_subplots
from scipy.stats import chi2, poisson, uniform, norm
from IPython.display import clear_output
import random
import estimator
from models import Target

# dt = 0.25
dt = 1.0
pd = 1.0

R = np.diag([10000.0, (0.1/3)**2,(0.1/3)**2]) #дисперсии, в deg
plt.rcParams['figure.figsize'] = [10, 6]
fig2 = make_subplots(rows=1, cols=1, specs=[[{'type': 'scatter3d'}]])

   

# =============== Блок 1 ===================
    
# ИНИЦИАЛИЗАЦИЯ МОДЕЛИ ДВИЖЕНИЯ
tg1 = Target()

# tg1.init_state({'x':10000.0, 'y':0.0, 'z':10000.0, 'vx':200.0, 'vy':0.0, 'vz':0.0, 'ax': 45.0, 'ay': 0.0, 'az': 29.4}) #for dt = 0.25, n = 50
tg1.init_state({'x':10000.0, 'y':0.0, 'z':10000.0, 'vx':200.0, 'vy':0.0, 'vz':0.0, 'ax': 45.0, 'ay': 0.0, 'az': 29.4})

def remove_zero_columns(arr):

    zero_columns = np.all(arr == 0, axis=0)
        
    # Удаляем столбцы, заполненные нулями
    arr = arr[:, ~zero_columns]
        
    return arr

def make_pass(X_true_data, pd): #при pd = 1 пропусков не будет
    xtmp = X_true_data.copy()
    
    for i in range (2, xtmp.shape[1]):
        if uniform.rvs() < (1 - pd):
           xtmp [:,i] = 0
    pass_columns = np.all(xtmp == 0, axis=0)
    return xtmp, pass_columns

def make_true (tg1,n): 

    x1=[];  y1=[]; z1=[]; vx1=[]; vy1=[]; vz1=[]; ax1 =[];  ay1 =[];   az1 =[]
    for i in range(n):
        state1 = tg1.CA(dt)
        x1.append(state1['x'])
        vx1.append(state1['vx'])
        ax1.append(state1['ax'])

        y1.append(state1['y'])
        vy1.append(state1['vy'])
        ay1.append(state1['ay'])

        z1.append(state1['z'])
        vz1.append(state1['vz'])
        az1.append(state1['az'])
        

    X_true_data_not_pass = np.array([x1,vx1,ax1,y1,vy1,ay1,z1,vz1,az1])             
    # plt.plot(x1,y1,'b')
    # plt.grid(True)
    # plt.show()
    return (X_true_data_not_pass)

n = 18
# n = 50 # for dt = 0.25
X_true_data_not_pass = make_true(tg1,n)

# print("x_true",X_true_data_not_pass)

X_true_data_with_pass, pass_index = make_pass(X_true_data_not_pass, pd)


with_pass = remove_zero_columns(X_true_data_with_pass)
# print("with_pass",with_pass)
# Создаем трехмерный scatter plot для массива
# scatter2 = go.Scatter3d(x=with_pass[0], y=with_pass[3], z=with_pass[6], mode='markers+lines', marker=dict(size=3, color='blue'), name='X_true_data')
# fig2.add_trace(scatter2)

# # Обновляем параметры макета
# fig2.update_layout(scene=dict(aspectmode="cube"))
# fig2.show()


# # ================= Блок 2 =================
# # создание истинной зашумленной траектории
process_var = 10
Qp = np.diag([process_var, process_var, process_var]) 
G = np.array([[dt**2/2,  0.0,          0.0],
             [dt,        0.0,          0.0],
             [1.0,       0.0,          0.0],
             [0.0,       dt**2/2,      0.0],
             [0.0,       dt,           0.0],
             [0.0,       1.0,          0.0],
             [0.0,       0.0,      dt**2/2],
             [0.0,       0.0,      dt     ],
             [0.0,       0.0,      1.0    ]])

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
plt.plot(X_true_plus_ProcNoise_not_pass[0],X_true_plus_ProcNoise_not_pass[6], label='X_true_plus_ProcNoise+pass', linestyle='-', marker='o')
plt.legend()
# print("X_true+ procces_noise", X_true_plus_ProcNoise)


# # ================= Блок 3 =====================

def Zsph2cart(Z):
    
    x = Z[0] * np.cos(np.deg2rad(Z[1])) * np.cos(np.deg2rad(Z[2]))
    y = Z[0] * np.sin(np.deg2rad(Z[1])) * np.cos(np.deg2rad(Z[2]))
    z = Z[0] * np.sin(np.deg2rad(Z[2]))
    Z_cart = np.vstack((x,y,z))
    return Z_cart

def do_measurement(X_plusProcNoise,R, pass_index):

    X_plusProcNoise[:, pass_index] = 0
    r_true_with_noise = np.sqrt(np.array(X_plusProcNoise[0])**2 + np.array(X_plusProcNoise[3])**2 + np.array(X_plusProcNoise[6])**2)
    az_true_with_noise = np.rad2deg(np.arctan2(np.array(X_plusProcNoise[3]),np.array(X_plusProcNoise[0]))) # Азимут
    um_true_with_noise = np.rad2deg(np.arctan2(np.array(X_plusProcNoise[6]),np.sqrt(np.array(X_plusProcNoise[0])**2+np.array(X_plusProcNoise[3])**2))) # Угол места

    #Wrong [index]
    # vr_with_noise = (np.array(X_plusProcNoise[0])*np.array(X_plusProcNoise[1]) + np.array(X_plusProcNoise[2])*np.array(X_plusProcNoise[3]) + np.array(X_plusProcNoise[4])* np.array(X_plusProcNoise[5])) / r_true_with_noise 
    # vr_with_noise = np.nan_to_num(vr_with_noise, nan=0.)

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
# print('Z=', Z)
dictData = {}
dictData['DeltaTime'] = dt
dictData['Measurement'] = Z.tolist()
dictData['MeasurementNoise'] = R.tolist()
dictData['ProcessNoise'] = Qp.tolist()

import json

msg_json = json.dumps(dictData)
with open("turnData.json", "w") as json_file:
    json_file.write(msg_json)


Z_not_pass = remove_zero_columns(Z_not_pass)
Zc = Zsph2cart(Z_not_pass)
# print('Zc=', Zc)
plt.plot(Zc[0],Zc[2], label='do_meas', linestyle='-', marker='x')
plt.legend()

# # ================= Блок 4 ===================

def estimate (Z):
    
    Z_cart = Zsph2cart(Z)
    X0 = np.vstack([Z_cart[0,0], 0.0, 0.0, Z_cart[1,0], 0.0, 0.0, Z_cart[2,0], 0.0, 0.0]) # инициализируем вектор состояния, равный первому измерению
    # print('X0=',X0)
    point = estimator.Points()
    point.alpha = 1e-3
    point.beta = 2
    point.kappa = 3 - X0.shape[0]
    ukf = estimator.BindTrackUkf_CA(X0,dt,Qp,R,point) #инициал. фильтра
    X_c = np.empty((len(X0), 0))
    for i in range (Z.shape[1]-1):
        if np.all(Z[:,i+1] == 0):
            X = ukf.step(dt)
            X_c = np.append(X_c,X,axis=1)
            continue
        # print('Z=',Z[:,i+1])
        X = ukf.step(Z[:,i+1])
        # print('X=',X)
        X_c = np.append(X_c,X,axis=1)
    print("X_Estimeted=",X_c)
    return X_c 

X_c = estimate(Z)
# print("X_Estimeted=",X_c)

# def err1(X_c,X_true_plus_ProcNoise):

#     er = X_c[:,:] - X_true_plus_ProcNoise [:,1:]
#     return er

# e = err1(X_c, X_true_plus_ProcNoise)
# # print("e =", e[0]) #Распечатка ошибок по Х одной трассы


# #==================Отрисовка==================
# Z_cart = Zsph2cart(Z)
plt.figure()
plt.plot(X_c[0], X_c[6], label='Correct', marker='o')
plt.plot(X_true_plus_ProcNoise[0],X_true_plus_ProcNoise[6], label='truth', marker='x')
plt.plot(Zc[0], Zc[2], label='Meas',marker='o')
plt.legend()


# # # ================= Блок 5 ===================
# # СБОР СТАТИСТИКИ

def calc_err(X):

    Xn = add_process_noise(X, Q)
    X_pass, pass_id = make_pass(Xn,pd)
    Zn = do_measurement(X_pass, R, pass_id)
    X_c = estimate(Zn)

    err = X_c[:,:] - Xn [:,1:] # ошибка вычисляется со второго столбца.

    # print("ошибка в статистике calc_err",err[0])

    return err

from tqdm import tqdm

def calc_std_err(X):
    num_iterations = 1
    var_err = np.zeros((X.shape[0], X.shape[1]-1))

    for i in tqdm(range(num_iterations)):
        err = calc_err(X)
        var_err += err ** 2

    var_err /= num_iterations
    return np.sqrt(var_err)


tg3G = Target()
tg3G.init_state({'x':10000.0, 'y':0.0, 'z':10000.0, 'vx':200.0, 'vy':0.0, 'vz':0.0, 'ax': 45.0, 'ay': 0.0, 'az': 29.4})

n=18 # for dt = 1
# n=50 # for dt = 0.25
X_true_data_not_pass_3G = make_true(tg3G,n)
std_err_3G = calc_std_err(X_true_data_not_pass_3G)



plt.figure(num="3G")
plt.subplot(9, 1, 1)
plt.plot((np.arange(len(std_err_3G[0, :])))*dt, std_err_3G[0, :])
plt.xlabel('Time,s')
plt.ylabel('std_x, met')
plt.grid(True)
plt.subplot(9, 1, 2)
plt.plot((np.arange(len(std_err_3G[1, :])))*dt, std_err_3G[1,:])
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_vx, m/s')

plt.subplot(9, 1, 3)
plt.plot((np.arange(len(std_err_3G[2, :])))*dt, std_err_3G[2, :])
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_ax, m/s^2')

plt.subplot(9, 1, 4)
plt.plot((np.arange(len(std_err_3G[3, :])))*dt, std_err_3G[3, :])
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_y, met')

plt.subplot(9, 1, 5)
plt.plot((np.arange(len(std_err_3G[4,:])))*dt, std_err_3G[4, :])
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_vy, m/s')

plt.subplot(9, 1, 6)
plt.plot((np.arange(len(std_err_3G[5, :])))*dt, std_err_3G[5, :])
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_ay, m/s^2')

plt.subplot(9, 1, 7)
plt.plot((np.arange(len(std_err_3G[6, :])))*dt, std_err_3G[6, :])
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_z, m')

plt.subplot(9, 1, 8)
plt.plot((np.arange(len(std_err_3G[7, :])))*dt, std_err_3G[7, :])
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_vx, m/s')

plt.subplot(9, 1, 9)
plt.plot((np.arange(len(std_err_3G[8, :])))*dt, std_err_3G[8, :])
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_az, m/s^2')

plt.grid(True)
plt.subplots_adjust(wspace=12.0, hspace=1.0)



plt.show()

