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

dt = 6.0
pd = 1.0

R_with_vr = np.diag([10000.0, (0.1/3.0)**2, (0.1/3.0)**2, 25.0])
R_without_vr = np.diag([10000.0, (0.1/3.0)**2, (0.1/3.0)**2])

plt.rcParams['figure.figsize'] = [10, 6]
fig2 = make_subplots(rows=1, cols=1, specs=[[{'type': 'scatter3d'}]])

    
# =============== Блок 1 ===================
    
# ИНИЦИАЛИЗАЦИЯ МОДЕЛИ ДВИЖЕНИЯ
tg1 = Target()
tg1.init_state({'x':50000.0,'y':0.0,'z':0.0, 'vx':-200.0,'vy':0.0,'vz':0.0})
# tg1.init_state({'x':0.0,'y':0.0,'z':0.0, 'vx':200.0,'vy':0.0,'vz':0.0})

def remove_zero_columns(arr):

    zero_columns = np.all(arr == 0, axis=0)
    arr = arr[:, ~zero_columns]
        
    return arr

def make_pass(X_true_data, pd):
    xtmp = X_true_data.copy()
    
    for i in range (2,xtmp.shape[1]):
        if uniform.rvs() < (1 - pd):
           xtmp [:,i] = 0
    pass_columns = np.all(xtmp == 0, axis=0)
    return xtmp, pass_columns


def make_true (tg1): 
    n = 90 # количество измерении

    x1=[];y1=[];z1=[]; vx1=[];vy1=[];vz1=[]; w=[]  
    for i in range(n):
        state1 = tg1.CV(dt)
        x1.append(state1['x'])
        y1.append(state1['y'])
        z1.append(state1['z'])
        vx1.append(state1['vx'])
        vy1.append(state1['vy'])
        vz1.append(state1['vz'])

    X_true_data_not_pass = np.array([x1,vx1,y1,vy1,z1,vz1])             
    # plt.plot(x1,y1,'b')
    # plt.grid(True)
    return (X_true_data_not_pass)

X_true_data_not_pass = make_true(tg1)
X_true_data_with_pass, pass_index = make_pass(X_true_data_not_pass, pd)

#==================Отрисовка======================
with_pass = remove_zero_columns(X_true_data_with_pass)
# Создаем трехмерный scatter plot для массива
# scatter2 = go.Scatter3d(x=with_pass[0], y=with_pass[2], z=with_pass[4], mode='markers+lines', marker=dict(size=3, color='blue'), name='X_true_data')
# fig2.add_trace(scatter2)

# # Обновляем параметры макета
# fig2.update_layout(scene=dict(aspectmode="cube"))
# fig2.show()


# ================= Блок 2 =================
# создание истинной зашумленной траектории
process_var = 0.0001
Qp = np.diag([process_var, process_var, process_var])
G = np.array([[dt**2/2,  0.0,         0.0],
             [dt,       0.0,          0.0],
             [0.0,      dt**2/2,      0.0],
             [0.0,      dt,           0.0],
             [0.0,      0.0,      dt**2/2],
             [0.0,      0.0,      dt     ]])
Q = G@Qp@G.T


def add_process_noise(X, Var):
    X_true_plus_ProcNoise = X + np.sqrt(Var) @ np.random.normal(loc=0, scale=1.0, size=(X.shape[0], X.shape[1]))
    return X_true_plus_ProcNoise

X_true_plus_ProcNoise = add_process_noise(X_true_data_not_pass, Q) # движение цели по модели
X_true_plus_ProcNoise_with_pass = add_process_noise(X_true_data_with_pass, Q) # движение цели, но с пропусками


# ================= Блок 3 =====================

def Zsph2cart(Z):
    
    # x = Z[0] * np.cos(Z[1]) * np.cos(Z[2])
    # y = Z[0] * np.sin(Z[1]) * np.cos(Z[2])
    # z = Z[0] * np.sin(Z[2])
    # Z_cart = np.vstack((x,y,z))
    # return Z_cart
    x = Z[0] * np.cos(np.deg2rad(Z[1])) * np.cos(np.deg2rad(Z[2]))
    y = Z[0] * np.sin(np.deg2rad(Z[1])) * np.cos(np.deg2rad(Z[2]))
    z = Z[0] * np.sin(np.deg2rad(Z[2]))
    Z_cart = np.vstack((x,y,z))
    return Z_cart

def do_measurement(X_plusProcNoise,R_with_vr, pass_index):

    X_plusProcNoise[:, pass_index] = 0
    r_true_with_noise = np.sqrt(np.array(X_plusProcNoise[0])**2 + np.array(X_plusProcNoise[2])**2 + np.array(X_plusProcNoise[4])**2)
    az_true_with_noise = np.rad2deg(np.arctan2(np.array(X_plusProcNoise[2]),np.array(X_plusProcNoise[0]))) # Азимут
    um_true_with_noise = np.rad2deg(np.arctan2(np.array(X_plusProcNoise[4]),np.sqrt(np.array(X_plusProcNoise[0])**2+np.array(X_plusProcNoise[2])**2))) # Угол места
  
    vr_with_noise = (np.array(X_plusProcNoise[0])*np.array(X_plusProcNoise[1]) + np.array(X_plusProcNoise[2])*np.array(X_plusProcNoise[3]) + np.array(X_plusProcNoise[4])* np.array(X_plusProcNoise[5])) / r_true_with_noise 
    vr_with_noise = np.nan_to_num(vr_with_noise, nan=0.)

    Zm = np.zeros((R_with_vr.shape[0], X_plusProcNoise.shape[1]))
    for i in range(Zm.shape[1]):
        Zm[0,i] = r_true_with_noise[i]
        Zm[1,i] = az_true_with_noise[i]
        Zm[2,i] = um_true_with_noise[i]
        Zm[3,i] = vr_with_noise[i]

    Z_plus_noise_vr = Zm + np.sqrt(R_with_vr) @ np.random.normal(loc=0, scale=math.sqrt(1.0), size=(Zm.shape[0], Zm.shape[1]))
    Z_plus_noise_vr[:, pass_index] = 0

    Z_plus_noise_without_vr = Z_plus_noise_vr[:-1]

    return Z_plus_noise_without_vr, Z_plus_noise_vr

Z, Zvr = do_measurement (X_true_plus_ProcNoise_with_pass, R_with_vr, pass_index)


Ztmp = remove_zero_columns(Z)
Zc = Zsph2cart(Ztmp)

print ("Zc",Zc)
# ================= Блок 4 ===================
def estimate (Z):

    Z_cart = Zsph2cart(Z)
    X0 = np.vstack([Z_cart[0,0], 0., Z_cart[1,0], 0., Z_cart[2,0], 0.]) # инициализируем вектор состояния, равный первому измерению
    # print ("X0", X0)

    point = estimator.Points()
    point.alpha = 1e-3
    point.beta = 2
    point.kappa = 3 - X0.shape[0]

    ukf = estimator.BindTrackUkf_CV(X0,Qp,R_without_vr,point) #инициал. фильтра
    
    X_c = np.empty((len(X0), 0))
    for i in range (Z.shape[1]-1):
        if np.all(Z[:,i+1] == 0):
            X = ukf.step(dt)
            X_c = np.append(X_c,X,axis=1)
            continue
        X = ukf.step(dt, Z[:,i+1])
        X_c = np.append(X_c,X,axis=1)
    return X_c 

def estimate_with_vr (Zvr):

    Z_cart = Zsph2cart(Zvr)
    X0 = np.vstack([Z_cart[0,0], 0., Z_cart[1,0], 0., Z_cart[2,0], 0.]) # инициализируем вектор состояния, равный первому измерению
    # print ("X0", X0)
    point = estimator.Points()
    point.alpha = 1e-3
    point.beta = 2
    point.kappa = 3 - X0.shape[0]
    ukf = estimator.BindTrackUkf_CV(X0,Qp,R_with_vr,point) #инициал. фильтра
    
    X_c = np.empty((len(X0), 0))
    for i in range (Zvr.shape[1]-1):
        if np.all(Zvr[:,i+1] == 0):
            X = ukf.step(dt)
            X_c = np.append(X_c,X,axis=1)
            continue
        X = ukf.step(dt, Zvr[:,i+1])
        X_c = np.append(X_c,X,axis=1)
    return X_c 


X_c = estimate(Z)
X_c_with_Vr = estimate_with_vr(Zvr)
# print ("X_c",X_c)
# print ("X_c_with_Vr",X_c_with_Vr)
def err1(X_c,X_true_plus_ProcNoise):

    er = X_c[:,:] - X_true_plus_ProcNoise [:,1:]
    return er

e = err1(X_c, X_true_plus_ProcNoise)
# print("e =", e[0]) #Распечатка ошибок по Х одной трассы


#==================Отрисовка==================
# Z_cart = Zsph2cart(Z)
plt.figure()
plt.plot(X_c[0], X_c[2], label='Correct', marker='o')
plt.plot(X_c_with_Vr[0], X_c_with_Vr[2], label='CorrectVr', marker='o', color='r')
plt.plot(X_true_plus_ProcNoise[0],X_true_plus_ProcNoise[2], label='truth', marker='x')
plt.plot(Zc[0], Zc[1], label='Meas',marker='o')
plt.legend()


# ================= Блок 5 ===================
# # СБОР СТАТИСТИКИ

def calc_err(X):

    Xn = add_process_noise(X, Q)
    X_pass, pass_id = make_pass(Xn,pd)
    Zn,Zvr = do_measurement(X_pass, R_with_vr, pass_id)
    X_c = estimate(Zn)
    X_c_with_Vr = estimate_with_vr(Zvr)

    err = X_c[:,:] - Xn [:,1:] # ошибка вычисляется со второго столбца.
    err_with_vr = X_c_with_Vr[:,:] - Xn [:,1:]
    # print("ошибка в статистике calc_err",err[0])

    return err, err_with_vr

from tqdm import tqdm

def calc_std_err(X):
    num_iterations = 1
    var_err = np.zeros((X.shape[0], X.shape[1]-1))
    var_err_with_vr = np.zeros((X.shape[0], X.shape[1]-1))

    for i in tqdm(range(num_iterations)):
        err,err_with_vr = calc_err(X)
        var_err += err ** 2
        var_err_with_vr += err_with_vr**2

    var_err /= num_iterations
    var_err_with_vr /= num_iterations
    return np.sqrt(var_err), np.sqrt(var_err_with_vr)

std_err, std_err_with_vr = calc_std_err(X_true_data_not_pass)

plt.figure(num="not Vr")
plt.subplot(6, 1, 1)
plt.plot((np.arange(len(std_err[0, :]))+2)*dt, std_err[0, :])
plt.xlabel('Time,s')
plt.ylabel('std_x, met')
plt.grid(True)
plt.subplot(6, 1, 2)
plt.plot((np.arange(len(std_err[1, :]))+2)*dt, std_err[1, :])
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_vx, m/s')
plt.subplot(6, 1, 3)
plt.plot((np.arange(len(std_err[2, :]))+2)*dt, std_err[2, :])
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_y, met')
plt.subplot(6, 1, 4)
plt.plot((np.arange(len(std_err[3, :]))+2)*dt, std_err[3, :])
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_vy, m/s')
plt.subplot(6, 1, 5)
plt.plot((np.arange(len(std_err[4, :]))+2)*dt, std_err[4, :])
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_z, met')
plt.subplot(6, 1, 6)
plt.plot((np.arange(len(std_err[5, :]))+2)*dt, std_err[5, :])
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_vz, m/s')
plt.subplots_adjust(wspace=12.0, hspace=1.0)


plt.figure(num="with Vr")
plt.subplot(6, 1, 1)
plt.plot((np.arange(len(std_err_with_vr[0, :]))+2)*dt, std_err_with_vr[0, :])
plt.xlabel('Time,s')
plt.ylabel('std_x, met')
plt.grid(True)
plt.subplot(6, 1, 2)
plt.plot((np.arange(len(std_err_with_vr[1, :]))+2)*dt, std_err_with_vr[1, :])
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_vx, m/s')
plt.subplot(6, 1, 3)
plt.plot((np.arange(len(std_err_with_vr[2, :]))+2)*dt, std_err_with_vr[2, :])
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_y, met')
plt.subplot(6, 1, 4)
plt.plot((np.arange(len(std_err_with_vr[3, :]))+2)*dt, std_err_with_vr[3, :])
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_vy, m/s')
plt.subplot(6, 1, 5)
plt.plot((np.arange(len(std_err_with_vr[4, :]))+2)*dt, std_err_with_vr[4, :])
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_z, met')
plt.subplot(6, 1, 6)
plt.plot((np.arange(len(std_err_with_vr[5, :]))+2)*dt, std_err_with_vr[5, :])
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_vz, m/s')
plt.subplots_adjust(wspace=12.0, hspace=1.0)

plt.figure(num="together")
plt.subplot(6, 1, 1)
plt.plot((np.arange(len(std_err[0, :]))+2)*dt, std_err[0, :], label='not Vr')
plt.plot((np.arange(len(std_err_with_vr[0, :]))+2)*dt, std_err_with_vr[0, :], label='Vr')
plt.xlabel('Time,s')
plt.ylabel('std_x, met')
plt.grid(True)
plt.subplot(6, 1, 2)
plt.plot((np.arange(len(std_err[1, :]))+2)*dt, std_err[1, :], label='not Vr')
plt.plot((np.arange(len(std_err_with_vr[1, :]))+2)*dt, std_err_with_vr[1, :], label='Vr')
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_vx, m/s')
plt.subplot(6, 1, 3)
plt.plot((np.arange(len(std_err[2, :]))+2)*dt, std_err[2, :], label='not Vr')
plt.plot((np.arange(len(std_err_with_vr[2, :]))+2)*dt, std_err_with_vr[2, :], label='Vr')
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_y, met')
plt.subplot(6, 1, 4)
plt.plot((np.arange(len(std_err[3, :]))+2)*dt, std_err[3, :], label='not Vr')
plt.plot((np.arange(len(std_err_with_vr[3, :]))+2)*dt, std_err_with_vr[3, :], label='Vr')
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_vy, m/s')
plt.subplot(6, 1, 5)
plt.plot((np.arange(len(std_err[4, :]))+2)*dt, std_err[4, :], label='not Vr')
plt.plot((np.arange(len(std_err_with_vr[4, :]))+2)*dt, std_err_with_vr[4, :], label='Vr')
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_z, met')
plt.subplot(6, 1, 6)
plt.plot((np.arange(len(std_err[5, :]))+2)*dt, std_err[5, :], label='not Vr')
plt.plot((np.arange(len(std_err_with_vr[5, :]))+2)*dt, std_err_with_vr[5, :], label='Vr')
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_vz, m/s')


plt.subplots_adjust(wspace=12.0, hspace=0.5)
plt.legend(bbox_to_anchor=(1, 1), loc='upper left')
plt.show()

