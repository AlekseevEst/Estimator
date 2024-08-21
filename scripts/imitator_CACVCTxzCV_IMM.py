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

dt = 1.0
pd = 1.0

R = np.diag([10000.0, (0.1/3)**2,(0.1/3)**2]) #дисперсии, в deg
plt.rcParams['figure.figsize'] = [10, 6]
fig2 = make_subplots(rows=1, cols=1, specs=[[{'type': 'scatter3d'}]])


# =============== Блок 1 ===================
    
# ИНИЦИАЛИЗАЦИЯ МОДЕЛИ ДВИЖЕНИЯ
tg1 = Target()
# init_state = {'x':100000.0, 'y':20000.0, 'z':10000.0, 'vx':-200.0, 'vy':0.0, 'vz':0.0, 'ax': 0.0, 'ay': 0.0, 'az': 19.4}
init_state = {'x':200000, 'y':20000, 'z':20000, 'vx':200,'vy':0.0,'vz':0.0, 'w':5.5}
tg1.init_state(init_state)
nCaPitch = 36
nCv = 300
nCt = 100
nCv2 = 200

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

def make_true (tg1): 

    x1=[];  y1=[]; z1=[]; vx1=[]; vy1=[]; vz1=[]
    # for i in range(nCaPitch):
    #     state1 = tg1.CA(dt)
    #     x1.append(state1['x'])
    #     vx1.append(state1['vx'])

    #     y1.append(state1['y'])
    #     vy1.append(state1['vy'])

    #     z1.append(state1['z'])
    #     vz1.append(state1['vz'])
    
    # tg1.init_state({'x':x1[-1], 'y':y1[-1], 'z':z1[-1], 'vx':vx1[-1], 'vy':0.0, 'vz':0.0, 'ax': 0.0, 'ay': 0.0, 'az': 0.0})

    # for i in range(nCv):
    #     state1 = tg1.CV(dt)
    #     x1.append(state1['x'])
    #     vx1.append(state1['vx'])

    #     y1.append(state1['y'])
    #     vy1.append(state1['vy'])

    #     z1.append(state1['z'])
    #     vz1.append(state1['vz'])

    # tg1.init_state({'x':x1[-1], 'y':y1[-1], 'z':z1[-1], 'vx':vx1[-1],'vy':0.0,'vz':0.0,'w':5.5})

    # for i in range(nCt):
    #     state1 = tg1.CTxz(dt)
    #     x1.append(state1['x'])
    #     vx1.append(state1['vx'])

    #     y1.append(state1['y'])
    #     vy1.append(state1['vy'])

    #     z1.append(state1['z'])
    #     vz1.append(state1['vz'])

    # tg1.init_state({'x':x1[-1], 'y':y1[-1], 'z':z1[-1], 'vx':200.0, 'vy':0.0, 'vz':-40.0, 'ax': 0.0, 'ay': 0.0, 'az': 0.0})
    # for i in range(nCv2):
    #     state1 = tg1.CV(dt)
    #     x1.append(state1['x'])
    #     vx1.append(state1['vx'])

    #     y1.append(state1['y'])
    #     vy1.append(state1['vy'])

    #     z1.append(state1['z'])
    #     vz1.append(state1['vz'])


    # tg1.init_state({'x':x1[-1], 'y':y1[-1], 'z':z1[-1], 'vx':vx1[-1],'vy':0.0,'vz':0.0,'w':5.5})
    for i in range(nCt):
        state1 = tg1.CTxz(dt)
        x1.append(state1['x'])
        vx1.append(state1['vx'])

        y1.append(state1['y'])
        vy1.append(state1['vy'])

        z1.append(state1['z'])
        vz1.append(state1['vz'])


    tg1.init_state({'x':x1[-1], 'y':y1[-1], 'z':z1[-1], 'vx':vx1[-1],'vy':0.0,'vz':0.0,'ax': 0.0, 'ay': 0.0, 'az': 0.0})
    for i in range(nCv):
        state1 = tg1.CV(dt)
        x1.append(state1['x'])
        vx1.append(state1['vx'])

        y1.append(state1['y'])
        vy1.append(state1['vy'])

        z1.append(state1['z'])
        vz1.append(state1['vz'])


    tg1.init_state({'x':x1[-1], 'y':y1[-1], 'z':z1[-1], 'vx':vx1[-1],'vy':0.0,'vz':0.0,'ax': 0.0, 'ay': 0.0, 'az': -10.0})
    for i in range(nCaPitch):
        state1 = tg1.CA(dt)
        x1.append(state1['x'])
        vx1.append(state1['vx'])

        y1.append(state1['y'])
        vy1.append(state1['vy'])

        z1.append(state1['z'])
        vz1.append(state1['vz'])


    tg1.init_state({'x':x1[-1], 'y':y1[-1], 'z':z1[-1], 'vx':vx1[-1],'vy':0.0,'vz':0.0,'ax': 0.0, 'ay': 0.0, 'az': 0.0})
    for i in range(nCv):
        state1 = tg1.CV(dt)
        x1.append(state1['x'])
        vx1.append(state1['vx'])

        y1.append(state1['y'])
        vy1.append(state1['vy'])

        z1.append(state1['z'])
        vz1.append(state1['vz'])

    tg1.init_state({'x':x1[-1], 'y':y1[-1], 'z':z1[-1], 'vx':vx1[-1],'vy':0.0,'vz':0.0,'ax': 0.0, 'ay': 0.0, 'az': 0.0, 'w': 5.5})
    for i in range(nCt):
        state1 = tg1.CTxz(dt)
        x1.append(state1['x'])
        vx1.append(state1['vx'])

        y1.append(state1['y'])
        vy1.append(state1['vy'])

        z1.append(state1['z'])
        vz1.append(state1['vz'])

    X_true_data_not_pass = np.array([x1,vx1,y1,vy1,z1,vz1])             

    return (X_true_data_not_pass)


X_true_data_not_pass = make_true(tg1)


X_true_data_with_pass, pass_index = make_pass(X_true_data_not_pass, pd)
with_pass = remove_zero_columns(X_true_data_with_pass)
# print("with_pass",with_pass)

# Создаем трехмерный scatter plot для массива
scatter2 = go.Scatter3d(x=with_pass[0], y=with_pass[2], z=with_pass[4], mode='markers+lines', marker=dict(size=3, color='blue'), name='X_true_data')
fig2.add_trace(scatter2)

# Обновляем параметры макета
fig2.update_layout(scene=dict(aspectmode="cube"))
fig2.show()


# # ================= Блок 2 =================
# # создание истинной зашумленной траектории
process_var = 10
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

#=================================================================
#построение графика истинной зашумленной траектории с пропусками

X_true_plus_ProcNoise_not_pass = X_true_plus_ProcNoise_with_pass
X_true_plus_ProcNoise_not_pass[:,pass_index] = 0
X_true_plus_ProcNoise_not_pass = X_true_plus_ProcNoise_with_pass[:, ~pass_index] # вырезаю пропуски, что построить график
# print("X_true_plus_ProcNoise_not_pass" , X_true_plus_ProcNoise_not_pass)

plt.plot(X_true_plus_ProcNoise_not_pass[0],X_true_plus_ProcNoise_not_pass[4], label='X_true_plus_ProcNoise+pass', linestyle='-', marker='o')
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
    r_true_with_noise = np.sqrt(np.array(X_plusProcNoise[0])**2 + np.array(X_plusProcNoise[2])**2 + np.array(X_plusProcNoise[4])**2)
    az_true_with_noise = np.rad2deg(np.arctan2(np.array(X_plusProcNoise[2]),np.array(X_plusProcNoise[0]))) # Азимут
    um_true_with_noise = np.rad2deg(np.arctan2(np.array(X_plusProcNoise[4]),np.sqrt(np.array(X_plusProcNoise[0])**2+np.array(X_plusProcNoise[2])**2))) # Угол места

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

    detection = estimator.Detection()

    r_meas = Z[0,0]
    az_meas = Z[1,0] 
    um_meas = Z[2,0] 
    
    meas = np.array([[r_meas],[az_meas],[um_meas]])

    detection.point = meas
    detection.timePoint = dt
   
    track = estimator.BindTrackUkfImm_CTxz(detection) #инициал. трассы
    
    X_c = np.empty((6, 0))
    m_i = np.empty((0, 3))

    for i in range (1, Z.shape[1]):
        print(i)
        m_i = np.append(m_i, track.get_m_i(), axis=0)
        r_meas = Z[0,i]
        az_meas = Z[1,i]
        um_meas = Z[2,i]
        meas = ([[r_meas],[az_meas],[um_meas]])

        detection.point = meas
        detection.timePoint = (i * dt) + dt
        
        if np.all(Z[:,i] == 0):
            X = track.step(detection.timePoint)
            X_c = np.append(X_c,X,axis=1)
            continue
        # print('Z=',Zsph2cart(Z[:,i]))

        X = track.step(detection)

        # print('X=',X)
        
        
        X_c = np.append(X_c,X,axis=1)
        # print('Xc=',X_c)
    # print("X_Estimeted=",X_c)
    # print ('m_i=',m_i)
        
    return X_c, m_i 

X_c, m_i = estimate(Z)
# # print("X_Estimeted=",X_c)

# # def err1(X_c,X_true_plus_ProcNoise):

# #     er = X_c[:,:] - X_true_plus_ProcNoise [:,1:]
# #     return er

# # e = err1(X_c, X_true_plus_ProcNoise)
# # # print("e =", e[0]) #Распечатка ошибок по Х одной трассы


# # #==================Отрисовка==================
# # Z_cart = Zsph2cart(Z)

plt.figure()
plt.plot(X_c[0], X_c[4], label='Correct', marker='o')
plt.plot(X_true_plus_ProcNoise[0],X_true_plus_ProcNoise[4], label='truth', marker='x')
plt.plot(Zc[0], Zc[2], label='Meas',marker='o')
plt.legend()

plt.figure()
plt.plot((np.arange(len(m_i[:, 0]))+1)*dt, m_i[:,0], label='m_i_CV', marker='o')
plt.plot((np.arange(len(m_i[:, 0]))+1)*dt,m_i[:,1], label='m_i_CT', marker='x')
plt.plot((np.arange(len(m_i[:, 0]))+1)*dt, m_i[:,2], label='m_i_CA',marker='o')
plt.legend()


# # # # ================= Блок 5 ===================
# # # СБОР СТАТИСТИКИ

def calc_err(X):

    Xn = add_process_noise(X, Q)
    X_pass, pass_id = make_pass(Xn,pd)
    Zn = do_measurement(X_pass, R, pass_id)
    X_c, m_i = estimate(Zn)

    err = X_c[:,:] - Xn [:6,1:] # ошибка вычисляется со второго столбца.

    # print("ошибка в статистике calc_err",err[0])

    return err

from tqdm import tqdm

def calc_std_err(X):
    num_iterations = 1
    var_err = np.zeros((X.shape[0], X.shape[1]-1)) ## изменения так как модель CV. а настроено было под CT (моделируется CT. обрезал последний w)

    for i in tqdm(range(num_iterations)):
        err = calc_err(X)
        var_err += err ** 2

    var_err /= num_iterations
    return np.sqrt(var_err)

tg2 = Target()
tg2.init_state(init_state)
X_true_data_not_pass_2 = make_true(tg2)
std_err_2 = calc_std_err(X_true_data_not_pass_2)


plt.figure(num="track")
plt.subplot(6, 1, 1)
plt.plot((np.arange(len(std_err_2[0, :]))+1)*dt, std_err_2[0, :])
plt.xlabel('Time,s')
plt.ylabel('std_x, met')
plt.grid(True)
plt.subplot(6, 1, 2)
plt.plot((np.arange(len(std_err_2[1, :]))+1)*dt, std_err_2[1, :])
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_vx, m/s')
plt.subplot(6, 1, 3)
plt.plot((np.arange(len(std_err_2[2, :]))+1)*dt, std_err_2[2, :])
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_y, met')
plt.subplot(6, 1, 4)
plt.plot((np.arange(len(std_err_2[3, :]))+1)*dt, std_err_2[3, :])
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_vy, m/s')
plt.subplot(6, 1, 5)
plt.plot((np.arange(len(std_err_2[4, :]))+1)*dt, std_err_2[4, :])
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_z, met')
plt.subplot(6, 1, 6)
plt.plot((np.arange(len(std_err_2[5, :]))+1)*dt, std_err_2[5, :])
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_vz, m/s')

plt.subplots_adjust(wspace=12.0, hspace=1.0)

plt.show()