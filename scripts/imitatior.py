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

dt = 6.0
pd = 0.9
R = np.diag([1.0, 1e-4, 1e-4]) # в рад задается
plt.rcParams['figure.figsize'] = [10, 6]
fig2 = make_subplots(rows=1, cols=1, specs=[[{'type': 'scatter3d'}]])

class Target():
    """
    target model description
    """
    def __init__(self):
        self.targetState = {                                   # state vector
            'x':0.0,'vx':0.0,'ax':0.0, # x states
            'y':0.0,'vy':0.0,'ay':0.0, # y states
            'z':0.0,'vz':0.0,'az':0.0, # z states
            'w':0.0} # angle velocity

        self.rand_gen = np.random.default_rng()           # set seed for random numbers
        
    def init_state(self,state):
        initKeys = state.keys()                               # init targetState components
        for k in initKeys:
            self.targetState[k] = state[k]                    # set values with same keys
    
    
    def CV(self):
        # move target with CV model due dt
        keys = ['x','vx','y','vy','z','vz']
        X = [self.targetState[i] for i in keys]
        F = [[1.0,  dt,     0.0,    0.0,   0.0,  0.0],
             [0.0,  1.0,    0.0,    0.0,   0.0,  0.0],
             [0.0,  0.0,    1.0,    dt,    0.0,  0.0],
             [0.0,  0.0,    0.0,    1.0,   0.0,  0.0],
             [0.0,  0.0,    0.0,    0.0,   1.0,   dt],
             [0.0,  0.0,    0.0,    0.0,   0.0,  1.0]]
        G = [[dt**2/2,  0.0,         0.0],
             [dt,       0.0,         0.0],
             [0.0,      dt**2/2,     0.0],
             [0.0,      dt,          0.0],
             [0.0,      0.0,     dt**2/2],
             [0.0,      0.0,     dt     ]]
        
    
        Xe = np.matmul(F,X)

        for (i,k) in enumerate(keys):
            self.targetState[k] = Xe[i]
        return {k:self.targetState[k] for k in self.targetState.keys()} # return new dictionary

    def CA(self):
        # move target with CA model due dt
        keys = ['x','vx','ax','y','vy','ay','z','vz','az']
        X = [self.targetState[i] for i in keys]
        F = [[1.0,  dt,     dt**2/2,  0.0,    0.0,    0.0,     0.0,    0.0,       0.0],
             [0.0,  1.0,    dt,       0.0,    0.0,    0.0,     0.0,    0.0,       0.0],
             [0.0,  0.0,    1.0,      0.0,    0.0,    0.0,     0.0,    0.0,       0.0],
             [0.0,  0.0,    0.0,      1.0,    dt,     dt**2/2, 0.0,    0.0,       0.0],
             [0.0,  0.0,    0.0,      0.0,    1.0,    dt,      0.0,    0.0,       0.0],
             [0.0,  0.0,    0.0,      0.0,    0.0,    1.0,     0.0,    0.0,       0.0],
             [0.0,  0.0,    0.0,      0.0,    0.0,    0.0,     1.0,    dt,    dt**2/2],
             [0.0,  0.0,    0.0,      0.0,    0.0,    0.0,     0.0,    1.0,        dt],
             [0.0,  0.0,    0.0,      0.0,    0.0,    0.0,     0.0,    0.0,       1.0]]
        
        G = [[dt**2/2,  0.0,         0.0],
             [dt,       0.0,         0.0],
             [1.0,      0.0,         0.0],
             [0.0,      dt**2/2,     0.0],
             [0.0,      dt,          0.0],
             [0.0,      1.0,         0.0],
             [0.0,      0.0,     dt**2/2],
             [0.0,      0.0,          dt],
             [0.0,      0.0,         1.0]]
        
        Xe = np.matmul(F,X)

        for (i,k) in enumerate(keys):
            self.targetState[k] = Xe[i]
        # return new dictionary
        return {k:self.targetState[k] for k in self.targetState.keys()}
    



# =============== Блок 1 ===================
    
# ИНИЦИАЛИЗАЦИЯ МОДЕЛИ ДВИЖЕНИЯ
tg1 = Target()
tg1.init_state({'x':500.0,'y':0.0,'z':0.0, 'vx':200.0,'vy':0.0,'vz':0.0}) # к ключам добавить ax и ay будет модель CA
#tg1.init_process_noise({'ax':0.0,'ay':0.0}) # при нулях ax и ay воздействия на цель не будет. будет прямая траектория.


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


def make_true (tg1): 
    n = 100 # количество измерении

    x1=[];y1=[];z1=[]; vx1=[];vy1=[];vz1=[] 
    for i in range(n):
        state1 = tg1.CV()
        x1.append(state1['x'])
        y1.append(state1['y'])
        z1.append(state1['z'])
        vx1.append(state1['vx'])
        vy1.append(state1['vy'])
        vz1.append(state1['vz'])

    X_true_data_not_pass = np.array([x1,vx1,y1,vy1,z1,vz1])        

    return (X_true_data_not_pass)

X_true_data_not_pass = make_true(tg1)

X_true_data_with_pass, pass_index = make_pass(X_true_data_not_pass, pd)


with_pass = remove_zero_columns(X_true_data_with_pass)
# Создаем трехмерный scatter plot для массива
scatter2 = go.Scatter3d(x=with_pass[0], y=with_pass[2], z=with_pass[4], mode='markers+lines', marker=dict(size=3, color='blue'), name='X_true_data')
fig2.add_trace(scatter2)

# Обновляем параметры макета
fig2.update_layout(scene=dict(aspectmode="cube"))
fig2.show()


# ================= Блок 2 =================
# создание истинной зашумленной траектории
process_var = 0.5
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
plt.plot(X_true_plus_ProcNoise_not_pass[0],X_true_plus_ProcNoise_not_pass[2], label='X_true_plus_ProcNoise+pass', linestyle='-', marker='o')
plt.legend()
# print("X_true+noise", X_true_plus_ProcNoise)


# ================= Блок 3 =====================

def Zsph2cart(Z):
    for i in range(Z.shape[1]):
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

    Zm = np.zeros((R.shape[0], X_plusProcNoise.shape[1]))
    for i in range(Zm.shape[1]):
        Zm[0,i] = r_true_with_noise[i]
        Zm[1,i] = az_true_with_noise[i]
        Zm[2,i] = um_true_with_noise[i]

    Z_plus_noise = Zm + np.sqrt(R) @ np.random.normal(loc=0, scale=math.sqrt(1.0), size=(Zm.shape[0], Zm.shape[1]))
    Z_plus_noise[:, pass_index] = 0

    return Z_plus_noise

Z = do_measurement (X_true_plus_ProcNoise_with_pass, R, pass_index)


#==================Отрисовка==================
Z_not_pass = Z
Z_not_pass = remove_zero_columns(Z_not_pass)
Zc = Zsph2cart(Z_not_pass)
plt.plot(Zc[0],Zc[1], label='do_meas', linestyle='-', marker='x')
plt.legend()

# ================= Блок 4 ===================
k = 1.0
def estimate (Z):

    Z_cart = Zsph2cart(Z)
    X0 = np.vstack([Z_cart[0,0], 0., Z_cart[1,0], 0., Z_cart[2,0], 0.]) # инициализируем вектор состояния, равный первому измерению

    ukf = estimator.BindTrackUkf(X0,dt,Qp,R,k) #инициал. фильтра
    X_c = np.empty((len(X0), 0))
    for i in range (Z.shape[1]-1):
        if np.all(Z[:,i+1] == 0):
            X = ukf.step(dt)
            X_c = np.append(X_c,X,axis=1)
            continue
        X = ukf.step(Z[:,i+1])
        X_c = np.append(X_c,X,axis=1)
    return X_c 

X_c = estimate(Z)

def err1(X_c,X_true_plus_ProcNoise):

    er = X_c[:,:] - X_true_plus_ProcNoise [:,1:]
    return er

e = err1(X_c, X_true_plus_ProcNoise)
# print("e =", e[0]) #Распечатка ошибок по Х одной трассы


#==================Отрисовка==================
Z_cart = Zsph2cart(Z)
plt.figure()
plt.plot(X_c[0], X_c[2], label='Correct', marker='o')
plt.plot(X_true_plus_ProcNoise[0],X_true_plus_ProcNoise[2], label='truth', marker='x')
plt.plot(Zc[0], Zc[1], label='Meas',marker='o')
plt.legend()


# ================= Блок 5 ===================
# СБОР СТАТИСТИКИ

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
    # print ("X", X)
    num_iterations = 2000
    var_err = np.zeros((X.shape[0], X.shape[1]-1))

    for i in tqdm(range(num_iterations)):
        err = calc_err(X)
        var_err += err ** 2

    var_err /= num_iterations
    return np.sqrt(var_err)

std_err = calc_std_err(X_true_data_not_pass)

plt.figure()
plt.subplot(6, 1, 1)
plt.plot((np.arange(len(std_err[0, :]))+2)*dt, std_err[0, :].T)
plt.xlabel('Time,s')
plt.ylabel('std_x, met')
plt.grid(True)
plt.subplot(6, 1, 2)
plt.plot((np.arange(len(std_err[1, :]))+2)*dt, std_err[1, :].T)
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_vx, m/s')
plt.subplot(6, 1, 3)
plt.plot((np.arange(len(std_err[2, :]))+2)*dt, std_err[2, :].T)
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_y, met')
plt.subplot(6, 1, 4)
plt.plot((np.arange(len(std_err[3, :]))+2)*dt, std_err[3, :].T)
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_vy, m/s')
plt.subplot(6, 1, 5)
plt.plot((np.arange(len(std_err[4, :]))+2)*dt, std_err[4, :].T)
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_z, met')
plt.subplot(6, 1, 6)
plt.plot((np.arange(len(std_err[5, :]))+2)*dt, std_err[5, :].T)
plt.grid(True)
plt.xlabel('Time,s')
plt.ylabel('std_vz, m/s')
plt.subplots_adjust(wspace=12.0, hspace=1.0)
plt.show()

