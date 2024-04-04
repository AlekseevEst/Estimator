import numpy as np

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
    
    
    def CV(self,dt):
        # move target with CV model due dt
        keys = ['x','vx','y','vy','z','vz']
        X = [self.targetState[i] for i in keys]
        F = [[1.0,  dt,     0.0,    0.0,   0.0,  0.0],
             [0.0,  1.0,    0.0,    0.0,   0.0,  0.0],
             [0.0,  0.0,    1.0,    dt,    0.0,  0.0],
             [0.0,  0.0,    0.0,    1.0,   0.0,  0.0],
             [0.0,  0.0,    0.0,    0.0,   1.0,   dt],
             [0.0,  0.0,    0.0,    0.0,   0.0,  1.0]]
        
    
        Xe = np.matmul(F,X)

        for (i,k) in enumerate(keys):
            self.targetState[k] = Xe[i]
        return {k:self.targetState[k] for k in self.targetState.keys()} # return new dictionary

    def CA(self,dt):
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
        
        
        Xe = np.matmul(F,X)

        for (i,k) in enumerate(keys):
            self.targetState[k] = Xe[i]
        # return new dictionary
        return {k:self.targetState[k] for k in self.targetState.keys()}
    
    def CT(self,dt):
        # move target with CT model due dt
        keys = ['x','vx','y','vy','z','vz','w']
        X = [self.targetState[i] for i in keys]
        w = self.targetState['w'] # save w value
        F = [[1.0,  1/w*np.sin(w*dt),     0.0,      -1/w*(1-np.cos(w*dt)),  0.0,    0.0,     0.0],
             [0.0,  np.cos(w*dt),         0.0,      -np.sin(w*dt),          0.0,    0.0,     0.0],
             [0.0,  1/w*(1-np.cos(w*dt)), 1.0,      1/w*np.sin(w*dt),       0.0,    0.0,     0.0],
             [0.0,  np.sin(w*dt),         0.0,       np.cos(w*dt),          0.0,    0.0,     0.0],
             [0.0,  0.0,                  0.0,      0.0,                    1.0,     dt,     0.0],
             [0.0,  0.0,                  0.0,      0.0,                    0.0,    1.0,     0.0],
             [0.0,  0.0,                  0.0,      0.0,                    0.0,    0.0,     1.0]]
        
        
        Xe = np.matmul(F,X)

        for (i,k) in enumerate(keys):
            self.targetState[k] = Xe[i]
        # return new dictionary
        return {k:self.targetState[k] for k in self.targetState.keys()}