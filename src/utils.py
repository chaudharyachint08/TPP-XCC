import numpy as np
from scipy import integrate

import matplotlib.pyplot as plt

# class numeric_integration:
#     def __init__(self, name):
#         '''
#         Below mentioned following numerical integration methods are supported
#         Trapezoidal : 
#         Simpson's 1/3 : 
#         Simpson's 3/8 : 
#         Gauss-Kronrod : 
#         '''

# class numerical_derivative:
#     'Should be performed using Tensorflow Numerical integration'

class SyntheticGenerativeProcess:
    def __init__(self, name, epsilon=1e-7, **model_params):
        '''
        poisson
            Constant Intensity Function, intensity set to 'mu'
        hawkes
            Constant part followed by exciting part coming from part occurences of events
            intensity set to 'mu', exciting-weight set to 'w'
        self-correcting
            DESCRIPTION
        autoregressive-condition-duration
            DESCRIPTION
        '''
        self.name, self.epsilon = name, epsilon
        if self.name=='poisson':
            self.param_mu = model_params['mu']
        if self.name=='hawkes':
            self.param_mu, self.param_alpha, self.param_beta = model_params['mu'], model_params['alpha'], model_params['beta']
        if self.name=='self-correcting':
            self.param_mu, self.param_alpha = model_params['mu'], model_params['alpha']
        if self.name=='auto-regressive-conditional-duration':
            self.param_mu, self.param_alpha = model_params['mu'], model_params['alpha']
    
    def intensity(self, timestamp, seq_timesteps=[]):
        'Modeling the ultimate intensity function used in literature for statistical TPPs'
        if timestamp==np.inf or np.inf in seq_timesteps:
            raise Exception('Infinity as timestamp reached')
        if self.name=='poisson':
            return self.param_mu
        if self.name=='hawkes':
            return self.param_mu + self.param_alpha * np.exp(-self.param_beta*(timestamp-np.array(seq_timesteps))).sum()
        if self.name=='self-correcting':
            return np.exp(self.param_mu*timestamp - self.param_alpha*len(seq_timesteps))
        if self.name=='auto-regressive-conditional-duration':
            seq_timesteps = np.array(seq_timesteps)
            return 1/(self.param_mu + self.param_alpha*(seq_timesteps[1:]-seq_timesteps[:-1]).sum())
        return intensity_value
    
    def survived_intensity(self, timestamp, seq_timesteps=[], survival_bound_function=lambda x:np.exp(-x)):
        'This function does cross product of intensity_function with survival_probability* at given time'
        intensity_value = self.intensity(timestamp, seq_timesteps)
        if timestamp==np.inf or np.inf in seq_timesteps:
            raise Exception('Infinity as timestamp reached')
        if seq_timesteps :
            time_gap = timestamp - seq_timesteps[-1]            
            if not time_gap :
                survival_value = 0
            elif time_gap < 0 :
                raise Exception(f'Current timestep {timestamp} cannot be lesser than last timestep {seq_timesteps[-1]}')
            elif time_gap > 0 :
                survival_value, _ = integrate.quad(lambda x:self.intensity(x,seq_timesteps), seq_timesteps[-1], timestamp)
        else:
            survival_value = 1
        survived_intensity_value = intensity_value * survival_bound_function(survival_value)
        return survived_intensity_value

    def generate(self, seq_len):
        if seq_len <= 0:
            raise Exception(f'Sequence length to be generate should be greater than 0')
        seq_timesteps = []
        while len(seq_timesteps)<seq_len:
            prev_timestamp    = seq_timesteps[-1] if seq_timesteps else 0
            next_timestamp, _ = integrate.quad(lambda x:x*self.survived_intensity(x, seq_timesteps), prev_timestamp, np.inf)
            seq_timesteps.append( next_timestamp )
        min_time, max_time = min(seq_timesteps), max(seq_timesteps)
        return np.array(seq_timesteps)
