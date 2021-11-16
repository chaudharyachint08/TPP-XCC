from scipy.integrate import integrate

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
    def __init__(self, name, epsilon=1e-7. **model_params):
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
            self.param_mu, self.param_w = model_params['mu'], model_params['w']
        if self.name=='self-correcting':
            pass
        if self.name=='poisson':
            pass
    
    def intensity(self, timestamp, seq_timesteps=[]):
        'Modeling the ultimate intensity function used in literature for statistical TPPs'
        if self.name=='poisson':
            return self.param_mu
        if self.name=='hawkes':
            pass
        if self.name=='self-correcting':
            pass
        if self.name=='poisson':
            pass
        return intensity_value
    
    def probability(self, timestamp, seq_timesteps=[]):
        'This function does cross product of intensity_function with survival_probability* at given time'
        intensity_value = self.intensity(timestamp, seq_timesteps)
        if seq_timesteps :
            time_gap = timestamp - seq_timesteps[-1]
            if not time_gap :
                survival_value = 0
            elif time_gap < 0 :
                raise Exception(f'Current timestep {timestamp} cannot be lesser than last timestep {seq_timesteps[-1]}')
            elif time_gap > 0 :
                survival_value = integrate.quad(self.intensity, seq_timesteps[-1], timestamp)
        else:
            survival_value = 1
        probability_value = intensity_value * survival_value
        return probability_value
        
        

    def generate(self, seq_len):
        if seq_len <= 0:
            raise Exception(f'Sequence length to be generate should be greater than 0')
        seq_timesteps = []
        while len(seq_timesteps)<seq_len:
            next_timestamp = integrate.quad(self.intensity, seq_timesteps[-1], timestamp)
            seq_timesteps.append( next_timestamp )
        min_time, max_time = min(seq_timesteps), max(seq_timesteps)
        return seq_timesteps
