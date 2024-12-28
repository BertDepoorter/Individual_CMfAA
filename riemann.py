'''
This file contains some high-level class that takes the input and calls the appropriate functions 
to solve the 1D Riemann problem for isothermal Hydro
'''

# necessary imports
import numpy as np
import matplotlib.pyplot as plt

# define class that implements various quantities to calculate

class riemann:
    '''
    general class containing functionality to calculate different aspects 
    of the 1D Riemann problem for isothermal hydrodynamics
    '''

    def __init__(self, rho_L, rho_R, T=1, mu=1):
        '''
        Initialize with Riemann initial conditions

        input:
        - rho_L (float): density for x<0
        - rho_R (float): density for x>0
        (optional)
        - T (float): temperature at which we keep the heat bath
        - mu (float): molecular mass in atomic mass units
        '''
        self.T = T
        self.mu = mu

        kB = 1  # Boltzmann constant. We put this to 1 for easy 
        m_p = 1 # Proton mass, here taken to be atomic mass unit

        self.R_gas = kB/self.mu/m_p
        self.rho_L = rho_L
        self.rho_R = rho_R
    
    def hugoniiot_locus(self, U, rho_range):
        '''
        function that calculates the hugoniot curves in state space for a given state U = (rho, m)

        input:
        - U (2-tuple): contains rho_0 and v_0 
        - rho_range (array): array containing a range of rho-values for which we calculate pos/neg Hugoniot locus

        output:
        - array containing the values for 2 diff Hugoniot curves
        '''
        ci = np.sqrt(self.R_gas*self.T)


        rho_hat, m_hat = U
        ci = np.sqrt(self.R_gas*self.T)

        # calculate hugoniot loci
        m_pos = rho/rho_hat*m_hat + np.sqrt(rho/rho_hat)*ci*(rho-rho_hat)
        m_neg = rho/rho_hat*m_hat - np.sqrt(rho/rho_hat)*ci*(rho-rho_hat)

        return np.array([m_pos, m_neg])



