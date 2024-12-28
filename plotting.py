'''
This class contains some plotting functionality for the solutions obtained using different methods

Different plotting methods can visualize various aspects of the problem. 
'''

# import necessary classes

# import general libraries
import matplotlib.pyplot as pl
import numpy as np
import scipy as sp



# idea: use default scipy methods to solve this type of problem and compare the results.

class Plotting:
    '''
    Class containing different functions necessary for plotting. 
    Only once intialize the class: define some important constants necessary for nearly all problems. 

    We just list different functions that take the solutions to various problems and then immediately plot them
    '''
    def __init__(self, mu, T):
        '''
        Initialize class with some necassary general variables.

        Input:
        - mu (float): molecular mass of the gas molecules 
        - T (float): temperature of the gas (isothermal so the same for all gas molecules.)
        '''
        kB = 1  # Boltzmann constant, set to 1
        m_p = 1 # mass of 1 proton
        self.R_gas = kB/mu/m_p
        self.T = T


    def Hugoniot_locus(self, U, rho_lim=None, v_lim=None, name='Hugoniot_locus'):
        '''
        Exercise 4.3: function for plotting the solutions to the equations
        governing the Hugoniot locus.

        input:
        - U_hat (2-tuple): contains state U hat
        (optional)
        - rho_lim (list): set limit of the rho array
        - m_lim (list): constrain the y-axis limit

        output:
        - figure shown as output

        '''
        # define elementary vanriables

        rho_hat, m_hat = U
        ci = np.sqrt(self.R_gas*self.T)

        # create range of rho
        if rho_lim != None:
            if type(rho_lim).__name__ != 'list':
                raise ValueError
            rho = np.linspace(rho_lim[0], rho_lim[1], 500)
        rho = np.linspace(0, 2, 500)
        m_pos = rho/rho_hat*m_hat + np.sqrt(rho/rho_hat)*ci*(rho-rho_hat)
        m_neg = rho/rho_hat*m_hat - np.sqrt(rho/rho_hat)*ci*(rho-rho_hat)

        fig, ax = plt.subplots(1,1)
        ax.plot(rho, m_pos, color='blue', label='Positive solution')
        ax.plot(rho, m_neg, color='blue', label='Negative solution')
        ax.set_ylim(v_lim)

        name_fig = 'Plots/'+name+'.png'
        fig.savefig(name_fig, dpi=300)
        plt.show()
