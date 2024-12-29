'''
This class contains some plotting functionality for the solutions obtained using different methods

Different plotting methods can visualize various aspects of the problem. 
'''

# import necessary classes
from riemann import riemann

# import general libraries
import matplotlib.pyplot as plt
import numpy as np
import scipy as sp



# idea: use default scipy methods to solve this type of problem and compare the results.

class plotting:
    '''
    Class containing different functions necessary for plotting. 
    Only once intialize the class: define some important constants necessary for nearly all problems. 

    We just list different functions that take the solutions to various problems and then immediately plot them
    '''
    def __init__(self, rho_L, rho_R, mu=1, T=1):
        '''
        Initialize class with some necassary general variables.

        Input:
        - rho_L (float): initial density x<0
        - rho_R (float): initial density x>0
        (optional)
        - mu (float): molecular mass of the gas molecules 
        - T (float): temperature of the gas (isothermal so the same for all gas molecules.)
        '''
        kB = 1  # Boltzmann constant, set to 1
        m_p = 1 # mass of 1 proton
        self.R_gas = kB/mu/m_p
        self.T = T

        # initialize Riemann class necessary to calculate stuff
        self.rho_R = rho_R
        self.rho_L = rho_L
        self.riemann = riemann(rho_L, rho_R)

    def plot_Hugoniot_locus(self, U, rho_lim=None, v_lim=None, name='Hugoniot_locus'):
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

        # create range of rho
        if rho_lim != None:
            if type(rho_lim).__name__ != 'list':
                raise ValueError
            rho = np.linspace(rho_lim[0], rho_lim[1], 500)
        rho = np.linspace(0, 2, 500)

        # get solution
        m = self.riemann.hugoniot_locus(U, rho)
        m_pos = m[0]
        m_neg = m[1]

        fig, ax = plt.subplots(1,1)
        ax.plot(rho, m_pos, color='blue', label='Positive solution')
        ax.plot(rho, m_neg, color='orange', label='Negative solution')
        ax.scatter(U[0], U[1], marker='*', s=70, label='Initial state')
        if v_lim != None:  
            ax.set_ylim(v_lim)

        ax.set_xlabel(r'Mass density $\rho$', fontsize=13)
        ax.set_ylabel(r"Momentum $m$", fontsize=13)
        ax.set_title('Hugoniot locus', fontsize=17)
        ax.legend()

        name_fig = 'Plots/'+name+'.png'
        fig.savefig(name_fig, dpi=300)
        plt.show()

    def plot_integral_curves(self, U, rho_lim=None, v_lim=None, name='Integral_curves'):
        '''
        Plot the integral curves for a given state U = (rho, m)

        input:
        - U (float): initial state for which we plot integral curves
        - rho_lim (list of 2 elements): set the limits for the rho axis
        - v_lim (list of 2 elements): set the limits for the y axis
        - name (str): set the name for the figure

        output:
        - figure
        '''
        # create range of rho
        if rho_lim != None:
            if type(rho_lim).__name__ != 'list':
                raise ValueError
            rho = np.linspace(rho_lim[0], rho_lim[1], 500)
        rho = np.linspace(0, 2, 500)

        v_1 = self.riemann.integral_curve(U, rho, wave_type='1-rarefaction')
        v_2 = self.riemann.integral_curve(U, rho, wave_type='2-rarefaction')

        fig, ax = plt.subplots(1, 1, figsize=(8,6))
        ax.plot(rho, v_1, linestyle='--', label='1-rarefaction')
        ax.plot(rho, v_2, linestyle='--', label='2-rarefaction')

        ax.legend()

        ax.set_xlabel(r'Mass density $\rho$', fontsize=13)
        ax.set_ylabel(r'velocity $v$', fontsize=13)
        ax.set_title('Integral curves', fontsize=17)

        name_plot = 'Plots/'+name+'.png'
        fig.savefig(name_plot, dpi=300)
        plt.show()
