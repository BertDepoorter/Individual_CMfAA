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

    def plot_Hugoniot_locus(self, U, rho_lim=None, v_lim=None, name='Hugoniot_locus', title='Hugoniot locus'):
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
        ax.set_title(title, fontsize=17)
        ax.legend()

        name_fig = 'Plots/'+name+'.png'
        fig.savefig(name_fig, dpi=300)
        plt.show()

    def plot_integral_curves(self, U, rho_lim=None, v_lim=None, name='Integral_curves', title='Integral curves'):
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

        m_1 = self.riemann.integral_curve(U, rho, wave_type='1-rarefaction')
        m_2 = self.riemann.integral_curve(U, rho, wave_type='2-rarefaction')

        fig, ax = plt.subplots(1, 1, figsize=(8,6))
        ax.plot(rho, m_1, linestyle='--', label='1-rarefaction')
        ax.plot(rho, m_2, linestyle='--', label='2-rarefaction')
        ax.scatter(U[0], U[1], s=70, marker='*', label='Initial state')

        ax.legend()

        ax.set_xlabel(r'Mass density $\rho$', fontsize=13)
        ax.set_ylabel(r'Momentum $m$', fontsize=13)
        ax.set_title(title, fontsize=17)

        name_plot = 'Plots/'+name+'.png'
        fig.savefig(name_plot, dpi=300)
        plt.show()



    def plot_integral_hugoniot(self, 
                               U, 
                               title='Hugoniot locus vs. integral curve', 
                               rho_lim=None, 
                               v_lim=None, 
                               name='hugoniot_integral'):
        '''
        Function that compares the hugoniot locus and integral curve solution for a
        initial state U = (rho, m)

        input:
        - U (2-tuple): 
        - title (str): title for the plot
        - rho_lim (list of 2): set limits on rho array
        - v_lim (list of 2): constrain the y-axis
        - name (str): set name for the saved figure
        '''

        fig, ax = plt.subplots(1,1, figsize=(8,6))

        # create rho array
        if rho_lim != None:
            rho = np.linspace(rho_lim[0], rho_lim[1], 500)
        else: 
            rho = np.linspace(0.0, 2, 500)
        
        # get Hugoniot locus
        m = self.riemann.hugoniot_locus(U, rho)
        m_pos = m[0]
        m_neg = m[1]

        # get integral curves
        m_1 = self.riemann.integral_curve(U, rho, wave_type='1-rarefaction')
        m_2 = self.riemann.integral_curve(U, rho, wave_type='2-rarefaction')

        # plot the curves
        ax.plot(rho, m_pos, linestyle='--', label='Hugoniot locus (pos)')
        ax.plot(rho, m_neg, linestyle='--', label='Hugoniot locus (neg)')

        ax.plot(rho, m_1, linestyle=':', label='1-rarefaction')
        ax.plot(rho, m_2, linestyle=':', label='2-rarefaction')

        # Plot initial state
        ax.scatter(U[0], U[1], s=70, label='Initial state', marker='*')

        # Make nice plot
        ax.legend()
        ax.set_title(title, fontsize=16)
        ax.set_xlabel(r'Mass density $\rho$', fontsize=13)
        ax.set_ylabel(r'Momentum $m$', fontsize=13)

        name_fig = 'Plots/'+name+'.png'
        fig.savefig(name_fig, dpi=300)
        plt.show()

    def plot_numerical_solution(self, 
                                method, 
                                limiter, 
                                timestepper,
                                t_end=2.0,
                                nx=100,
                                name='Final_state_density_momentum', 
                                title='Final state'
                                ):
        '''
        This is a function that plots the final result for density and momentum after numericallly solving Riemann problem. 

        input:
        - method (str): flux scheme to use
        - limiter (str): limit the slope reconstruction
        - timestepper (str): how to do the time integration
        (optional)
        - t_end (float): for how long to run the simulation 
            default is 2.0
        - nx (int): number of spatial points
            default is 100
        - name (str): name of the saved figure
            Default is 'final_state_density_momentum'
        - title (str): title for the plot
            default is 'Final state'
        '''
        params = {
            'method': method,
            'timestepper': timestepper,
            'limiter': limiter,
            't_end': t_end,
            'nx': nx,
            'name': name,
            'title': title
        }

        x, rho, momentum = self.riemann.solve_riemann_problem(**params)

        fig, ax = plt.subplots(2, 1, figsize=(8, 14))
        
        ax[0].plot(x, rho, label=r'density $\rho$')
        ax[0].set_xlabel('x', fontsize=13)
        ax[0].set_ylabel(r'density $\rho$', fontsize=13)
        ax[0].set_title('Final density profile', fontsize=16)
        
        ax[1].plot(x, momentum, color='orange', label=r'Momentum $m$')
        ax[1].set_xlabel('x', fontsize=13)
        ax[1].set_ylabel(r'Momentum $m$', fontsize=13)
        ax[1].set_title('Final momentum profile', fontsize=16)

        ax[0].legend(), ax[1].legend()
        name_fig = 'Plots/'+name+'.png'
        fig.savefig(name_fig, dpi=300)
        plt.show()
