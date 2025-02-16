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
    def __init__(self, rho_L, rho_R, v_L, v_R, mu=1, T=1):
        '''
        Initialize class with some necassary general variables.

        Input:
        - rho_L (float): initial density x<0
        - rho_R (float): initial density x>0
        - v_L (float): initial velocity x<0
        - v_R (float): initial velocity x>0
        (optional)
        - mu (float): molecular mass of the gas molecules 
        - T (float): temperature of the gas (isothermal so the same for all gas molecules.)
        '''
        kB = 1  # Boltzmann constant, set to 1
        m_p = 1 # mass of 1 proton
        self.R_gas = kB/mu/m_p
        self.T = T
        self.mu = mu

        # initialize Riemann class necessary to calculate stuff
        self.rho_R = rho_R
        self.rho_L = rho_L
        self.v_L = v_L
        self.v_R = v_R
        self.riemann = riemann(self.rho_L,
                               self.rho_R, 
                               self.v_L, 
                               self.v_R, 
                               T=self.T, 
                               mu=self.mu)

    def plot_Hugoniot_locus(self, U, rho_lim=None, v_lim=None, name='Hugoniot_locus', title='Hugoniot locus'):
        '''
        Exercise 4.3: function for plotting the solutions to the equations
        governing the Hugoniot locus.

        input:
        - U_hat (2-tuple): contains state U hat = (rho_hat, m_hat)
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

    def plot_Hugoniot_locus_derivative(self, U, 
                                       rho_lim=None, 
                                       v_lim=None, 
                                       name='Hugoniot_locus_derivative', 
                                       title='Hugoniot locus with derivatives'):
        '''
        Exercise 4.3: function for plotting the solutions to the equations
        governing the Hugoniot locus.

        input:
        - U_hat (2-tuple): contains state U hat = (rho_hat, m_hat)
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

        # get derivative
        derivatives = self.riemann.Hugoniot_derivative(U)
        line_1 = (rho-U[0])*derivatives[0] + U[1]
        line_2 = (rho-U[0])*derivatives[1] + U[1]

        fig, ax = plt.subplots(1,1)
        ax.plot(rho, m_pos, color='blue', label='Positive solution')
        ax.plot(rho, m_neg, color='orange', label='Negative solution')
        ax.scatter(U[0], U[1], marker='*', s=70, label='Initial state')
        if v_lim != None:  
            ax.set_ylim(v_lim)

        # Plot derivatives
        ax.plot(rho, line_1, linestyle='dashed', color='blue', alpha=0.6, label='Derivative (pos)')
        ax.plot(rho, line_2, linestyle='dashed', color='orange', alpha=0.6, label='Derivative (neg)')

        # make nice plot
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
        ax.plot(rho, m_pos, label='Hugoniot locus (pos)')
        ax.plot(rho, m_neg, label='Hugoniot locus (neg)')

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
                                title='Final state',
                                plot_momentum=False
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
        - plot_momentum (bool): plot momentum or not
            default is False
        '''
        params = {
            'method': method,
            'timestepper': timestepper,
            'limiter': limiter,
            't_end': t_end,
            'nx': nx
        }

        x, rho, momentum = self.riemann.solve_riemann_problem(**params)
        if plot_momentum  == True:
            fig, ax = plt.subplots(2, 1, figsize=(8, 14))
        
            ax[0].plot(x, rho, label=r'density $\rho$')
            ax[0].set_xlabel('x', fontsize=13)
            ax[0].set_ylabel(r'density $\rho$', fontsize=13)
            ax[0].set_title('Final density profile', fontsize=16)
            
            ax[1].plot(x, momentum, color='orange', label=r'Momentum $m$')
            ax[1].set_xlabel('x', fontsize=13)
            ax[1].set_ylabel(r'Momentum $m$', fontsize=13)
            ax[1].set_title('Final momentum profile', fontsize=16)
            fig.suptitle(title, fontsize=18)
            ax[0].legend(), ax[1].legend()
            name_fig = 'Plots/'+name+'.png'
            fig.savefig(name_fig, dpi=300)
            plt.show()
        else:
            fig, ax = plt.subplots(1, 1, figsize=(8, 6))
        
            ax.plot(x, rho, label=r'density $\rho$')
            ax.set_xlabel('x', fontsize=13)
            ax.set_ylabel(r'density $\rho$', fontsize=13)
            ax.set_title('Final density profile', fontsize=16)
        
            fig.suptitle(title, fontsize=18)
            ax.legend()
            name_fig = 'Plots/'+name+'.png'
            fig.savefig(name_fig, dpi=300)
            plt.show()

    def plot_numerical_state_space(self, 
                                   method, 
                                   limiter,
                                   timestepper,
                                   t_end=2.0,
                                   nx=100,
                                   name='numerical_state_space',
                                   title='Numerical solution in state space'
                                ):
        '''
        Function that plots the numerical solution in state space, i.e . in the (rho, m) plane

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

        output: 
        - figure (also saved in folder plots)
        '''
        params = {
            'method': method,
            'timestepper': timestepper,
            'limiter': limiter,
            't_end': t_end,
            'nx': nx
        }

        x, rho, momentum = self.riemann.solve_riemann_problem(**params)

        fig, ax = plt.subplots(1, 1, figsize=(8,6))
        ax.plot(rho, momentum, color='blue', label='Numerical solution')

        ax.legend()
        ax.set_title(title, fontsize=16)
        ax.set_xlabel(r'Mass density $\rho$', fontsize=13)
        ax.set_ylabel(r'Momentum $m$', fontsize=13)
        name_file = 'Plots/'+name+'.png'
        fig.savefig(name_file, dpi=300)

        plt.show()

    def Plot_Rankine_hugoniot(self, U, 
                              rho_lim=None,
                              title='Rankine-Hugoniot solution', 
                              name='Plot_Rankine_hugoniot'):
        '''
        Function that creates a plot of the 1-parameter solution to the Rankine-Hugoniot system. 
        The first plot shows the positive and negative solution for the momentum, 
        the second one shows the solution for shock speed. 

        input:
        - U (array): contains rho_hat and m_hat, the initial state
        (optional)
        - rho_range (2-tuple/array): contains limits on rho-axis
        - title (str): title for the plot
        - name (str): name for the saved figure
        '''

        # create rho array
        if rho_lim != None:
            if type(rho_lim).__name__ != 'list':
                raise ValueError
            rho = np.linspace(rho_lim[0], rho_lim[1], 500)
        rho = np.linspace(0, 2, 500)

        # get solutions
        m = self.riemann.hugoniot_locus(U, rho)
        m_pos = m[0]
        m_neg = m[1]

        s = self.riemann.shock_speed(U, rho)
        s_pos = s[0]
        s_neg = s[1]

        fig, ax = plt.subplots(2,1, figsize=(8, 14))
        ax[0].plot(rho, m_pos, color='blue', label='Positive solution')
        ax[0].plot(rho, m_neg, color='orange', label='Negative solution')
        ax[0].scatter(U[0], U[1], marker='*', s=70, label='Initial state')
        ax[0].set_title('Hugoniot locus: momentum', fontsize=15)
        ax[0].set_xlabel(r'Mass density $\rho$', fontsize=12)
        ax[0].set_ylabel(r'Momentum $m$', fontsize=12)

        ax[1].plot(rho, s_pos, color='blue', label='Positive solution')
        ax[1].plot(rho, s_neg, color='orange', label='Negative solution')
        ax[1].set_title('Hugoniot locus: shock speed', fontsize=15)
        ax[1].set_xlabel(r'Mass density $\rho$', fontsize=12)
        ax[1].set_ylabel(r'shock speed $s$', fontsize=12)
        ax[0].legend(), ax[1].legend()

        name_fig = 'Plots/'+name+'.png'
        fig.savefig(name_fig, dpi=300)
        plt.show()

    def Plot_Rankine_hugoniot_derivative(self, U, 
                              rho_lim=None,
                              title='Rankine-Hugoniot solution', 
                              name='Plot_Rankine_hugoniot'):
        '''
        Function that creates a plot of the 1-parameter solution to the Rankine-Hugoniot system. 
        The first plot shows the positive and negative solution for the momentum, 
        the second one shows the solution for shock speed. 

        input:
        - U (array): contains rho_hat and m_hat, the initial state
        (optional)
        - rho_range (2-tuple/array): contains limits on rho-axis
        - title (str): title for the plot
        - name (str): name for the saved figure
        '''

        # create rho array
        if rho_lim != None:
            if type(rho_lim).__name__ != 'list':
                raise ValueError
            rho = np.linspace(rho_lim[0], rho_lim[1], 500)
        rho = np.linspace(0, 2, 500)

        # get solutions
        m = self.riemann.hugoniot_locus(U, rho)
        m_pos = m[0]
        m_neg = m[1]

        s = self.riemann.shock_speed(U, rho)
        s_pos = s[0]
        s_neg = s[1]

        # get derivative
        derivatives = self.riemann.Hugoniot_derivative(U)
        line_1 = (rho-U[0])*derivatives[0] + U[1]
        line_2 = (rho-U[0])*derivatives[1] + U[1]


        fig, ax = plt.subplots(2,1, figsize=(6, 10))
        ax[0].plot(rho, m_pos, color='blue', label='Positive solution')
        ax[0].plot(rho, m_neg, color='orange', label='Negative solution')
         # Plot derivatives
        ax[0].plot(rho, line_1, linestyle='dashed', color='blue', alpha=0.6, label='Derivative (pos)')
        ax[0].plot(rho, line_2, linestyle='dashed', color='orange', alpha=0.6, label='Derivative (neg)')
        ax[0].scatter(U[0], U[1], marker='*', s=70, label='Initial state')
        ax[0].set_title('Hugoniot locus: momentum', fontsize=15)
        ax[0].set_xlabel(r'Mass density $\rho$', fontsize=14)
        ax[0].set_ylabel(r'Momentum $m$', fontsize=14)

        ax[1].plot(rho, s_pos, color='blue', label='Positive solution')
        ax[1].plot(rho, s_neg, color='orange', label='Negative solution')
        ax[1].set_title('Hugoniot locus: shock speed', fontsize=15)
        ax[1].set_xlabel(r'Mass density $\rho$', fontsize=14)
        ax[1].set_ylabel(r'shock speed $s$', fontsize=14)
        ax[0].legend(), ax[1].legend()
        fig.suptitle(title, fontsize=18)
        plt.tight_layout()
        name_fig = 'Plots/'+name+'.png'
        fig.savefig(name_fig, dpi=300)
        plt.show()

    def plot_2_hugoniot_loci(self, U_l, U_r, 
                             rho_range=None, 
                             intermediate=False,
                             name='plot_2_hugoniot_loci', 
                             title=" Hugoniot locus for left and right states"):
        '''
        Plot the hugoniot locus for 2 different states U_l and U_r on the same plot

        input:
        - U_l (array): array of left state
        - U_r (array): array of right state
        (optional)
        - rho_range(array): array with endpoints of self-specified rho range
            default is none
        - intermediate (bool): calculate intermediate state U_m if allowed
            default is False 
        - name (str): name for saved figure
        - title (str): title for the figure

        output
        - figure with two Hugoniot loci
        '''
        # create rho array
        if rho_range != None:
            rho = np.linspace(rho_range[0], rho_range[1], 500)
        else:
            rho = np.linspace(0, 2, 500)

        # get solutions
        m_l = self.riemann.hugoniot_locus(U_l, rho)
        m_r = self.riemann.hugoniot_locus(U_r, rho)

        fig, ax = plt.subplots(1,1, figsize=(8,6))
        ax.plot(rho, m_l[0], color='blue', label='Left state: positive')
        ax.plot(rho, m_r[0], linestyle='dashed', color='blue', label='Right state: positive')
        ax.plot(rho, m_l[1], color='orange', label='Left state: negative')
        ax.plot(rho, m_r[1], linestyle='dashed', color='orange', label='Right state: negative')

        ax.scatter(U_l[0], U_l[1], label='Initial left state', s=70, marker='*')
        ax.scatter(U_r[0], U_r[1], label='Initial right state', s=70, marker='*')

        # get intermediate state
        if intermediate == True:
            U_m = self.riemann.get_intermediat_state(U_l, U_r)
            ax.scatter(U_m[0], U_m[1], label='Intermediate state', s=80, marker='*', color='red')

        ax.annotate(r"U_L", (U_l[0], 1.2*U_l[1]))
        ax.annotate(r"U_R", (U_r[0], 1.1*U_r[1]+0.25))

        ax.legend()
        ax.set_title(title, fontsize=17)
        ax.set_xlabel(r'Mass density $\rho$', fontsize=13)
        ax.set_ylabel(r'Momentum $m$', fontsize=13)

        name_fig = 'Plots/'+name+'.png'
        fig.savefig(name_fig, dpi=300)
        plt.show()

    def full_solution(self,
                      name='Full_solution_state_space',
                      title='Riemann Problem Solution in State Space',
                      figure=True):
        '''
        Function that plots the final state space solution to the Riemann problem. 
        Integral curves and Hugoniot loci are plotted, with the actual solution emphasized in bold.

        input: 
        - name (str): name for saved figure
        - title (str): title for the plot. 
        (optional)
        - figure (bool): indicate whether or not to create & plot a figure
        output:
        - Full state space solution figure
        '''
        rho_L = self.rho_L
        rho_R = self.rho_R

        m_L = self.v_L*self.rho_L
        m_R = self.v_R*self.rho_R

        U_L = np.array([rho_L, m_L])
        U_R = np.array([rho_R, m_R])

        if figure == True:
            # Create Plot
            fig, ax = plt.subplots(1,1, figsize=(8,6))
                                
            # Plot the initial states
            ax.scatter(U_L[0], U_L[1], label='Left initial state', marker='*', s=70)
            ax.scatter(U_R[0], U_R[1], label='Right initial state', s=70, marker='*')

            # create rho arrays (lower and higher density for both states)
            L_lower = np.linspace(0.04, rho_L, 500)
            L_higher = np.linspace(rho_L, 2, 500)

            R_lower = np.linspace(0.04, rho_R, 500)
            R_higher = np.linspace(rho_R, 2, 500)

            # print decision tree solution
            sol = self.riemann.solve_riemann()
            print(sol)

            if sol['1-wave'] == 'shock' and (sol['2-wave']== 'shock'):
                # get intermediate state for 2 shocks
                U_m = self.riemann.get_intermediat_state(U_L, U_R)
                ax.scatter(U_m[0], U_m[1], s=70, marker='*', color='red', label='Intermediate state')
            if sol['1-wave'] == 'rarefaction' and (sol['2-wave']== 'rarefaction'):
                # get intermediate state for 2 shocks
                U_m = self.riemann.rarefaction_exact_intermediate()
                ax.scatter(U_m[0], U_m[1], s=70, marker='*', color='red', label='Intermediate state')

            R1 = self.riemann.integral_curve(U_L, L_lower, '1-rarefaction')
            R2 = self.riemann.integral_curve(U_R, R_lower, '2-rarefaction')
            S1 = self.riemann.hugoniot_locus(U_L, L_higher)[1]
            S2 = self.riemann.hugoniot_locus(U_R, R_higher)[0]

            ax.plot(L_lower, R1, color='blue', label='R1', linestyle='dashed')
            ax.plot(L_higher, S1, color='blue', label='S1')
            ax.plot(R_lower, R2, color='orange', label='R2', linestyle='dashed')
            ax.plot(R_higher, S2, color='orange', label='S2')

            ax.set_xlabel(r'Mass density $\rho$', fontsize=13)
            ax.set_ylabel(r'Momentum $m$', fontsize=13)
            ax.set_title(title, fontsize=17)
            ax.legend()

            name_fig = 'Plots/'+name+'.png'
            fig.savefig(name_fig, dpi=300)
            plt.show()
        elif figure == False:
            # only print qualitative solution to Riemann problem
            sol = self.riemann.solve_riemann()
            print(sol)
        else:
            raise ValueError("keyword 'figure' should be bool: True or False. ")
