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

    def __init__(self, rho_L, rho_R, v_L, v_R, T=1, mu=1):
        '''
        Initialize with Riemann initial conditions

        input:
        - rho_L (float): density for x<0
        - rho_R (float): density for x>0
        - v_L (float): velocity for x<0
        - v_R (float): velocity for x>0
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
        self.v_L = v_L
        self.v_R = v_R

        self.c_s = np.sqrt(self.T*self.R_gas)
    
    def hugoniot_locus(self, U, rho_range):
        '''
        function that calculates the hugoniot curves in state space for a given state U = (rho, m)

        input:
        - U (2-tuple): contains rho_0 and v_0 
        - rho_range (array): array containing a range of rho-values for which we calculate pos/neg Hugoniot locus

        output:
        - array containing the values for 2 diff Hugoniot curves
        '''
        rho_hat, m_hat = U

        # calculate hugoniot loci
        m_pos = rho_range/rho_hat*m_hat + np.sqrt(rho_range/rho_hat)*self.c_s*(rho_range-rho_hat)
        m_neg = rho_range/rho_hat*m_hat - np.sqrt(rho_range/rho_hat)*self.c_s*(rho_range-rho_hat)

        return np.array([m_pos, m_neg])

    def Hugoniot_derivative(self, U):
        '''
        Calculate slope of derivatives to the Hugoniot curves

        input:
        - U (array): contains the state of which we want to calculate local derivatives

        output:
        - array: containing two derivatives
        '''
        if U[0] == 0.0:
            return ValueError("Derivative cannot be calcuated when density is zero.")
        else:
            v_hat = U[1]/U[0]
            der_pos = v_hat + self.c_s
            der_neg = v_hat - self.c_s
            return np.array([der_pos, der_neg])
        
    def shock_speed(self, U, rho_range):
        '''
        Function that returns the shock pseed of Rankine_hugoniot system

        input: 
        - U (array): contains an intial state, array has 2 elements
        - rho_range (array): for which rho_values do we calculate shock speed

        output:
        - s (array): contains shock speed for each given rho_value
            2-dimensional: there's positive and negative solution
        '''
        rho_hat, m_hat = U[0], U[1]

        s_positive = m_hat/rho_hat + np.sqrt(rho_range/rho_hat)*self.c_s
        s_negative = m_hat/rho_hat - np.sqrt(rho_range/rho_hat)*self.c_s

        return np.array([s_positive, s_negative])


            
    def integral_curve(self, U, rho_range, wave_type):
        '''
        Compute integral curves (rarefactions) for given state U

        input:
        - U (2-tuple): contains initial state
        - rho_range (array): values for which to calculate rarefaction
        - wave_type (str): specify wave type
            options include '1-rarefaction', '2-rarefaction'

        output:
        - array of (rho, v) values
        '''
        rho_hat, m_hat = U
        v_hat = m_hat/rho_hat

        # check input type
        if (wave_type != '1-rarefaction') and (wave_type != '2-rarefaction'):
            raise ValueError("Invalid wave type. Choose '1-rarefaction' or '2-rarefaction'.")

        # proceed to calculation 
        if wave_type == '1-rarefaction':
            m = v_hat*rho_range - self.c_s * rho_range * np.log(rho_range / rho_hat)
        elif wave_type == '2-rarefaction':
            m = v_hat*rho_range + self.c_s * rho_range * np.log(rho_range / rho_hat)
        return m

    def solve_riemann_problem(self, 
                            method='tvdlf', 
                            limiter='minmod', 
                            timestepper='onestep', 
                            nx=100,
                            cfl=0.5, 
                            t_end=2
                            ):
        '''
        Function that solves the 1D Riemann problem numerically. 

        input: 
        - method (str): which flux discretization method to use
            Options: 'tvdlf', 'upwind', 'maccormack'
        - limiter (str): which limiter for the flux reconstruction to use. 
            Options: 'minmod', 'woodward'
        - timestepper (str): which numerical scheme to use for the numerical time integraton of fluxes
            Options: 'onestep', 'twostep'
        - nx (int): number of discrete spatial cells
        - cfl (foat): determines maximal time step to enforce TVD
        - t_end (float): time for which we simulate behaviour of the system

        output:
        - rho (array): final values for density
        - momentum (array): final values for momentum
        '''
        methods = ['tvdlf', 'upwind', 'maccormack']
        if method not in methods:
            raise ValueError(("Invalid flux scheme. Choose 'tvdlf', 'upwind' or 'maccormack'."))
        
        timesteppers = ['onestep', 'twostep']
        if timestepper not in timesteppers:
            raise ValueError(("Invalid time stepper. Choose 'onestep' or 'twostep'."))
        
        limiters = ['minmod', 'woodward']
        if limiter not in limiters:
            raise ValueError("Invalid limiter. Choose from 'minmod' or 'woodward'.")

        # initial discretization setup
        nx = 100
        x = np.linspace(-1, 1, nx)
        dx = x[1] - x[0]
        v_max = max(abs(self.v_L), abs(self.v_R))+self.c_s
        dt = cfl * dx / v_max  # CFL condition

        # initialize U: density and momentum vector
        U = np.zeros((nx, 2))
        mid = nx//2

        # initialize left and right state with Riemann initial conditions
        U[:mid, 0] = self.rho_L
        U[:mid, 1] = self.v_L * self.rho_L
        U[mid:, 0] = self.rho_R
        U[mid:, 1] = self.v_R * self.rho_R

        # do time integration (spatial flux calculation is done under the hood)
        t = 0
        while t < t_end:
            if timestepper == 'onestep':
                U = self.timestep(U, dx, dt, method, limiter)
            elif timestepper == 'twostep':
                # employ 2 step predictor corrector method
                U_predictor =self.timestep(U, dx, dt, method, limiter)
                U = 1/2 * (U + self.timestep(U_predictor, dx, dt, method, limiter))
            t += dt
        
        # return final value
        rho = U[:, 0]
        momentum = U[:, 1]
        return x, rho, momentum
    
    def flux(self, U):
        '''
        Function that calculates flux based on a given state U[i]

        input:
        - U (2-dim array): contains rho, momentum for given spatial grid point

        output:
        - outward flux based on analytical calculation
        '''
        rho = U[0]
        v = U[1]/U[0]
        return np.array([rho * v, rho * v**2 + self.c_s**2 * rho])
    
    def limit_slope(self, slope_L, slope_R, limiter):
        '''
        Function that limits the slope reconstruction.

        input: 
        - slope_L (float): slope on the left boundary
        - slope_R (float): slope on the irght boundary
        - limiter (str): which limiting procedure to use
            Options: 'minmod', 'woodward'

        output:
        - limited slope of the conserved variables.    
        '''
        if limiter == 'minmod':
            slope_lim = np.sign(slope_L) * np.maximum(0, np.minimum(abs(slope_L), abs(slope_R)))
        elif limiter == 'woodward':
            slope_lim =np.maximum(0, np.minimum((slope_L + slope_R) / 2, 2 * slope_L, 2 * slope_R))
        else:
            raise NotImplementedError("Method not implemented. Choose from 'minmod', 'woodward'. ")
        return slope_lim

    def timestep(self, U, dx, dt, method, limiter):
        '''
        Function that integrates over the time using different flux schemes and limiter methods

        input:
        - U (array): density and momentum over discretized grid
        - dx (float): spatial discretization interval
        - dt (float): time discretization interval
        - method (str): flux scheme method
            Options: 'tvdlf', 'maccormack', 'upwind'
        - limiter (str): limiter for flux reconstruction using interpolation method
            Options: 'minmod', 'woodward

        output:
        - U (array) after one timestep dt for all spatial cells
        '''
        nx = len(U)
        F_half = np.zeros((nx-1, 2))

        # compute fluxes at cell interfaces:
        for i in range(nx-1):
            # Compute slopes from limiters
            slope_L = U[i] - U[i-1]
            slope_R = U[i+1] - U[i]
            limited_slope = self.limit_slope(slope_L, slope_R, limiter)

            U_L = U[i]+limited_slope/2
            U_R = U[i+1] - limited_slope/2
            
            # Do the actual flux computation
            if method == 'tvdlf':
                # compute max speed
                max_speed = max(abs(U_L[1]/U_L[0])+self.c_s,
                                 abs(U_R[1]/U_R[0])+self.c_s)
                F_half[i] = 0.5*(self.flux(U_L)+ self.flux(U_R) - max_speed*(U_R - U_L))

            elif method == 'maccormack':
                # use predictor - corrector type method
                # predictor step: forward difference
                predictor = U_L - dt/dx * (self.flux(U_R)-self.flux(U_L))
                # corrector step: backward difference
                F_half[i] = 1/2 * (self.flux(U_L)+self.flux(predictor))
        
            elif method == 'upwind':
                # use no information from neighbouring cells
                # directly update using the flux from only this cell
                F_half[i] = self.flux(U_L)

        # Update the conserved variables
        U_new = U.copy()
        for i in range(1, nx-1):
            # Update all interior grid points using flux difference (conservation law)
            U_new[i] -= dt/dx *(F_half[i]-F_half[i-1])
        return U_new
    

    def get_intermediat_state(self, U_L, U_R):
        '''
        Function that calculates the intermediate state, connected through allowed 
        shocks to the left and right states

        input:
        - U_L (array): the left state
        - U_R (array): the right state

        output: 
        - U_M (array): the intermediate state
        '''
        rho_L, m_L = U_L[0], U_L[1]
        rho_R, m_R = U_R[0], U_R[1]

        # intermediate quantities
        rho_reduced = np.sqrt(rho_L*rho_R)/(np.sqrt(rho_L)+ np.sqrt(rho_R))
        v_diff = (m_L/rho_L - m_R/rho_R)/self.c_s


        # get inetrmediate solutions
        rho_m =  ( 0.5*v_diff*rho_reduced + 1/2*np.sqrt(v_diff**2 * rho_reduced**2 + 4*np.sqrt(rho_L*rho_R)))** 2
        rho_inter = (np.sqrt(rho_m)+ np.sqrt(rho_L*rho_R/rho_m))
        v_m = (m_L/rho_L- m_R/rho_R)/2 + self.c_s* (np.sqrt(rho_L) - np.sqrt(rho_R))/np.sqrt(rho_R*rho_L) * rho_inter
        return np.array([rho_m, v_m*rho_m])
    

    def solve_riemann(self):
        '''
        Function that solves the riemann problem for 1D oisothermal hydro based on 
        Hugoniot loci and integral curves. 

        input:
        - none: we use left and righ state given with class initialization. 

        output: 
        dictionary describing the final solution
        '''
        sol = {}
        rho_L = self.rho_L
        rho_R = self.rho_R

        v_L = self.v_L
        v_R = self.v_R

        U_L = np.array([rho_L, rho_L*v_L])
        U_R = np.array([rho_R, rho_R*v_R])

        lam_L = self.get_lambda([rho_L, rho_L*v_L])
        lam1_L = lam_L[0]
        lam1_L = lam_L[1]

        lam_R = self.get_lambda([rho_R, rho_R*v_R])
        lam1_R = lam_R[0]
        lam1_R = lam_R[1]

        # check trivial case
        if (rho_L == rho_R) and (v_L == v_R):
            wave_1 = 'absent'
            wave_2 = 'absent'

        # check degenerate cases
        cond1 = self.integral_curve(U_L, rho_R, '1-rarefaction')
        cond3 = self.integral_curve(U_R, rho_L, '2-rarefaction')
  
        # R on Hugoniot locus of L
        if (self.hugoniot_locus(U_L, rho_R)[0] == rho_R*v_R) or (self.hugoniot_locus(U_L, rho_R)[0] == rho_R*v_R):
            wave_1 = 'shock'
            wave_2 = 'absent'

        # R on integral curve of L
        
        elif (cond1 == rho_R*v_R):
            wave_1 = 'rarefaction'
            wave_2 = 'absent'

        # L on Hugoniot locus of R
        elif (self.hugoniot_locus(U_R, rho_L)[0] == rho_L*v_L) or (self.hugoniot_locus(U_R, rho_L)[0] == rho_L*v_L):
            wave_1 = 'absent'
            wave_2 = 'shock'

        # L on integral curve of R
        elif (cond3 == rho_L*v_L):
            wave_1 = 'absent'
            wave_2 = 'rarefaction'

        # check generic cases
        
        #check whether R lies above or below R1 / S1 curve
        else:
            if rho_R > rho_L:
                if self.hugoniot_locus(U_L, rho_R)[1] < v_R*rho_R:
                    # get value of R2(U_R, rho_L)
                    val = self.integral_curve(U_R, rho_L, '2-rarefaction')
                    if rho_L*v_L > val:
                        # intersects R1
                        wave_1 = 'rarefaction'
                        wave_2 = 'rarefaction'
                    else:
                        wave_1 = 'shock'
                        wave_2 = 'rarefaction'
                else:
                    wave_1 = 'shock'
                    wave_2 = 'shock'
            else: 
                if self.integral_curve(U_L, rho_R, '1-rarefaction') > v_R*rho_R:
                    # get value of R2(U_R, rho_L)
                    val = self.hugoniot_locus(U_R, rho_L)[0]
                    if rho_L*v_L > val:
                        # intersects R1
                        wave_1 = 'shock'
                        wave_2 = 'shock'
                    else:
                        wave_1 = 'rarefaction'
                        wave_2 = 'shock'
                else:
                    wave_1 = 'rarefaction'
                    wave_2 = 'rarefaction'


        # return solution  as dictionary
        sol['1-wave'] = wave_1
        sol['2-wave'] = wave_2
        return sol
    

    def get_lambda(self, U):
        '''
        function that calculates the eigenvalues of F_U for given state

        input:
        - U (array): given state vector

        output: 
        - eigenvalue array [lambda 1, lambda 2]
        '''
        lambda_1 = U[1]/U[0] + self.c_s
        lambda_2 = U[1]/U[0] - self.c_s

        return np.array([lambda_1, lambda_2])
    
    def rarefaction_exact_intermediate(self):
        '''
        Function that calculates the exact solution when intermediate state is
        connected to left and right state via rarefaction waves

        input:
        - none

        output: 
        - state vector U
        '''

        rho_L = self.rho_L
        v_L = self.v_L
        rho_R = self.rho_R
        v_R = self.v_R

        rho_m = np.sqrt(rho_L* rho_R)*np.exp((v_L-v_R)/2/self.c_s)
        v_m = (v_L + v_R)/2 + self.c_s/2* np.log(rho_L/rho_R)
        return np.array([rho_m, v_m * rho_m])
    
    