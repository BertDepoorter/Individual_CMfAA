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
        rho_hat, v_hat = U
       

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
                            cfl=0.8, 
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
        x = np.linspace(0, 1, nx)
        dx = x[1] - x[0]
        dt = cfl * dx / self.c_s  # CFL condition

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
                U = timestep(U, dx, dt, method, limiter)
            elif timestepper == 'twostep':
                # employ 2 step predictor corrector method
                U_predictor = timestep(U, dx, dt, method, limiter)
                U = 1/2 * (U + timestep(U_predictor, dx, dt, method, limiter))
            t += dt
        
        # return final value
        rho = U[:, 0]
        momentum = U[:, 1]
        return rho, momentum
    
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
        return np.array([rho * v, rho * v**2 + c_s**2 * rho])
    
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
        c_s = np.sqrt(self.R_gas*self.T)

        # compute fluxes at cell interfaces:
        for i in range(nx-1):
            # Compute slopes from limiters
            slope_L = U[i] - U[i-1]
            slope_R = U[i+1] - U[i]
            limited_slope = limit_slope(slope_L, slope_R, limiter)

            # Do the actual flux computation
            if method == 'tvdlf':
                U_L = U[i]+limited_slope/2
                U_R = U[i+1] - limited_slope/2
                # compute max speed
                max_speed = max(abs(U_L[i,1]/U_L[i,0])+self.c_s,
                                 abs(U_R[i+1,1]/U_R[i+1,0])+self.c_s)
                F_half[i] = 0.5*(flux(U_L)+ flux(U_R) - max_speed*(U_R - U_L))

            elif method == 'maccormack':
                # use predictor - corrector type method
                # predictor step: forward difference
                predictor = U[i] - dt/dx * (flux(U[i+1])-flux(U[i]))
                # corrector step: backward difference
                F_half[i] = 1/2 * (flux(U[i])+flux(predictor))
        
            elif method == 'upwind':
                # use no information from neighbouring cells
                # directly update using the flux from only this cell
                F_half[i] = flux(U[i])

        # Update the conserved variables
        U_new = U.copy()
        for i in range(1, nx-1):
            # Update all interior grid points using flux difference (conservation law)
            U_new[i] -= dt/dx *(F_half[i]-F_half[i-1])
        return U_new
    
    
    
        
            





















    def _TVDLF_flux_step(self, U_range):
        pass

    def solve_riemann_problem_numerically(method, limiter, timestepper):
        """
        Solves the 1D Riemann problem numerically using specified methods.

        input:
        - method (str): flux discretization method ('tvdlf', 'upwind', 'maccormack')
        - limiter (str): flux reconstruction limiter ('minmod', 'woodward')
        - timestepper (str): time integration scheme ('onestep', 'twostep')

        output:
        - Numerical solution for the Riemann problem.
        """
        def flux(U):
            """Compute the flux for the isothermal hydrodynamics system."""
            rho, momentum = U
            v = momentum / rho
            return np.array([rho * v, rho * v**2 + c_s**2 * rho])

        def apply_limiter(slope_L, slope_R, limiter):
            """Apply the chosen limiter to compute the reconstructed slope."""
            if limiter == 'minmod':
                return np.sign(slope_L) * np.maximum(0, np.minimum(abs(slope_L), abs(slope_R)))
            elif limiter == 'woodward':
                return np.maximum(0, np.minimum((slope_L + slope_R) / 2, 2 * slope_L, 2 * slope_R))
            else:
                raise ValueError("Invalid limiter. Choose 'minmod' or 'woodward'.")

        def timestep(U, dx, dt, method):
            """Advance the solution one time step using the chosen method."""
            nx = len(U)
            F_half = np.zeros((nx - 1, 2))

            # Compute fluxes at cell interfaces
            for i in range(nx - 1):
                if method == 'tvdlf':
                    max_lambda = max(abs(U[i][1] / U[i][0]) + c_s, abs(U[i + 1][1] / U[i + 1][0]) + c_s)
                    F_half[i] = 0.5 * (flux(U[i]) + flux(U[i + 1]) - max_lambda * (U[i + 1] - U[i]))
                elif method == 'upwind':
                    F_half[i] = flux(U[i])
                elif method == 'maccormack':
                    # Predictor step (forward difference)
                    predictor = U[i] - dt / dx * (flux(U[i + 1]) - flux(U[i]))
                    # Corrector step (backward difference)
                    F_half[i] = 0.5 * (flux(U[i]) + flux(predictor))
                else:
                    raise ValueError("Invalid method. Choose 'tvdlf', 'upwind', or 'maccormack'.")

            # Update conserved variables
            U_new = U.copy()
            for i in range(1, nx - 1):
                U_new[i] -= dt / dx * (F_half[i] - F_half[i - 1])
            return U_new

        # Initial setup
        nx = 100
        x = np.linspace(0, 1, nx)
        dx = x[1] - x[0]
        dt = 0.5 * dx / c_s  # CFL condition

        # Initialize U (density and momentum)
        U = np.zeros((nx, 2))
        mid = nx // 2
        U[:mid, 0] = 1.0  # Left density
        U[:mid, 1] = 1.0 * 0.5  # Left momentum
        U[mid:, 0] = 0.5  # Right density
        U[mid:, 1] = 0.5 * -0.5  # Right momentum

        # Time integration
        t = 0
        t_end = 0.2
        while t < t_end:
            if timestepper == 'onestep':
                U = timestep(U, dx, dt, method)
            elif timestepper == 'twostep':
                # Perform a predictor-corrector step
                U_predictor = timestep(U, dx, dt, method)
                U = 0.5 * (U + timestep(U_predictor, dx, dt, method))
            else:
                raise ValueError("Invalid timestepper. Choose 'onestep' or 'twostep'.")
            t += dt

        # Extract final density and velocity
        rho = U[:, 0]
        v = U[:, 1] / rho
        return x, rho, v
