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
    
    def hugoniot_locus(self, U, rho_range):
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
        m_pos = rho_range/rho_hat*m_hat + np.sqrt(rho_range/rho_hat)*ci*(rho_range-rho_hat)
        m_neg = rho_range/rho_hat*m_hat - np.sqrt(rho_range/rho_hat)*ci*(rho_range-rho_hat)

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
        c_i = np.sqrt(self.R_gas*self.T)

        # check input type
        if (wave_type != '1-rarefaction') and (wave_type != '2-rarefaction'):
            raise ValueError("Invalid wave type. Choose '1-rarefaction' or '2-rarefaction'.")

        # proceed to calculation 
        if wave_type == '1-rarefaction':
            v = v_hat - c_i * np.log(rho_range / rho_hat)
        elif wave_type == '2-rarefaction':
            v = v_hat + c_i * np.log(rho_range / rho_hat)
        return v

    def solve_riemann_problem(self, method, limiter, timestepper):
        '''
        Function that solves the 1D Riemann problem numerically. 

        input: 
        - method (str): which flux discretization method to use
            Options: 'tvdlf', 'upwind', 'maccormack'
        - limiter (str): which limiter for the flux reconstruction to use. 
            Options: 'minmod', 'woodward'
        - timestepper (str): which numerical scheme to use for the numerical time integraton of fluxes
            Options: 'onestep', 'twostep'
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
        
        # here mke sure the right subfunctions are called

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
