# This file contains a class equipped with general methods for solving HD equations

class solver: 
    '''
    This class contains different numerical solvers
    '''

    def __init__(self, method='tvdlf', dimension=1):
        '''
        Initialize the class

        input:
        - method (str): numerical method that has to be used in order to solve the Riemann problem in isothermal hydro. 
        - dimension (int): dimension of the problem.
        '''
        self.method=method
        self.dim = dimension

    def __call__(self):
        pass