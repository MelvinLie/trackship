import numpy as np

class Particle():

    def __init__(self, mass, charge, initial_state):
        '''Default constructor
        
        :param mass:
            The particle mass.

        :param charge:
            The particle charge in electron charges.

        :param initial_state:
            The particles initial state.

        '''

        self.mass = mass
        self.charge = -1.602176634e-19*charge
        self.initial_state = initial_state