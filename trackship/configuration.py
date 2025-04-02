import numpy as np

from .magnet import Magnet

class Configuration():

    def __init__(self):
        '''Default constructor
        
        '''

        self.magnets = []
        

    def add_magnet(self, cp_core, cp_yoke, b_goal, mag_type, mirror=True):
        '''Add a magnet to the configuration.
        
        :param cp_core:
            The corner points of the core.

        :param cp_yoke:
            The corner points of the yoke.

        :param b_goal:
            The desired field.

        :param mag_type:
            The magnet type.

        :param  mirror:
            Set this flag to mirror the magnet in the yz plane.

            
        :return:
            None.
        '''

        self.magnets.append(Magnet(cp_core, cp_yoke, b_goal, mag_type))

        if mirror:
            cp_core_c = cp_core.copy()
            cp_core_c[:, 0] *= -1.0
            cp_yoke_c = cp_yoke.copy()
            cp_yoke_c[:, 0] *= -1.0
            self.magnets.append(Magnet(cp_core_c, cp_yoke_c, b_goal, mag_type))

        return None

    def plot(self, ax):
        '''Plot the configuration.
        
        :param ax:
            The matplotlib axes object.
            
        :return:
            None.
        '''

        for m in self.magnets:

            m.plot(ax)

        return None

    def compute_total_surface(self):
        '''Compute the total surface.
            
        :return:
            The value of the total surface.
        '''
        S = 0.0

        for m in self.magnets:
            S += m.compute_core_surface()
            S += m.compute_yoke_surface()

        return S
    
    def get_data_arrays(self):
        '''Get the data as arrays for the handover to c.

        :return None:

        '''

        # tne number of magnets
        num_magnets = len(self.magnets)

        # the number of cornerpoints for the cores and yokes
        num_core = np.zeros((num_magnets,), dtype=np.int64)
        num_yoke = np.zeros((num_magnets,), dtype=np.int64)

        # the data arrays
        x_core = np.zeros((0, ))
        y_core = np.zeros((0, ))
        z_core = np.zeros((0, ))
        x_yoke = np.zeros((0, ))
        y_yoke = np.zeros((0, ))
        z_yoke = np.zeros((0, ))
        b_goal = np.zeros((num_magnets, ))
        types = np.zeros((num_magnets, ), dtype=np.int64)

        # fill it
        for i, m in enumerate(self.magnets):
            num_core[i] = m.cp_core.shape[0]
            num_yoke[i] = m.cp_yoke.shape[0]
            b_goal[i] = m.b_goal
            types[i] = m.type

            x_core = np.append(x_core, m.cp_core[:, 0])
            y_core = np.append(y_core, m.cp_core[:, 1])
            z_core = np.append(z_core, m.cp_core[:, 2])

            x_yoke = np.append(x_yoke, m.cp_yoke[:, 0])
            y_yoke = np.append(y_yoke, m.cp_yoke[:, 1])
            z_yoke = np.append(z_yoke, m.cp_yoke[:, 2])

        return num_magnets, num_core, x_core, y_core, z_core,\
                num_yoke, x_yoke, y_yoke, z_yoke, b_goal, types
    