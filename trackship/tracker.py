import numpy as np

from . import track_c_mod as track_c
from . tracks import Track

class RK4_ParticleTracker():

    def __init__(self, config):
        '''Default constructor
        
        :param config:
            The configuration.

        :param particle_name:
            The name of the particle. So far only Muons (default) are implemented.

        :param y0:
            The initial conditions. Default np.zeros((0, 6)).
        '''

        # set the configuration
        self.config = config

        # the speed of light
        self.c = 299792458  # m/s

        # set the limits
        self.x_max = 1e12
        self.x_min = -1e12
        self.y_max = 1e12
        self.y_min = -1e12
        self.z_max = 1e12
        self.z_min = -1e12

    def set_configuration(self, config):
        '''Set (overwrite) the configuration.
        
        :param config:
            A configuration object.
            
        :return None.
        
        '''
        self.config = config


    def set_x_limits(self, x_min, x_max):
        '''Set the tracking limits in x.

        :param x_min:
            The minimum x dimension.

        :param x_max:
            The maximum x dimension.
        
        return None.
        '''
        self.x_max = x_max
        self.x_min = x_min

    def set_y_limits(self, y_min, y_max):
        '''Set the tracking limits in y.

        :param y_min:
            The minimum y dimension.

        :param y_max:
            The maximum y dimension.
        
        return None.
        '''
        self.y_max = y_max
        self.y_min = y_min

    def set_z_limits(self, z_min, z_max):
        '''Set the tracking limits in z.

        :param z_min:
            The minimum z dimension.

        :param z_max:
            The maximum z dimension.
        
        return None.
        '''
        self.z_max = z_max
        self.z_min = z_min

    def track(self, particles, max_steps, h=1e-10, subsampling_rate=10,
              z_pos_sp=80.0, z_lim=31.0, w_sp=2.0):
        '''Track the particles.
      
        :param max_steps:
            The maximum number of steps.

        :param h:
            The step size. Default 1e-10.
            
        :param subsampling_rate:
            The subsampling rate.

        :param z_pos_sp:
            The z position of the scoring plane.

        :param z_lim:
            The z limit to determine where a particle reached end of the muon shield.

        :param w_sp:
            The width of the scoring plane.

        :return:
            None.
        '''
        
        # the number of stored samples
        num_t_stored = np.int32(max_steps/subsampling_rate)

        track_data, num_steps = track_c.track_rk4(self.config, particles, max_steps, h, x_lim=(self.x_min, self.x_max),
                          y_lim=(self.y_min, self.y_max), z_lim=(self.z_min, self.z_max),
                          subsampling=subsampling_rate)
        

        tracks = [Track(track_data[i, :, :], num_steps[i], particles[i],
                    z_lim=z_lim,
                    z_scoring_plane=z_pos_sp,
                    w_scoring_plane=w_sp) for i in range(track_data.shape[0])]

        return tracks