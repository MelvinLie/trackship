import numpy as np


class Track():

    def __init__(self, track_data, num_steps, particle, z_lim=31.0, z_scoring_plane=80.0, w_scoring_plane=2.0):
        '''Default constructor
        
        :param track_data:
            The track data.

        :param num_steps:
            The number of steps for the track data.

        :param particle:
            The particle.
            
        :param z_lim:
            The the limit in z. For particles reaching this limit we will continue their tracks
            to the scoring plane.
        
        :param z_scoring_plane:
            The z position of the scoring plane.

        :param w_scoring_plane:
            The width of the scoring plane.

        '''

        # the speed of light
        c = 299792458.0
        c_sq = c*c

        # get the state vectors
        self.states = track_data[:num_steps, :]

        # flag a hit
        self.hit = False

        # determine if we need to continue tracking to the target station
        if self.states[-1, 2] >= z_lim:

            # compute the momentum norm
            p_norm = np.linalg.norm(self.states[-1, 3:])

            # an auxiliary factor
            fac = c/np.sqrt(particle.mass*particle.mass*c_sq + p_norm*p_norm)

            # the derivative dz/dt
            dzdt = self.states[-1, 5]*fac

            # the time when the particle hits the scoring plane
            t_sp = (z_scoring_plane - self.states[-1, 2])/dzdt

            # the (x, z) position when the particle hits the scoring plane
            x_sp = self.states[-1, 3]*fac*t_sp + self.states[-1, 0]
            z_sp = self.states[-1, 5]*fac*t_sp + self.states[-1, 2]

            if (abs(x_sp) < w_scoring_plane):
                self.hit = True

            # append this datapoint
            self.states = np.append(self.states, np.array([[x_sp, 0., z_sp,
                                                            self.states[-1, 3],
                                                            self.states[-1, 4],
                                                            self.states[-1, 5]]]), axis=0)
            
    def plot(self, ax, color='mark_hit'):
        '''Plot this track into a matplotlib axes object.
        
        :param ax:
            The matplotlib axes.
            
        :param color:
            The color. Default 'mark_hit'. If 'mark_hit', then we color hits black and
            all other tracks gray.
        :return:
            None
        '''

        if isinstance(color, str):
            if (color == 'mark_hit'):
                if (self.hit == False):
                    ax.plot(self.states[:, 0], self.states[:, 2], color=[0.3, 0.3, 0.3, 0.2])
                else:
                    ax.plot(self.states[:, 0], self.states[:, 2], color='k')
            else:
                print('Color string {} unknown!')
        else:
            ax.plot(self.states[:, 0], self.states[:, 2], color=color)