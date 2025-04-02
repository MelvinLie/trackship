# distutils: sources = [track_c/ctools.c, track_c/geometry.c, track_c/magnet2D.c, track_c/my_vectors.c, track_c/rktools.c]
# distutils: include_dirs = [track_c/]

import numpy as np
cimport numpy as np
from cython cimport view
from libc.stdlib cimport malloc, free
import copy

cimport ctools

def track_rk4(config, particles, n_max, h, x_lim=(-1e12, 1e12), y_lim=(-1e12, 1e12), z_lim=(-1e12, 1e12), subsampling=10):
    '''Track a collection of particles.

    :param config:
        The configuration.

    :param particles:
        The particles to track.

    :param n_max:
        The maximum number of steps.

    :param h:
        The step size.

    :param x_lim:
        The x limits.

    :param y_lim:
        The y limits.

    :param z_lim:
        The z limits.

    :param subsampling:
        The subsampling rate. Default 100.

    :return:
        The particle tracks.
    '''

    # ===========================================
    # Python Part
    # ===========================================

    
    # the total number of particles
    num_particles = len(particles)

    # the config information
    num_magnets, num_core, x_core, y_core, z_core,\
                num_yoke, x_yoke, y_yoke, z_yoke, b_goal, mag_types = config.get_data_arrays()

    # the number of particles
    num_particles = len(particles)

    # get the particle information
    p_m = np.zeros((num_particles, ))
    p_q = np.zeros((num_particles, ))
    p_y0 = np.zeros((6*num_particles, ))

    for i, p in enumerate(particles):
        p_m[i] = p.mass
        p_q[i] = p.charge

        p_y0[6*i:6*(i+1)] = p.initial_state

    # the number of stored samples
    num_t_stored = np.int32(n_max/subsampling)

    # the return data
    # tracks = np.zeros((num_particles*num_t_stored*6, ))

    # ===========================================
    # Convert python --> C
    # ===========================================
    cdef np.ndarray[int, ndim=1, mode = 'c'] num_core_buff = np.ascontiguousarray(num_core.flatten(), dtype = np.int32)
    cdef int* num_core_c = <int*> num_core_buff.data

    cdef np.ndarray[int, ndim=1, mode = 'c'] num_yoke_buff = np.ascontiguousarray(num_yoke.flatten(), dtype = np.int32)
    cdef int* num_yoke_c = <int*> num_yoke_buff.data

    cdef np.ndarray[double, ndim=1, mode = 'c'] x_core_buff = np.ascontiguousarray(x_core.flatten(), dtype = np.double)
    cdef double* x_core_c = <double*> x_core_buff.data
    
    cdef np.ndarray[double, ndim=1, mode = 'c'] y_core_buff = np.ascontiguousarray(y_core.flatten(), dtype = np.double)
    cdef double* y_core_c = <double*> y_core_buff.data

    cdef np.ndarray[double, ndim=1, mode = 'c'] z_core_buff = np.ascontiguousarray(z_core.flatten(), dtype = np.double)
    cdef double* z_core_c = <double*> z_core_buff.data

    cdef np.ndarray[double, ndim=1, mode = 'c'] x_yoke_buff = np.ascontiguousarray(x_yoke.flatten(), dtype = np.double)
    cdef double* x_yoke_c = <double*> x_yoke_buff.data
    
    cdef np.ndarray[double, ndim=1, mode = 'c'] y_yoke_buff = np.ascontiguousarray(y_yoke.flatten(), dtype = np.double)
    cdef double* y_yoke_c = <double*> y_yoke_buff.data

    cdef np.ndarray[double, ndim=1, mode = 'c'] z_yoke_buff = np.ascontiguousarray(z_yoke.flatten(), dtype = np.double)
    cdef double* z_yoke_c = <double*> z_yoke_buff.data

    cdef np.ndarray[double, ndim=1, mode = 'c'] b_goal_buff = np.ascontiguousarray(b_goal.flatten(), dtype = np.double)
    cdef double* b_goal_c = <double*> b_goal_buff.data

    cdef np.ndarray[int, ndim=1, mode = 'c'] types_buff = np.ascontiguousarray(mag_types.flatten(), dtype = np.int32)
    cdef int* types_c = <int*> types_buff.data

    cdef ctools.configuration config_c = ctools.make_configuration(num_magnets, num_core_c, num_yoke_c,
                                        x_core_c, y_core_c, z_core_c,
                                        x_yoke_c, y_yoke_c, z_yoke_c,
                                        b_goal_c, types_c)

    cdef np.ndarray[double, ndim=1, mode = 'c'] p_m_buff = np.ascontiguousarray(p_m.flatten(), dtype = np.double)
    cdef double* p_m_c = <double*> p_m_buff.data

    cdef np.ndarray[double, ndim=1, mode = 'c'] p_q_buff = np.ascontiguousarray(p_q.flatten(), dtype = np.double)
    cdef double* p_q_c = <double*> p_q_buff.data

    cdef np.ndarray[double, ndim=1, mode = 'c'] p_y0_buff = np.ascontiguousarray(p_y0.flatten(), dtype = np.double)
    cdef double* p_y0_c = <double*> p_y0_buff.data

    cdef ctools.particle *prtcls_c = ctools.make_particles(num_particles, p_m_c, p_q_c, p_y0_c)

    cdef ctools.track_limits lim_c =  ctools.make_track_limits(x_lim[0], x_lim[1],
                                                                y_lim[0], y_lim[1],
                                                                z_lim[0], z_lim[1])

    cdef int* num_steps_c = <int *> malloc(num_particles * sizeof(int))

    # the memory for the track data
    cdef double* tracks_c = <double *> malloc(6*num_particles*num_t_stored*sizeof(double))
    # cdef np.ndarray[double, ndim=1, mode = 'c'] tracks_buff = np.ascontiguousarray(track_data.flatten(), dtype = np.double)
    # cdef double* tracks_c = <double*> tracks_buff.data

    # ===========================================
    # Run C code
    # ===========================================
    ctools.track_particles(tracks_c, num_particles, h, n_max*h, subsampling, &lim_c, &config_c, prtcls_c, num_steps_c)

    # ===========================================
    # Convert C -> python
    # ===========================================
    cdef view.array tracks_array = <double[:num_t_stored*num_particles*6]> tracks_c
    track_data = np.asarray(tracks_array)
    track_data.shape = (num_particles, num_t_stored, 6)

    track_py = copy.deepcopy(track_data)

    cdef view.array num_steps_array = <int[:num_particles]> num_steps_c    
    num_steps = np.asarray(num_steps_array)

    # ===========================================
    # Clean up
    # ===========================================
    free(prtcls_c)
    free(tracks_c)

    # ===========================================
    # Return
    # ===========================================
    return track_py, num_steps