import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

from .configuration import Configuration

def type_1_generator(x, args):
    '''This is a 3 magnet generator.
    The magnets in section 1 and 2 are considered as
    single magnets.
    The generator has only 6 design variables. 
    
    :param x:
        The design variables:
        x[0]: The core width at the exit of section 1.
        x[1]: The void width at the exit of section 1.
        x[2]: The core width at the entrance of section 2.
        x[3]: The void width at the entrance of section 2.
        x[4]: The core width at the exit of section 2.
        x[5]: The void width at the exit of section 2.  

    :param args:
        A dictionary with the following parameters:
    l_mhs(m):
        The length of the magnetized hadron stopper in m.
    w_mhs(m):
        The core width of the magnetized hadron stopper in m.
    v_mhs(m):
        The void width of the magnetized hadron stopper in m.
    b_mhs(T):
        The B field in the magnetized hadon stopper in T.
    b_s1(T):
        The B field in section 1 in T.
    b_s2(T):
        The B field in section 2 in T.
    l_1(m):
        The length of section 1 in m.
    dz(m):
        The distance between section 1 and 2 in m.
    l_tot(m):
        The total length of the muon shieldin m.

    :return:
        The configuration given the input parameters.
    '''

    l_mhs = args["l_mhs(m)"]
    w_mhs = args["w_mhs(m)"]
    v_mhs = args["v_mhs(m)"]
    b_mhs = args["b_mhs(T)"]
    l_tot = args["l_tot(m)"]
    b_s1 = args["b_s1(T)"]
    b_s2 = args["b_s2(T)"]
    l_1 = args["l_1(m)"]
    dz = args["dz(m)"]

    # the length of section 2
    l_2 = l_tot - l_1 - l_mhs -2*dz

    # initialize the configuration
    config = Configuration()

    # the magnetized hadron stopper
    config.add_magnet(np.array([[0., 0., 0.],
                            [w_mhs, 0., 0.],
                            [w_mhs, 0., l_mhs],
                            [0., 0., l_mhs]]),
                  np.array([[v_mhs + w_mhs, 0., 0.],
                            [v_mhs + 2*w_mhs, 0., 0.],
                            [v_mhs + 2*w_mhs, 0., l_mhs],
                            [v_mhs + w_mhs, 0., l_mhs]]),
                            b_mhs, 1)

    # section 1
    z_pos = l_mhs + dz

    config.add_magnet(np.array([[0., 0., z_pos],
                            [w_mhs, 0., z_pos],
                            [x[0], 0., z_pos + l_1],
                            [0., 0., z_pos + l_1]]),
                  np.array([[v_mhs + w_mhs, 0., z_pos],
                            [v_mhs + 2*w_mhs, 0., z_pos],
                            [x[1] + 2*x[0], 0., z_pos + l_1],
                            [x[1] + x[0], 0., z_pos + l_1]]),
                            b_s1, 1)

    # section 2
    z_pos += l_1 + dz

    config.add_magnet(np.array([[0., 0., z_pos],
                            [x[2], 0., z_pos],
                            [x[4], 0., z_pos + l_2],
                            [0., 0., z_pos + l_2]]),
                  np.array([[x[2]+x[3], 0., z_pos],
                            [2*x[2]+x[3], 0., z_pos],
                            [2*x[4]+x[5], 0., z_pos + l_2],
                            [x[4]+x[5], 0., z_pos + l_2]]),
                            b_s2, 3)

    return config

def type_1_discretizer(x, num_mag_1, num_mag_2, args, plot=False):
    '''This is a 3 magnet generator.
    The magnets in section 1 and 2 are considered as
    single magnets.
    The generator has only 6 design variables. 
    
    :param x:
        The design variables:
        x[0]: The core width at the exit of section 1.
        x[1]: The void width at the exit of section 1.
        x[2]: The core width at the entrance of section 2.
        x[3]: The void width at the entrance of section 2.
        x[4]: The core width at the exit of section 2.
        x[5]: The void width at the exit of section 2.  

    :param num_mag_1:
        The number of magnets in section 1.

    :param num_mag_2_
        The number of magnets in section 2.

    :param args:
        A dictionary with the following parameters:
    l_mhs(m):
        The length of the magnetized hadron stopper in m.
    w_mhs(m):
        The core width of the magnetized hadron stopper in m.
    v_mhs(m):
        The void width of the magnetized hadron stopper in m.
    b_mhs(T):
        The B field in the magnetized hadon stopper in T.
    b_s1(T):
        The B field in section 1 in T.
    b_s2(T):
        The B field in section 2 in T.
    l_1(m):
        The length of section 1 in m.
    dz(m):
        The distance between section 1 and 2 in m.
    l_tot(m):
        The total length of the muon shieldin m.

    :return:
        The configuration given the input parameters.
    '''

    # make the configuration
    config = type_1_generator(x, args)

    # get the additional arguments for convienience
    l_mhs = args["l_mhs(m)"]
    w_mhs = args["w_mhs(m)"]
    v_mhs = args["v_mhs(m)"]
    b_mhs = args["b_mhs(T)"]
    l_tot = args["l_tot(m)"]
    b_s1 = args["b_s1(T)"]
    b_s2 = args["b_s2(T)"]
    l_1 = args["l_1(m)"]
    dz = args["dz(m)"]

    # the length of section 2
    l_2 = l_tot - l_1 - l_mhs - 2*dz

    # the core width follows the curve
    # x_c(z) = m_c.z + b_c
    m_c_1 = (x[0] - w_mhs)/l_1
    b_c_1 = w_mhs - m_c_1*(l_mhs + dz)

    m_c_2 = (x[4] - x[2])/l_2
    b_c_2 = x[2] - m_c_2*(l_mhs + l_1 + 2*dz)

    # the void follows the curve
    # x_v(z) = m_v.z + b_v
    m_v_1 = (x[1] + x[0] - v_mhs - w_mhs)/l_1
    b_v_1 = v_mhs + w_mhs - m_v_1*(l_mhs + dz)

    m_v_2 = (x[4] + x[5] - x[2] - x[3])/l_2
    b_v_2 = x[2] + x[3] - m_v_2*(l_mhs + l_1 + 2*dz)

    # we discretize this configuration
    config_disc = Configuration()

    # the magnetized hadron stopper
    config_disc.add_magnet(np.array([[0., 0., 0.],
                                [w_mhs, 0., 0.],
                                [w_mhs, 0., l_mhs],
                                [0., 0., l_mhs]]),
                    np.array([[v_mhs + w_mhs, 0., 0.],
                                [v_mhs + 2*w_mhs, 0., 0.],
                                [v_mhs + 2*w_mhs, 0., l_mhs],
                                [v_mhs + w_mhs, 0., l_mhs]]),
                                b_mhs, 1)

    l_mag_1 = (l_1 - (num_mag_1 - 1)*dz)/num_mag_1
    z_1 = np.array([l_mhs + dz,
                    l_mhs + l_mag_1 + 2*dz,
                    l_mhs + 2*l_mag_1 + 3*dz,
                    l_mhs + 3*l_mag_1 + 4*dz])

    kp_core = np.zeros((4, 3, num_mag_1))
    kp_yoke = np.zeros((4, 3, num_mag_1))

    # the total B field needs to be adjusted
    b_disc_1 = b_s1*l_1/(l_1 - 3*dz)

    print(f"The B field in section 1 is increased to {b_disc_1:.2f} T")

    for i in range(num_mag_1):
        w_core_1 = m_c_1*z_1[i] + b_c_1
        w_core_2 = m_c_1*(z_1[i] + l_mag_1) + b_c_1

        x_void_1 = m_v_1*z_1[i] + b_v_1
        x_void_2 = m_v_1*(z_1[i] + l_mag_1) + b_v_1

        kp_core[:, :, i] = np.array([[0., 0., z_1[i]],
                                    [w_core_1, 0., z_1[i]],
                                    [w_core_2, 0., z_1[i] + l_mag_1],
                                    [0., 0., z_1[i] + l_mag_1]])

        kp_yoke[:, :, i] = np.array([[x_void_1, 0., z_1[i]],
                                    [x_void_1 + w_core_1, 0., z_1[i]],
                                    [x_void_2 + w_core_2, 0., z_1[i] + l_mag_1],
                                    [x_void_2, 0., z_1[i] + l_mag_1]])

        config_disc.add_magnet(kp_core[:, :, i],
                            kp_yoke[:, :, i],
                            b_disc_1, True)
        
    l_mag_2 = (l_2 - (num_mag_2 - 1)*dz)/num_mag_2
    z_2 = np.array([l_mhs + 2*dz + l_1,
                    l_mhs + 2*dz + l_1 + l_mag_2 + dz,
                    l_mhs + 2*dz + l_1 + 2*l_mag_2 + 2*dz,
                    l_mhs + 2*dz + l_1 + 3*l_mag_2 + 3*dz])

    kp_core = np.zeros((4, 3, num_mag_2))
    kp_yoke = np.zeros((4, 3, num_mag_2))

    # the total B field needs to be adjusted
    b_disc_2 = b_s2*l_2/(l_2 - 2*dz)

    print(f"The B field in section 2 is increased to {b_disc_2:.2f} T")

    for i in range(num_mag_2):
        w_core_1 = m_c_2*z_2[i] + b_c_2
        w_core_2 = m_c_2*(z_2[i] + l_mag_2) + b_c_2

        x_void_1 = m_v_2*z_2[i] + b_v_2
        x_void_2 = m_v_2*(z_2[i] + l_mag_2) + b_v_2

        kp_core[:, :, i] = np.array([[0., 0., z_2[i]],
                                    [w_core_1, 0., z_2[i]],
                                    [w_core_2, 0., z_2[i] + l_mag_2],
                                    [0., 0., z_2[i] + l_mag_2]])

        kp_yoke[:, :, i] = np.array([[x_void_1, 0., z_2[i]],
                                    [x_void_1 + w_core_1, 0., z_2[i]],
                                    [x_void_2 + w_core_2, 0., z_2[i] + l_mag_2],
                                    [x_void_2, 0., z_2[i] + l_mag_2]])

        config_disc.add_magnet(kp_core[:, :, i],
                            kp_yoke[:, :, i],
                            b_disc_2, 3, True)

    if plot:
        # to plot the slope
        z_hr = np.linspace(0, 40, 1000)

        fig = plt.figure()
        ax = fig.add_subplot(121)
        ax.set_title('Original')
        ax.plot(m_c_1*z_hr + b_c_1, z_hr, '--', color='k')
        ax.plot(m_v_1*z_hr + b_v_1, z_hr, '--', color='k')
        ax.plot(m_c_2*z_hr + b_c_2, z_hr, '--', color='k')
        ax.plot(m_v_2*z_hr + b_v_2, z_hr, '--', color='k')
        config.plot(ax)

        ax.set_xlabel('$x$ in m')
        ax.set_ylabel('$z$ in m')
        ax = fig.add_subplot(122)
        ax.set_title('Discretized')
        ax.plot(m_c_1*z_hr + b_c_1, z_hr, '--', color='k')
        ax.plot(m_v_1*z_hr + b_v_1, z_hr, '--', color='k')
        ax.plot(m_c_2*z_hr + b_c_2, z_hr, '--', color='k')
        ax.plot(m_v_2*z_hr + b_v_2, z_hr, '--', color='k')
        config_disc.plot(ax)

        ax.set_xlabel('$x$ in m')
        ax.set_ylabel('$z$ in m')

        plt.show()

    # the number of magnets (discretized)
    num_mag_disc = len(config_disc.magnets)

    out_columns = ['x_core_1(m)', 'z_core_1(m)',
                'x_core_2(m)', 'z_core_2(m)',
                'x_core_3(m)', 'z_core_3(m)',
                'x_core_4(m)', 'z_core_4(m)',
                'x_yoke_1(m)', 'z_yoke_1(m)',
                'x_yoke_2(m)', 'z_yoke_2(m)',
                'x_yoke_3(m)', 'z_yoke_3(m)',
                'x_yoke_4(m)', 'z_yoke_4(m)',
                'B_goal(T)', 'type']

    out_data = np.zeros((num_mag_disc, 18))

    for i, m in enumerate(config_disc.magnets):
        out_data[i, :] = np.array([m.cp_core[0, 0], m.cp_core[0, 2],
                                m.cp_core[1, 0], m.cp_core[1, 2],
                                m.cp_core[2, 0], m.cp_core[2, 2],
                                m.cp_core[3, 0], m.cp_core[3, 2],
                                m.cp_yoke[0, 0], m.cp_yoke[0, 2],
                                m.cp_yoke[1, 0], m.cp_yoke[1, 2],
                                m.cp_yoke[2, 0], m.cp_yoke[2, 2],
                                m.cp_yoke[3, 0], m.cp_yoke[3, 2],
                                m.b_goal, m.type])

    return config_disc, pd.DataFrame(data=out_data, columns=out_columns)
