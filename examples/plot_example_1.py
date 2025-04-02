"""Example: The optimization of an idealized muon shield based on
   two conical sections. Here we take care that we are not exceeding
   the limit in z (29.7 m)
   ==========
"""

import os
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import datetime

from scipy.optimize import minimize

import trackship as ts

# %%
# Script parameters
# =================

# the date and time
date_time_str = '{date:%Y-%m-%d_%H-%M-%S}'.format(date=datetime.datetime.now())

# the log filename
LOG_FN = 'opt_log_{}.txt'.format(date_time_str)

with open(LOG_FN, 'w') as f:
    f.write('x1,x2,x3,x4,x5,x6,s,a,c,f\n')

# what is 1 GEV
GEV = 5.344286e-19

# arguments for the objective function
obj_args = {
                "x_track_min(m)": -4.0,
                "x_track_max(m)": 4.0,
                "z_track_min(m)": -1.0,
                "z_track_max(m)": 35.0,
                "max_steps": 50000,
                "subsampling_rate": 5,
                "w_scoring_plane(m)": 2.0,
                "x_lim(m)": 3.0,
                "log_filename": LOG_FN
            }

# arguments for the generator
gen_args = {
                "l_mhs(m)": 2.3,
                "w_mhs(m)": 0.3,
                "v_mhs(m)": 0.3,
                "b_mhs(T)": 1.9,
                "b_s1(T)": 1.8,
                "b_s2(T)": 1.7,
                "l_tot(m)": 29.7,
                "l_1(m)": 1.9*17.0/1.8,
                "dz(m)": 0.2
}

# the optimization method
METHOD = 'Powell'
# METHOD = 'Nelder-Mead'

plt.ion()

fig = plt.figure()
ax = fig.add_subplot(111)

# %%
# Make the particles
# ==================
NUM_P = 100
NUM_PHI = 5
PHI_MAX = 1*np.pi/180.0

NUM_MUONS = NUM_P*NUM_PHI

pz_0 = np.logspace(0.0, np.log10(350), NUM_P)*GEV
phi_0 = np.linspace(0.00, PHI_MAX, NUM_PHI)

particles = []

all_phi = np.zeros((NUM_P*NUM_PHI, ))
all_p = np.zeros((NUM_P*NUM_PHI, ))

for i, pphi in enumerate(phi_0):
    for j, pp in enumerate(pz_0):

        particles.append(ts.Particle(105.66*1.79e-30,
                        1.0,
                        np.array([0., 0., 0., np.sin(pphi)*pp, 0., np.cos(pphi)*pp])))

        all_phi[i*NUM_P+j] = pphi
        all_p[i*NUM_P+j] = pp

# fig = plt.figure()
# ax = fig.add_subplot(111)
# ax.hist2d(all_p/GEV, np.sin(all_phi)*all_p/GEV, norm=mpl.colors.LogNorm(), bins=100)
# ax.set_xlabel('$|P|$ in GEV')
# ax.set_ylabel('$P_t$ in GEV')
# plt.show()

# %%
# Launch optimization
# ================

x0 = np.array([0.3, 0.3, 0.3, 0.3, 0.3, 0.3])

bounds = [(0.1, 2.0), (0.1, 2.0),
          (0.1, 1.0), (0.1, 2.0),
          (0.3, 2.0), (0.3, 3.0)]

res = minimize(ts.obj_fcn_1, x0=x0, bounds=bounds,
               method=METHOD, args=(ts.type_1_generator,
                                    particles,
                                    obj_args,
                                    gen_args,
                                    ax))

print(res)

x_opt = res.x

np.save('x_opt_{}'.format(date_time_str), x_opt)